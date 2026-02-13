"""
Cross-Species Susceptibility Assessment via Ensemble-Based Distance Analysis
========================================

This script performs susceptibility assessment for cross-species molecular docking
data using an ensemble-based approach.

Methodology:
1. Loads docking metrics from summary file (with ensemble structure)
2. Identifies reference species and fits scaler on reference only
3. Removes reference outliers via robust Mahalanobis
4. Calculates centroid-based Mahalanobis distance for all species
5. Applies multi-criteria assessment and generates confidence scores
6. Generates species-level susceptibility predictions

"""

import os
os.environ["OMP_NUM_THREADS"] = '1'
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import ListedColormap
import pandas as pd
import seaborn as sns
from scipy.spatial.distance import mahalanobis, euclidean
from scipy.stats import chi2
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.covariance import EmpiricalCovariance
import warnings
warnings.filterwarnings('ignore')


def get_ref_info(df):
    """
    Prompt user to identify the reference species/model.
    The reference species is the same as the reference model that was used
    to calculate the metrics. This is the species that we know is susceptible.
    
    This must be called BEFORE load_data to enable reference-only scaling.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe containing docking data (unscaled)
        
    Returns:
    --------
    ref_species : str
        Reference species name
    """
    # Prompt the user to identify the reference species/model
    ref_species = input("\nEnter the species/PDB name of the reference (case sensitive): ")
    
    # Check that this species exists in the dataframe
    while ref_species not in df['species'].unique():
        print(f"Error: '{ref_species}' not found in species list.")
        print("Available species:")
        for sp in df['species'].unique():
            print(f"  - {sp}")
        ref_species = input("Please enter a valid species name: ")
    
    n_ref_poses = ((df['species'] == ref_species)).sum()
    print(f"\nUsing {n_ref_poses} poses from {ref_species} as reference")
    print("\nIMPORTANT: Scaler will be fit on reference species only")
    
    return ref_species


def load_data(file_path, ref_species):
    """
    Load and preprocess molecular docking data in ensemble format.
    
    The scaler is fit ONLY on reference data. 
    This ensures that standardization reflects the reference distribution, not 
    the combined dataset. This prevents test species from influencing the metric 
    space definition.
    
    Parameters:
    -----------
    file_path : str
        Path to the CSV file containing docking metrics
        Expected columns: binding_model, species, ensemble, binding_affinity, 
                         ppsscore, lig_rmsd, plif_tanimoto
    ref_species : str
        Reference species name - scaler will be fit only on this species
        
    Returns:
    --------
    df : pd.DataFrame
        Preprocessed dataframe with standardized metrics
    scaler : StandardScaler
        Fitted scaler for the metric columns (fit on reference only)
    metric_cols : list
        Names of the metric columns
    """
    # Load data
    df = pd.read_csv(file_path)
    
    # Verify expected columns
    expected_cols = ['binding_model', 'species', 'ensemble', 'binding_affinity', 
                     'ppsscore', 'lig_rmsd', 'plif_tanimoto']
    
    if not all(col in df.columns for col in expected_cols):
        print("Warning: Expected columns not found. Found columns:")
        print(df.columns.tolist())
        print("\nExpected columns:")
        print(expected_cols)
        raise ValueError("CSV file must contain the expected columns")
    
    print(f"\nLoaded {len(df)} binding models")
    print(f"  Species: {df['species'].nunique()} ({', '.join(df['species'].unique())})")
    print(f"  Ensemble models per species: {df.groupby('species')['ensemble'].nunique().tolist()}")
    
    # Identify self-docking poses
    n_self_docking = (df['ensemble'] == 0).sum()
    print(f"  Self-docking poses: {n_self_docking}")
    
    # Store original values for reference
    df['binding_affinity_orig'] = df['binding_affinity']
    df['lig_rmsd_orig'] = df['lig_rmsd']
    
    # Invert binding affinity and ligand RMSD so higher values = better
    # This makes all metrics directionally consistent
    df["binding_affinity"] = df["binding_affinity"] * -1
    df["lig_rmsd"] = df["lig_rmsd"] * -1
    
    # Define metric columns
    metric_cols = ['binding_affinity', 'ppsscore', 'lig_rmsd', 'plif_tanimoto']
    
    # Fit scaler on reference data ONLY
    # This ensures standardization reflects reference distribution, not test species
    print("\n" + "="*70)
    print("REFERENCE-ONLY STANDARDIZATION")
    print("="*70)
    
    ref_mask = df['species'] == ref_species
    if not ref_mask.any():
        raise ValueError(f"Reference species '{ref_species}' not found in data")
    
    print(f"\nFitting scaler on reference species only: {ref_species}")
    print(f"  Reference samples: {ref_mask.sum()}")
    
    # Extract reference data for fitting
    ref_data = df.loc[ref_mask, metric_cols]
    
    print(f"\nReference distribution statistics (before scaling):")
    print(ref_data.describe().T[['mean', 'std', 'min', 'max']])
    
    # Fit scaler on reference only
    scaler = StandardScaler()
    scaler.fit(ref_data)
    
    print(f"\nScaler parameters (from reference):")
    for i, col in enumerate(metric_cols):
        print(f"  {col}: mean={scaler.mean_[i]:.3f}, std={scaler.scale_[i]:.3f}")
    
    # Transform ALL data using the reference-fitted scaler
    df[metric_cols] = scaler.transform(df[metric_cols])
    
    print(f"\nAll data transformed using reference-based scaler")
    print(f"  Note: Test species may have values outside [-3, +3] range")
    print(f"  This is expected and correct - they're scaled relative to reference")
    
    # Check for missing values
    if df[metric_cols].isnull().any().any():
        print("\nWarning: Missing values detected in metric columns")
        print(df[metric_cols].isnull().sum())
        print("Rows with missing values will be excluded from analysis")
        df = df.dropna(subset=metric_cols)
    
    # Report extreme values (not "outliers" in test species - this is expected)
    extreme_mask = (np.abs(df[metric_cols]) > 4).any(axis=1)
    if extreme_mask.sum() > 0:
        print(f"\nNote: {extreme_mask.sum()} data points with |z-score| > 4")
        print("For test species, this indicates deviation from reference distribution")
        by_species = df[extreme_mask].groupby('species').size()
        print("\nExtreme values by species:")
        for species, count in by_species.items():
            print(f"  {species}: {count}")
    
    return df, scaler, metric_cols


def correlation_analysis(df, metric_cols, output_dir):
    """
    Generate correlation analysis to check for multicollinearity.
    
    This helps determine whether to use Mahalanobis distance (accounting for
    correlations) or Euclidean distance (assuming independence).
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe containing docking metrics
    metric_cols : list
        Names of metric columns
    output_dir : str
        Directory to save output plots
        
    Returns:
    --------
    corr : pd.DataFrame
        Correlation matrix
    """
    corr = df[metric_cols].corr()
    
    # Generate correlation heatmap
    print("\nGenerating correlation plot...")
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.set_theme(style="white")
    sns.set_style("whitegrid")
    
    sns.heatmap(corr, annot=True, fmt='.2f', cmap="Greens", 
                square=True, cbar_kws={"shrink": 0.8},
                xticklabels=['Binding Affinity', 'PPS-Score', 'Ligand RMSD', 'PLIF Tc'],
                yticklabels=['Binding Affinity', 'PPS-Score', 'Ligand RMSD', 'PLIF Tc'])
    
    plt.title("Correlation Matrix of Docking Metrics", fontsize=14, pad=20)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "correlation_plot.png"), dpi=300, bbox_inches='tight')
    print(f"Correlation plot saved to: {os.path.join(output_dir, 'correlation_plot.png')}")
    plt.close()
    
    return corr


def calc_centroid_distance_with_outlier_removal(df, metric_cols, ref_species):
    """
    Calculate centroid-based Mahalanobis distance with robust outlier removal.
    =======================================
    This method identifies and removes outliers from the reference set using
    Mahalanobis distance, then calculates the centroid and covariance on the
    cleaned reference data. All models are then scored based on their distance
    from this robust centroid.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe containing all binding models (already scaled on reference)
    metric_cols : list
        Names of metric columns
    ref_species : str
        Name of reference species
        
    Returns:
    --------
    df : pd.DataFrame
        Dataframe with distance columns added
    ref_metrics_clean : pd.DataFrame
        Cleaned reference metrics (outliers removed)
    ref_centroid : np.array
        Reference centroid in metric space
    cov_matrix : np.array
        Covariance matrix (if Mahalanobis)
    distance_col : str
        Name of the distance column
    """
    # Extract reference metrics (already scaled)
    ref_metrics_df = df[df['species'] == ref_species][metric_cols].copy()
    ref_metrics_df.reset_index(drop=True, inplace=True)
    print(f"\nInitial reference set: {len(ref_metrics_df)} poses from {ref_species}")
    
    # Fit covariance on reference data
    print("\nFitting covariance on reference data...")
    cov = EmpiricalCovariance().fit(ref_metrics_df)
    cov_matrix = cov.covariance_
    cov_inv = cov.precision_
    
    # Calculate reference centroid
    ref_centroid = ref_metrics_df.mean(axis=0).values
    
    print(f"\nInitial reference centroid (standardized space):")
    for i, col in enumerate(metric_cols):
        print(f"  {col}: {ref_centroid[i]:.3f}")
    
    # Calculate Mahalanobis distance from each reference point to centroid
    ref_distances = np.array([
        mahalanobis(ref_metrics_df.iloc[i].values, ref_centroid, cov_inv)
        for i in range(len(ref_metrics_df))
    ])
    
    # Determine outlier threshold using chi-square distribution
    # For n dimensions, 97.5th percentile provides robust outlier detection
    # This corresponds to ~2 standard deviations in normal distribution
    n_dims = len(metric_cols)
    threshold = np.sqrt(chi2.ppf(0.975, df=n_dims))
    
    print(f"\nOutlier detection parameters:")
    print(f"  Dimensions: {n_dims}")
    print(f"  Chi-square threshold (97.5%): {threshold:.3f}")
    print(f"  (Squared threshold: {threshold**2:.3f})")
    
    print(f"\nReference Mahalanobis distances from centroid:")
    print(f"  Mean: {ref_distances.mean():.3f}")
    print(f"  Std: {ref_distances.std():.3f}")
    print(f"  Min: {ref_distances.min():.3f}")
    print(f"  Max: {ref_distances.max():.3f}")
    print(f"  Median: {np.median(ref_distances):.3f}")
    
    # Identify outliers
    outlier_mask = ref_distances > threshold
    n_outliers = outlier_mask.sum()
    
    if n_outliers > 0:
        print(f"\n  *** OUTLIERS DETECTED: {n_outliers} reference pose(s) ***")
        print("\n  Outlier details:")
        outlier_indices = np.where(outlier_mask)[0]
        for idx in outlier_indices:
            print(f"    Index {idx}: Mahalanobis distance = {ref_distances[idx]:.3f}")
            print(f"      Metric values: {ref_metrics_df.iloc[idx].values}")
        
        # Remove outliers
        ref_metrics_clean = ref_metrics_df[~outlier_mask].copy()
        ref_metrics_clean.reset_index(drop=True, inplace=True)
        print(f"\n  Cleaned reference set: {len(ref_metrics_clean)} poses (removed {n_outliers})")
        
        # Refit covariance on cleaned data
        print(f"\n  Refitting covariance on cleaned reference data...")
        cov = EmpiricalCovariance().fit(ref_metrics_clean)
        cov_matrix = cov.covariance_
        cov_inv = cov.precision_
        ref_centroid = ref_metrics_clean.mean(axis=0).values
        
        print(f"\n  Updated reference centroid (after outlier removal):")
        for i, col in enumerate(metric_cols):
            print(f"    {col}: {ref_centroid[i]:.3f}")
        
    else:
        print(f"\n  No outliers detected - all {len(ref_metrics_df)} reference poses retained")
        ref_metrics_clean = ref_metrics_df.copy()
    
    print(f"\nFinal reference centroid (standardized space):")
    for i, col in enumerate(metric_cols):
        print(f"  {col}: {ref_centroid[i]:.3f}")
    
    print(f"\nCovariance matrix (from cleaned reference):")
    cov_df = pd.DataFrame(cov_matrix, 
                            index=metric_cols, 
                            columns=metric_cols)
    print(cov_df.to_string(float_format='%.4f'))
    
    print(f"\nCorrelation structure in covariance:")
    corr_from_cov = np.corrcoef(ref_metrics_clean.T)
    corr_df = pd.DataFrame(corr_from_cov,
                            index=metric_cols,
                            columns=metric_cols)
    print(corr_df.to_string(float_format='%.3f'))
    
    # Calculate Mahalanobis distance from centroid for ALL points
    print(f"\nCalculating Mahalanobis distance from centroid for all {len(df)} models...")
    
    def mahalanobis_from_centroid(row):
        return mahalanobis(row[metric_cols].values, ref_centroid, cov_inv)
    
    df['mahalanobis_distance'] = df.apply(mahalanobis_from_centroid, axis=1)
    
    # Create similarity score (inverted distance for consistency)
    df['similarity_score'] = -df['mahalanobis_distance']
    distance_col = 'mahalanobis_distance'
    
    print(f"\nMahalanobis distance statistics:")
    print(f"\n  Reference species ({ref_species}):")
    ref_stats = df[df['species']==ref_species]['mahalanobis_distance']
    print(f"    Mean: {ref_stats.mean():.3f}")
    print(f"    Std: {ref_stats.std():.3f}")
    print(f"    Range: [{ref_stats.min():.3f}, {ref_stats.max():.3f}]")
    
    print(f"\n  All species combined:")
    print(f"    Mean: {df['mahalanobis_distance'].mean():.3f}")
    print(f"    Std: {df['mahalanobis_distance'].std():.3f}")
    print(f"    Range: [{df['mahalanobis_distance'].min():.3f}, {df['mahalanobis_distance'].max():.3f}]")
    
    print(f"\n  Test species (non-reference):")
    test_stats = df[df['species']!=ref_species]['mahalanobis_distance']
    print(f"    Mean: {test_stats.mean():.3f}")
    print(f"    Std: {test_stats.std():.3f}")
    print(f"    Range: [{test_stats.min():.3f}, {test_stats.max():.3f}]")
    
    # Interpretation guide
    print(f"\n  Interpretation (chi-square with {n_dims} df):")
    print(f"    Distance < {threshold:.3f}: Within 95% reference distribution")
    print(f"    Distance > {threshold:.3f}: Outside typical reference range")
    print(f"\nSimilarity score (inverted distance):")
    print(f"  Higher values = more similar to reference")
    print(f"  Mean: {df['similarity_score'].mean():.3f} ± {df['similarity_score'].std():.3f}")
    print(f"  Range: [{df['similarity_score'].min():.3f}, {df['similarity_score'].max():.3f}]")

    return df, ref_metrics_clean, ref_centroid, cov_matrix, distance_col


def multicriteria_eval(df, metric_cols, ref_metrics_df, tolerance_std=1.5):
    """
    Apply multi-criteria evaluation comparing to reference set.
    
    For each metric, we check if the test value falls within tolerance of
    the reference range (considering ALL reference poses, not just centroid).
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe containing all models
    metric_cols : list
        Names of metric columns
    ref_metrics_df : pd.DataFrame
        Dataframe with all reference pose metrics (cleaned)
    tolerance_std : float
        Tolerance in standard deviations
        
    Returns:
    --------
    df : pd.DataFrame
        Dataframe with criteria columns added
    criteria_cols : list
        Names of criteria columns
    """
    print(f"\nTolerance: {tolerance_std} standard deviations")
    print("A model passes a criterion if it falls within tolerance of reference range")
    
    criteria_cols = []
    
    for metric in metric_cols:
        # Get range from reference poses (cleaned)
        ref_min = ref_metrics_df[metric].min()
        ref_max = ref_metrics_df[metric].max()
        ref_mean = ref_metrics_df[metric].mean()
        ref_std = ref_metrics_df[metric].std() if len(ref_metrics_df) > 1 else 1.0
        
        # Define acceptable range
        lower_bound = ref_min - tolerance_std * ref_std
        upper_bound = ref_max + tolerance_std * ref_std
        
        # Check if within tolerance
        criterion_name = f"{metric}_ok"
        df[criterion_name] = ((df[metric] >= lower_bound) & 
                              (df[metric] <= upper_bound)).astype(int)
        criteria_cols.append(criterion_name)
        
        n_pass = df[criterion_name].sum()
        pct_pass = 100 * n_pass / len(df)
        
        print(f"\n{metric}:")
        print(f"  Reference mean ± std: {ref_mean:.3f} ± {ref_std:.3f}")
        print(f"  Reference range: [{ref_min:.3f}, {ref_max:.3f}]")
        print(f"  Acceptable range: [{lower_bound:.3f}, {upper_bound:.3f}]")
        print(f"  Models passing: {n_pass}/{len(df)} ({pct_pass:.1f}%)")
    
    # Calculate overall criteria score (fraction of criteria passed)
    df['criteria_score'] = df[criteria_cols].mean(axis=1)
    
    print(f"\nOverall Criteria Score:")
    print(f"  Mean: {df['criteria_score'].mean():.3f} ± {df['criteria_score'].std():.3f}")
    print(f"  Range: {df['criteria_score'].min():.3f} to {df['criteria_score'].max():.3f}")
    
    return df, criteria_cols


def calc_comb_confidence(df, dist_weight=0.5, crit_weight=0.5):
    """
    Calculate combined confidence score from distance and criteria.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with similarity_score and criteria_score
    dist_weight : float
        Weight for distance-based similarity
    crit_weight : float
        Weight for criteria-based assessment
        
    Returns:
    --------
    df : pd.DataFrame
        Dataframe with confidence_score added
    """
    print(f"\nWeights:")
    print(f"  Distance similarity: {dist_weight:.2f}")
    print(f"  Criteria assessment: {crit_weight:.2f}")
    
    # Normalize both scores to 0-1 range before combining
    # Similarity score is already inverted (higher = better)
    sim_min = df['similarity_score'].min()
    sim_max = df['similarity_score'].max()
    sim_normalized = (df['similarity_score'] - sim_min) / (sim_max - sim_min) if sim_max != sim_min else df['similarity_score']
    
    # Criteria score is already 0-1
    df['confidence_score'] = (dist_weight * sim_normalized + 
                             crit_weight * df['criteria_score'])
    
    print(f"\nConfidence Score Statistics:")
    print(f"  Mean: {df['confidence_score'].mean():.3f} ± {df['confidence_score'].std():.3f}")
    print(f"  Range: {df['confidence_score'].min():.3f} to {df['confidence_score'].max():.3f}")
    
    return df


def classify_susceptibility(df, conf_threshold=0.5, min_criteria=0.5):
    """
    Classify models as susceptible, uncertain, or not susceptible.
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with confidence_score and criteria_score
    conf_threshold : float
        Minimum confidence score for susceptibility
    min_criteria : float
        Minimum criteria score for susceptibility
        
    Returns:
    --------
    df : pd.DataFrame
        Dataframe with susceptible and uncertain columns added
    """
    print(f"\nThresholds:")
    print(f"  Confidence threshold: {conf_threshold:.2f}")
    print(f"  Minimum criteria score: {min_criteria:.2f}")
    
    # Classify based on both confidence and criteria
    df['susceptible'] = ((df['confidence_score'] >= conf_threshold) & 
                        (df['criteria_score'] >= min_criteria)).astype(int)
    
    # Mark uncertain cases (moderate confidence)
    df['uncertain'] = ((df['confidence_score'] >= 0.3) & 
                      (df['confidence_score'] < conf_threshold)).astype(int)
    
    n_susceptible = df['susceptible'].sum()
    n_uncertain = df['uncertain'].sum()
    n_not_susceptible = len(df) - n_susceptible - n_uncertain
    
    print(f"\nModel-level Classification:")
    print(f"  Susceptible: {n_susceptible} ({100*n_susceptible/len(df):.1f}%)")
    print(f"  Uncertain: {n_uncertain} ({100*n_uncertain/len(df):.1f}%)")
    print(f"  Not susceptible: {n_not_susceptible} ({100*n_not_susceptible/len(df):.1f}%)")
    
    return df


def generate_species_summary(df, ref_species):
    """
    Generate summary statistics at the species level.
    
    For each species, we aggregate across all ensemble models to get:
    - Best model (highest confidence)
    - Mean and std of confidence across all models
    - Final susceptibility call based on best model
    
    Parameters:
    -----------
    df : pd.DataFrame
        Dataframe with all models
    ref_species : str
        Reference species name
        
    Returns:
    --------
    summary : pd.DataFrame
        Species-level summary
    """
    # Group by species and get the best model (highest confidence)
    species_best = df.loc[df.groupby('species')['confidence_score'].idxmax()]
    
    # Calculate mean/std for each species
    species_stats = df.groupby('species').agg({
        'confidence_score': ['mean', 'std', 'count'],
        'similarity_score': 'mean',
        'criteria_score': 'mean'
    }).reset_index()
    
    species_stats.columns = ['species', 'mean_confidence', 'std_confidence', 'n_models',
                            'mean_similarity', 'mean_criteria']
    
    # Merge with best model info
    summary_data = []
    
    for _, best_row in species_best.iterrows():
        species = best_row['species']
        stats = species_stats[species_stats['species'] == species].iloc[0]
        
        summary_data.append({
            'species': species,
            'is_reference': species == ref_species,
            'n_models_evaluated': int(stats['n_models']),
            'confidence_score': best_row['confidence_score'],
            'mean_confidence_all_models': stats['mean_confidence'],
            'std_confidence_all_models': stats['std_confidence'],
            'similarity_score': best_row['similarity_score'],
            'criteria_score': best_row['criteria_score'],
            'susceptible': 'Yes' if best_row['susceptible'] == 1 else 'No',
            'uncertain': 'Yes' if best_row['uncertain'] == 1 else 'No',
            'binding_affinity': best_row['binding_affinity_orig'],
            'ppsscore': best_row['ppsscore'],
            'lig_rmsd': best_row['lig_rmsd_orig'],
            'plif_tanimoto': best_row['plif_tanimoto']
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df = summary_df.sort_values('confidence_score', ascending=False)
    
    print("\nSpecies Susceptibility Assessment:")
    print("="*70)
    for _, row in summary_df.iterrows():
        status = "REFERENCE" if row['is_reference'] else row['susceptible']
        print(f"\n{row['species']} - {status}")
        print(f"  Confidence score: {row['confidence_score']:.3f}")
        print(f"  Similarity to reference: {row['similarity_score']:.3f}")
        print(f"  Criteria score: {row['criteria_score']:.3f}")
        if not row['is_reference']:
            print(f"  Assessment: {row['susceptible']}" + 
                  (f" (Uncertain: {row['uncertain']})" if row['uncertain'] == 'Yes' else ""))
    
    return summary_df


def pca_vis(df, metric_cols, ref_species, ref_centroid, output_dir):
    """
    Generate PCA projection showing all models with reference centroid.
    
    Parameters:
    -----------
    df : pd.DataFrame
        All models
    metric_cols : list
        Metric column names
    ref_species : str
        Reference species name
    ref_centroid : np.array
        Reference centroid in original metric space
    output_dir : str
        Output directory
        
    Returns:
    --------
    df : pd.DataFrame
        Dataframe with PCA coordinates added
    """
    # Perform PCA
    pca = PCA(n_components=2)
    pca_coords = pca.fit_transform(df[metric_cols])
    
    df['PCA1'] = pca_coords[:, 0]
    df['PCA2'] = pca_coords[:, 1]
    
    # Transform reference centroid to PCA space
    ref_centroid_pca = pca.transform(ref_centroid.reshape(1, -1))[0]
    
    explained_var = pca.explained_variance_ratio_
    print(f"\nPCA explained variance:")
    print(f"  PC1: {explained_var[0]*100:.1f}%")
    print(f"  PC2: {explained_var[1]*100:.1f}%")
    print(f"  Total: {sum(explained_var)*100:.1f}%")
    
    # Separate reference and test data
    ref_data = df[df['species'] == ref_species]
    non_ref = df[df['species'] != ref_species]
    
    # Calculate distances from test points to reference centroid in PCA space
    distances = np.sqrt((non_ref['PCA1'] - ref_centroid_pca[0])**2 + 
                       (non_ref['PCA2'] - ref_centroid_pca[1])**2)
    
    print(f"\nDistance statistics (test points to reference centroid in PCA space):")
    print(f"  Min distance: {distances.min():.3f}")
    print(f"  Max distance: {distances.max():.3f}")
    print(f"  Mean distance: {distances.mean():.3f}")
    
    # Normalize confidence scores to size range (20 to 200)
    conf_min, conf_max = non_ref['confidence_score'].min(), non_ref['confidence_score'].max()
    sizes = 20 + (non_ref['confidence_score'] - conf_min) / (conf_max - conf_min) * 180
    
    print(f"\nConfidence score range: [{conf_min:.3f}, {conf_max:.3f}]")
    
    # Create custom colormap
    custom_cmap = LinearSegmentedColormap.from_list(
        'distance_cmap', ["#AA0743", "#0C76D2"])
    
    # Create visualization
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Plot test species with distance-based color and confidence-based size
    scatter = ax.scatter(non_ref['PCA1'], non_ref['PCA2'], 
                        c=distances,
                        cmap=custom_cmap, s=sizes, alpha=0.8, 
                        edgecolors='black', linewidth=0.5,
                        label='Test Species')
    
    # Calculate covariance for reference ellipse
    ref_coords = np.column_stack([ref_data['PCA1'], ref_data['PCA2']])
    cov_matrix = np.cov(ref_coords.T)
    
    # Eigenvalue decomposition for ellipse parameters
    eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)
    order = eigenvalues.argsort()[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    
    # Calculate ellipse angle and dimensions (2 std devs = ~95% coverage)
    angle = np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]))
    width, height = 2 * 2 * np.sqrt(eigenvalues)  # 2 std devs
    
    # Draw ellipse centered on PCA-transformed centroid
    ellipse = Ellipse(ref_centroid_pca, width, height, angle=angle,
                     facecolor='#004D40', edgecolor="#076757",
                     linewidth=2, linestyle='-', alpha=0.5, zorder=4)
    ax.add_patch(ellipse)
    
    # Plot reference species points
    ax.scatter(ref_data['PCA1'], ref_data['PCA2'],
              c='#004D40', s=120, alpha=0.8, marker='^',
              edgecolors='#00251A', linewidth=1.5,
              label=f'Reference ({ref_species})', zorder=5)
    
    # Plot reference centroid
    ax.scatter(ref_centroid_pca[0], ref_centroid_pca[1],
              c="#000000", s=120, alpha=1.0, marker='s',
              label='Reference centroid', zorder=6)
    
    # Add colorbar for distance
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('Distance from Reference Centroid', fontsize=14)
    
    # Labels and title
    ax.set_xlabel(f'PC1 ({explained_var[0]*100:.1f}%)', fontsize=14)
    ax.set_ylabel(f'PC2 ({explained_var[1]*100:.1f}%)', fontsize=14)
    ax.legend(fontsize=12, loc='best')
    ax.grid(False)
    
    plt.tight_layout()
    output_path = os.path.join(output_dir, "pca_projection_centroid.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nPCA projection saved to: {output_path}")
    plt.close()
    
    return df


def save_results(df, species_summary, criteria_cols, metric_cols, output_dir):
    """
    Save detailed results to CSV files.
    
    Parameters:
    -----------
    df : pd.DataFrame
        All models with analysis results
    species_summary : pd.DataFrame
        Species-level summary
    criteria_cols : list
        Criteria column names
    metric_cols : list
        Metric column names
    output_dir : str
        Output directory
    """
    # Determine which distance column exists
    if 'mahalanobis_distance' in df.columns:
        distance_col = 'mahalanobis_distance'
    else:
        distance_col = 'euclidean_distance'
    
    # 1. All models with full analysis
    output_cols = ['binding_model', 'species', 'ensemble', 
                  'binding_affinity_orig', 'ppsscore', 'lig_rmsd_orig', 'plif_tanimoto',
                  distance_col, 'similarity_score', 'criteria_score', 'confidence_score',
                  'susceptible', 'uncertain', 'PCA1', 'PCA2']
    
    all_models_file = os.path.join(output_dir, "all_models_assessment.csv")
    df[output_cols].to_csv(all_models_file, index=False, float_format='%.4f')
    print(f"All models results: {all_models_file}")
    
    # 2. Species summary
    species_file = os.path.join(output_dir, "species_susceptibility_summary.csv")
    species_summary.to_csv(species_file, index=False, float_format='%.4f')
    print(f"Species summary: {species_file}")
    
    # 3. Simple susceptibility list
    list_file = os.path.join(output_dir, "susceptibility_list.txt")
    with open(list_file, 'w') as f:
        f.write("SPECIES SUSCEPTIBILITY ASSESSMENT\n")
        f.write("="*70 + "\n\n")
        f.write("Methodology: Centroid-based Mahalanobis distance with outlier removal\n")
        f.write("Reference-only scaler fitting and covariance estimation\n\n")
        
        for _, row in species_summary.iterrows():
            if not row['is_reference']:
                f.write(f"{row['species']}: {row['susceptible']}\n")
                f.write(f"  Confidence: {row['confidence_score']:.3f}\n")
                f.write(f"  Uncertainty: {row['uncertain']}\n\n")
    
    print(f"Susceptibility list: {list_file}")


def main():
    """
    Main execution function for ensemble-based susceptibility assessment.
    """
    print("\n" + "="*70)
    print("ENSEMBLE-BASED SUSCEPTIBILITY ASSESSMENT")
    print("="*70)
    print("\nMethodological overview:")
    print("  1. Load data and identify reference species")
    print("  2. Fit scaler on reference species only and transform all data")
    print("  3. Correlation analysis to check metric relationships")
    print("  4. Remove outliers from reference set using Mahalanobis distance")
    print("  5. Calculate Mahalanobis distance from reference centroid for all models")
    print("  6. Multi-criteria evaluation based on reference pose distribution")
    print("  7. Combine distance and criteria into overall confidence score")
    print("  8. Classify susceptibility based on confidence and criteria")
    print("  9. Generate species-level summary and PCA visualization")
    print("  10. Save all results to output directory")
    
    # =========================================================================
    # STEP 1: Load data (preliminary - for reference identification)
    # =========================================================================
    summary_file = input("\nEnter the file path for the summary CSV file: ")
    summary_file = summary_file.replace('"', '').strip()
    
    if not os.path.exists(summary_file):
        print(f"Error: File not found: {summary_file}")
        return
    
    output_dir = os.path.dirname(summary_file) if os.path.dirname(summary_file) else '.'
    
    # Load data without scaling to identify reference
    df_temp = pd.read_csv(summary_file)
    
    # =========================================================================
    # STEP 2: Identify reference species
    # =========================================================================
    print("\n" + "="*70)
    print("REFERENCE SPECIES IDENTIFICATION")
    print("="*70)
    ref_species = get_ref_info(df_temp)

    # Reload data with reference-only scaling
    df, scaler, metric_cols = load_data(summary_file, ref_species)
    
    # =========================================================================
    # STEP 3: Correlation analysis
    # =========================================================================
    print("\n" + "="*70)
    print("CORRELATION ANALYSIS")
    print("="*70)  
    
    corr_matrix = correlation_analysis(df, metric_cols, output_dir)
    
    # =========================================================================
    # STEP 4: Calculate centroid-based distances with outlier removal
    # =========================================================================
    print("\n" + "="*70)
    print("CENTROID-BASED DISTANCE WITH OUTLIER REMOVAL")
    print("="*70)

    df, ref_metrics_clean, ref_centroid, cov_matrix, distance_col = \
        calc_centroid_distance_with_outlier_removal(
            df, metric_cols, ref_species
        )
    
    # =========================================================================
    # STEP 5: Multi-criteria assessment
    # =========================================================================    
    print("\n" + "="*70)
    print("MULTI-CRITERIA ASSESSMENT")
    print("="*70)
    
    tolerance = input("\nEnter tolerance in standard deviations [default=1.5]: ").strip()
    tolerance = float(tolerance) if tolerance else 1.5
    
    df, criteria_cols = multicriteria_eval(
        df, metric_cols, ref_metrics_clean, tolerance_std=tolerance
    )
    
    # =========================================================================
    # STEP 6: Calculate confidence scores
    # =========================================================================    
    print("\n" + "="*70)
    print("CONFIDENCE SCORE CALCULATION")
    print("="*70)

    dist_weight = input("\nEnter weight for distance similarity [default=0.5]: ").strip()
    dist_weight = float(dist_weight) if dist_weight else 0.5
    crit_weight = 1.0 - dist_weight
    
    df = calc_comb_confidence(df, dist_weight, crit_weight)
    
    # =========================================================================
    # STEP 7: Classify susceptibility (model level)
    # =========================================================================    
    print("\n" + "="*70)
    print("SUSCEPTIBILITY CLASSIFICATION")
    print("="*70)

    conf_thresh = input("\nEnter confidence threshold for susceptibility [default=0.5]: ").strip()
    conf_thresh = float(conf_thresh) if conf_thresh else 0.5
    
    min_crit = input("Enter minimum criteria score [default=0.5]: ").strip()
    min_crit = float(min_crit) if min_crit else 0.5
    
    df = classify_susceptibility(df, conf_thresh, min_crit)
    
    # =========================================================================
    # STEP 8: Generate species-level summary
    # =========================================================================    
    print("\n" + "="*70)
    print("SPECIES-LEVEL SUMMARY")
    print("="*70)
    species_summary = generate_species_summary(df, ref_species)
    
    # =========================================================================
    # STEP 9: Generate PCA visualization
    # =========================================================================    
    print("\n" + "="*70)
    print("GENERATING PCA PROJECTION")
    print("="*70)
    
    df = pca_vis(df, metric_cols, ref_species, ref_centroid, output_dir)
    
    # =========================================================================
    # STEP 10: Save results
    # =========================================================================
    print("\n" + "="*70)
    print("SAVING RESULTS")
    print("="*70)

    save_results(df, species_summary, criteria_cols, metric_cols, output_dir)
    
    # =========================================================================
    # COMPLETION
    # =========================================================================
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print(f"\nAll results have been saved to: {output_dir}")
    print("\nGenerated files:")
    print("  - correlation_plot.png")
    print("  - pca_projection_centroid.png")
    print("  - all_models_assessment.csv")
    print("  - species_susceptibility_summary.csv")
    print("  - susceptibility_list.txt")


if __name__ == "__main__":
    main()