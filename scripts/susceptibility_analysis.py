"""
Cross-Species Susceptibility Assessment via Centroid-Based Distance Analysis
============================================================================

This script performs susceptibility assessment for cross-species molecular docking
data using a centroid-based ensemble approach.

Methodology:
1. Load docking metrics and identify reference species
2. Remove outliers from reference and test species using species-specific Mahalanobis distances
3. Configure confidence thresholds (hard or permissive)
4. Evaluate ensemble-level confidence scores
5. Generate PCA visualization
6. Save comprehensive results
"""

import os
os.environ["OMP_NUM_THREADS"] = '1'
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D
from adjustText import adjust_text
from scipy.spatial.distance import mahalanobis
from scipy.stats import chi2
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.covariance import EmpiricalCovariance
import warnings
warnings.filterwarnings('ignore')

# ============================================================================
# CONFIGURATION
# ============================================================================

# Metric directionality for biological interpretation
# -1: smaller values are worse (e.g., binding affinity, similarity scores)
# +1: larger values are worse (e.g., RMSD)
METRIC_DIRECTION = {
    'binding_affinity': -1,  # More negative = better binding
    'lig_rmsd': +1,          # Smaller = better fit
    'plif_tanimoto': -1,     # Higher = better similarity
    'ppsscore': -1           # Higher = better score
}

# Default confidence thresholds
DEFAULT_THRESHOLDS = {
    'plif_tanimoto': 0.5,
    'ppsscore': 0.5,
    'binding_affinity': -6.0,
    'lig_rmsd': 2.0
}

# Outlier detection threshold (chi-square percentile)
OUTLIER_THRESHOLD_PERCENTILE = 0.975


# ============================================================================
# STEP 1: DATA LOADING AND REFERENCE IDENTIFICATION
# ============================================================================

def get_reference_species(df):
    """
    Prompt user to identify the reference species.
    
    Returns:
        str: Name of reference species
    """
    print("\n" + "="*70)
    print("REFERENCE SPECIES IDENTIFICATION")
    print("="*70)
    
    print("\nAvailable species:")
    for species in sorted(df['species'].unique()):
        n = (df['species'] == species).sum()
        print(f"  - {species}")
    
    while True:
        ref_species = input("\nEnter reference species name (case-sensitive): ").strip()
        
        if ref_species in df['species'].unique():
            n_ref = (df['species'] == ref_species).sum()
            print(f"\n✓ Using {ref_species} as reference (n={n_ref} models)")
            return ref_species
        else:
            print(f"✗ Error: '{ref_species}' not found. Please try again.")


def load_and_preprocess_data(file_path, ref_species):
    """
    Load molecular docking data and standardize using REFERENCE-ONLY scaling.
    
    Args:
        file_path: Path to CSV file
        ref_species: Name of reference species
        
    Returns:
        df: Processed DataFrame
        scaler: Fitted StandardScaler
        metric_cols: List of metric column names
    """
    print("\n" + "="*70)
    print("DATA LOADING AND PREPROCESSING")
    print("="*70)
    
    # Load data
    df = pd.read_csv(file_path)
    
    # Validate required columns
    required_cols = ['binding_model', 'species', 'ensemble', 'binding_affinity', 
                     'ppsscore', 'lig_rmsd', 'plif_tanimoto']
    
    missing_cols = set(required_cols) - set(df.columns)
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    print(f"\nLoaded {len(df)} models from {df['species'].nunique()} species")
    
    # Store original metric values
    metric_cols = ['binding_affinity', 'ppsscore', 'lig_rmsd', 'plif_tanimoto']
    for col in metric_cols:
        df[f'{col}_orig'] = df[col]
    
    # Reference-only standardization
    print(f"\nStandardizing using REFERENCE-ONLY scaling ({ref_species})")
    
    ref_mask = df['species'] == ref_species
    ref_data = df.loc[ref_mask, metric_cols]
    
    print(f"  Reference samples: {ref_mask.sum()}")
    print(f"\nReference statistics (before scaling):")
    print(ref_data.describe().T[['mean', 'std', 'min', 'max']].to_string())
    
    # Fit scaler on reference only
    scaler = StandardScaler()
    scaler.fit(ref_data)
    
    # Transform all data using reference scaler
    df[metric_cols] = scaler.transform(df[metric_cols])
    
    print(f"\n✓ All data standardized using reference parameters")
    
    # Check for missing values
    if df[metric_cols].isnull().any().any():
        n_missing = df[metric_cols].isnull().any(axis=1).sum()
        print(f"\n⚠ Warning: Removing {n_missing} rows with missing values")
        df = df.dropna(subset=metric_cols)
    
    return df, scaler, metric_cols


# ============================================================================
# STEP 2: OUTLIER REMOVAL
# ============================================================================

def remove_outliers_by_species(df, metric_cols, outlier_threshold=OUTLIER_THRESHOLD_PERCENTILE):
    """
    Remove outliers from each species using species-specific Mahalanobis distances.
    
    Each species' centroid and covariance are computed independently, and outliers
    are identified relative to their own species distribution.
    
    Args:
        df: DataFrame with standardized metrics
        metric_cols: List of metric column names
        outlier_threshold: Chi-square percentile for outlier cutoff
        
    Returns:
        cleaned_indices: Index of non-outlier models
        outlier_info: List of dictionaries with outlier details
    """
    print("\n" + "="*70)
    print("OUTLIER REMOVAL (SPECIES-SPECIFIC MAHALANOBIS)")
    print("="*70)
    
    n_dims = len(metric_cols)
    threshold = np.sqrt(chi2.ppf(outlier_threshold, df=n_dims))
    
    print(f"\nOutlier detection parameters:")
    print(f"  Dimensions: {n_dims}")
    print(f"  Chi-square threshold ({outlier_threshold*100:.1f}%): {threshold:.3f}")
    
    outlier_info = []
    all_outlier_mask = pd.Series(False, index=df.index)
    
    for species in df['species'].unique():
        species_mask = df['species'] == species
        species_data = df.loc[species_mask, metric_cols]
        n_species = len(species_data)
        
        print(f"\n{species} (n={n_species}):")
        
        # Skip if insufficient data for covariance estimation
        if n_species < n_dims + 1:
            print(f"  ⚠ Skipping (need ≥{n_dims+1} samples for covariance)")
            continue
        
        # Compute species-specific centroid and covariance
        try:
            cov = EmpiricalCovariance().fit(species_data)
            species_centroid = species_data.mean(axis=0).values
            cov_inv = cov.precision_
            
            # Compute Mahalanobis distances from species centroid
            distances = np.array([
                mahalanobis(species_data.iloc[i].values, species_centroid, cov_inv)
                for i in range(n_species)
            ])
            
            # Identify outliers
            outlier_mask = distances > threshold
            n_outliers = outlier_mask.sum()
            
            print(f"  Distance range: [{distances.min():.3f}, {distances.max():.3f}]")
            print(f"  Mean ± SD: {distances.mean():.3f} ± {distances.std():.3f}")
            
            if n_outliers > 0:
                print(f"  ✗ Outliers detected: {n_outliers}")
                
                # Track outliers
                species_indices = df[species_mask].index
                for idx_local in np.where(outlier_mask)[0]:
                    idx_global = species_indices[idx_local]
                    outlier_info.append({
                        'species': species,
                        'binding_model': df.loc[idx_global, 'binding_model'],
                        'ensemble': int(df.loc[idx_global, 'ensemble']),
                        'mahalanobis_distance': distances[idx_local],
                        'threshold': threshold
                    })
                
                # Mark for removal
                all_outlier_mask.loc[species_indices[outlier_mask]] = True
            else:
                print(f"  ✓ No outliers detected")
                
        except np.linalg.LinAlgError:
            print(f"  ⚠ Covariance singular - skipping outlier detection")
            continue
    
    # Get cleaned indices
    cleaned_indices = df.index[~all_outlier_mask]
    
    print(f"\n" + "="*70)
    print(f"OUTLIER REMOVAL SUMMARY")
    print(f"="*70)
    print(f"  Total models: {len(df)}")
    print(f"  Outliers removed: {all_outlier_mask.sum()}")
    print(f"  Models retained: {len(cleaned_indices)}")
    
    return cleaned_indices, outlier_info


# ============================================================================
# STEP 3: THRESHOLD CONFIGURATION
# ============================================================================

def get_threshold_configuration():
    """
    Get confidence threshold configuration from user.
    
    Returns:
        thresholds: Dictionary of metric thresholds
        use_permissive: Boolean for permissive thresholding
        permissive_factor: Standard deviation multiplier (if permissive thresholding)
    """
    print("\n" + "="*70)
    print("CONFIDENCE THRESHOLD CONFIGURATION")
    print("="*70)
    
    # Hard thresholds
    print("\nDefault thresholds:")
    for metric, thresh in DEFAULT_THRESHOLDS.items():
        print(f"  {metric}: {thresh}")
    
    use_defaults = input("\nUse defaults? (y/n): ").strip().lower() in ['y', 'yes']
    
    if use_defaults:
        thresholds = DEFAULT_THRESHOLDS.copy()
    else:
        print("\nEnter custom thresholds (press Enter for default):")
        thresholds = {}
        for metric, default in DEFAULT_THRESHOLDS.items():
            value = input(f"  {metric} [{default}]: ").strip()
            thresholds[metric] = float(value) if value else default
    
    # Permissive thresholding option
    print("\n" + "-"*70)
    print("PERMISSIVE THRESHOLDING")
    print("-"*70)
    print("\nPermissive thresholding adjusts cutoffs based on reference variance,")
    print("making thresholds more permissive for metrics with high variability.")
    
    use_permissive = input("\nUse permissive thresholding? (y/n): ").strip().lower() in ['y', 'yes']
    
    if use_permissive:
        permissive_factor_input = input("  Tolerance (SD multiplier) [1.5]: ").strip()
        permissive_factor = float(permissive_factor_input) if permissive_factor_input else 1.5
        print(f"\n✓ Using permissive thresholds (±{permissive_factor} SD)")
    else:
        permissive_factor = None
        print("\n✓ Using hard thresholds")
    
    return thresholds, use_permissive, permissive_factor


# ============================================================================
# STEP 4: permissive THRESHOLD CALCULATION
# ============================================================================

def calculate_permissive_thresholds(thresholds, ref_data_orig, permissive_factor):
    """
    Calculate permissive thresholds based on reference variance.
    
    permissive thresholding makes cutoffs more permissive by incorporating
    reference standard deviation weighted by metric directionality.
    
    Args:
        thresholds: Dictionary of hard thresholds
        ref_data_orig: DataFrame of reference data (original scale)
        permissive_factor: Standard deviation multiplier
        
    Returns:
        effective_thresholds: Dictionary of adjusted thresholds
    """
    print("\n" + "="*70)
    print("CALCULATING PERMISSIVE THRESHOLDS")
    print("="*70)
    print(f"\nTolerance: {permissive_factor} standard deviations")
    
    effective_thresholds = {}
    
    for metric in thresholds.keys():
        ref_mean = ref_data_orig[f'{metric}_orig'].mean()
        ref_std = ref_data_orig[f'{metric}_orig'].std()
        hard_thresh = thresholds[metric]
        direction = METRIC_DIRECTION[metric]
        
        # Calculate permissive threshold
        # For "smaller is worse" (direction=-1): subtract SD to be more permissive
        # For "larger is worse" (direction=+1): add SD to be more permissive
        permissive_thresh = ref_mean - (direction * permissive_factor * ref_std)
        
        # Ensure permissive threshold is more permissive than hard threshold
        if direction == -1:  # Smaller is worse
            effective_thresholds[metric] = min(hard_thresh, permissive_thresh)
        else:  # Larger is worse
            effective_thresholds[metric] = max(hard_thresh, permissive_thresh)
        
        print(f"\n{metric}:")
        print(f"  Direction: {'larger is worse' if direction > 0 else 'smaller is worse'}")
        print(f"  Reference: {ref_mean:.3f} ± {ref_std:.3f}")
        print(f"  Hard threshold: {hard_thresh:.3f}")
        print(f"  Permissive threshold: {effective_thresholds[metric]:.3f}")
    
    return effective_thresholds


# ============================================================================
# STEP 5: CONFIDENCE SCORING
# ============================================================================

def evaluate_ensemble_confidence(df, thresholds, use_permissive, permissive_factor, 
                                ref_data_orig, cleaned_indices):
    """
    Evaluate confidence level for each ensemble based on threshold criteria.
    
    Confidence levels:
    - Strong: 3-4 metrics pass thresholds
    - Moderate: 2 metrics pass thresholds
    - Weak: 0-1 metrics pass thresholds
    
    Args:
        df: Full DataFrame
        thresholds: Dictionary of thresholds
        use_permissive: Boolean for permissive thresholding
        permissive_factor: SD multiplier (if permissive thresholding)
        ref_data_orig: Reference data in original scale
        cleaned_indices: Indices after outlier removal
        
    Returns:
        df: DataFrame with confidence scores added
    """
    print("\n" + "="*70)
    print("ENSEMBLE CONFIDENCE EVALUATION")
    print("="*70)
    
    # Determine effective thresholds
    if use_permissive:
        effective_thresholds = calculate_permissive_thresholds(
            thresholds, ref_data_orig, permissive_factor
        )
    else:
        effective_thresholds = thresholds.copy()
        print("\nUsing hard thresholds:")
        for metric, thresh in effective_thresholds.items():
            print(f"  {metric}: {thresh:.3f}")
    
    # Initialize confidence columns
    df['confidence_level'] = None
    df['n_metrics_pass'] = None
    df['plif_pass'] = None
    df['pps_pass'] = None
    df['ba_pass'] = None
    df['rmsd_pass'] = None
    
    # Evaluate each cleaned model
    print("\n" + "-"*70)
    print("EVALUATING MODELS")
    print("-"*70)
    
    df_clean = df.loc[cleaned_indices].copy()
    
    for idx in df_clean.index:
        row = df_clean.loc[idx]
        
        # Evaluate each metric
        plif_pass = row['plif_tanimoto_orig'] >= effective_thresholds['plif_tanimoto']
        pps_pass = row['ppsscore_orig'] >= effective_thresholds['ppsscore']
        ba_pass = row['binding_affinity_orig'] <= effective_thresholds['binding_affinity']
        rmsd_pass = row['lig_rmsd_orig'] <= effective_thresholds['lig_rmsd']
        
        n_pass = sum([plif_pass, pps_pass, ba_pass, rmsd_pass])
        
        # Assign confidence level
        if n_pass >= 3:
            confidence = 'Strong'
        elif n_pass == 2:
            confidence = 'Moderate'
        else:
            confidence = 'Weak'
        
        # Store results
        df.loc[idx, 'confidence_level'] = confidence
        df.loc[idx, 'n_metrics_pass'] = n_pass
        df.loc[idx, 'plif_pass'] = plif_pass
        df.loc[idx, 'pps_pass'] = pps_pass
        df.loc[idx, 'ba_pass'] = ba_pass
        df.loc[idx, 'rmsd_pass'] = rmsd_pass
    
    # Print summary by species
    print("\nConfidence distribution by species:")
    for species in df_clean['species'].unique():
        species_data = df_clean[df_clean['species'] == species]
        conf_counts = species_data['confidence_level'].value_counts()
        
        print(f"\n{species} (n={len(species_data)}):")
        for level in ['Strong', 'Moderate', 'Weak']:
            count = conf_counts.get(level, 0)
            pct = 100 * count / len(species_data) if len(species_data) > 0 else 0
            print(f"  {level}: {count} ({pct:.1f}%)")
    
    return df


def calculate_species_summary(df, cleaned_indices):
    """
    Calculate species-level confidence summary.
    
    Species confidence is determined by best ensemble performance.
    
    Args:
        df: DataFrame with ensemble confidence scores
        cleaned_indices: Indices after outlier removal
        
    Returns:
        species_summary: DataFrame with per-species statistics
    """
    print("\n" + "="*70)
    print("SPECIES-LEVEL SUMMARY")
    print("="*70)
    
    df_clean = df.loc[cleaned_indices].copy()
    
    species_summary = []
    
    for species in df_clean['species'].unique():
        species_data = df_clean[df_clean['species'] == species]
        
        # Get best ensemble
        best_idx = species_data['n_metrics_pass'].idxmax()
        best_ensemble = species_data.loc[best_idx]
        
        # Count confidence levels
        conf_counts = species_data['confidence_level'].value_counts()
        
        # Determine species-level confidence (best ensemble)
        if conf_counts.get('Strong', 0) > 0:
            species_confidence = 'Strong'
        elif conf_counts.get('Moderate', 0) > 0:
            species_confidence = 'Moderate'
        else:
            species_confidence = 'Weak'
        
        species_summary.append({
            'species': species,
            'n_ensembles': len(species_data),
            'best_ensemble': int(best_ensemble['ensemble']),
            'best_n_metrics_pass': int(best_ensemble['n_metrics_pass']),
            'species_confidence': species_confidence,
            'strong_count': int(conf_counts.get('Strong', 0)),
            'moderate_count': int(conf_counts.get('Moderate', 0)),
            'weak_count': int(conf_counts.get('Weak', 0))
        })
        
        print(f"\n{species}: {species_confidence.upper()}")
        print(f"  Best: Ensemble {best_ensemble['ensemble']} "
              f"({best_ensemble['n_metrics_pass']}/4 metrics)")
        print(f"  Distribution: Strong={conf_counts.get('Strong', 0)}, "
              f"Moderate={conf_counts.get('Moderate', 0)}, "
              f"Weak={conf_counts.get('Weak', 0)}")
    
    return pd.DataFrame(species_summary)


# ============================================================================
# STEP 6: PCA VISUALIZATION
# ============================================================================

def generate_pca_visualization(df, metric_cols, ref_species, species_summary,
                              output_dir, cleaned_indices):
    """
    Generate PCA visualization colored by species-level confidence.
    
    Args:
        df: DataFrame with all data
        metric_cols: List of metric columns
        ref_species: Name of reference species
        species_summary: DataFrame with species statistics
        output_dir: Output directory path
        cleaned_indices: Indices after outlier removal
        
    Returns:
        df: DataFrame with PCA coordinates added
    """
    print("\n" + "="*70)
    print("PCA VISUALIZATION")
    print("="*70)
    
    # Fit PCA on cleaned data, transform all data
    print("\nFitting PCA on cleaned data...")
    pca = PCA(n_components=2)
    pca.fit(df.loc[cleaned_indices, metric_cols])
    
    # Transform all data (including outliers for visualization)
    pca_coords = pca.transform(df[metric_cols])
    df['PCA1'] = pca_coords[:, 0]
    df['PCA2'] = pca_coords[:, 1]
    
    explained_var = pca.explained_variance_ratio_
    print(f"\nExplained variance:")
    print(f"  PC1: {explained_var[0]*100:.1f}%")
    print(f"  PC2: {explained_var[1]*100:.1f}%")
    print(f"  Total: {sum(explained_var)*100:.1f}%")
    
    # Prepare data subsets
    df_clean = df.loc[cleaned_indices]
    df_outliers = df.loc[~df.index.isin(cleaned_indices)]
    ref_data = df_clean[df_clean['species'] == ref_species]
    non_ref = df_clean[df_clean['species'] != ref_species]
    
    # Map species to confidence colors
    confidence_colors = {
        'Strong': '#006400',    # Dark green
        'Moderate': '#FFA500',  # Orange
        'Weak': '#8B0000'       # Dark red
    }
    
    species_conf_map = dict(zip(
        species_summary['species'],
        species_summary['species_confidence']
    ))
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Plot outliers
    if len(df_outliers) > 0:
        ax.scatter(df_outliers['PCA1'], df_outliers['PCA2'],
                  c='lightgray', s=50, alpha=0.4, marker='x',
                  label='Outliers (removed)', zorder=1)
    
    # Calculate species centroids first (needed for connecting lines)
    species_centroids = {}
    for species in non_ref['species'].unique():
        species_data = non_ref[non_ref['species'] == species]
        centroid_x = species_data['PCA1'].mean()
        centroid_y = species_data['PCA2'].mean()
        species_centroids[species] = (centroid_x, centroid_y)
    
    # Plot connecting lines from each point to its centroid
    for species in non_ref['species'].unique():
        species_data = non_ref[non_ref['species'] == species]
        centroid_x, centroid_y = species_centroids[species]
        conf_color = confidence_colors[species_conf_map[species]]
        
        for idx, row in species_data.iterrows():
            ax.plot([row['PCA1'], centroid_x], 
                   [row['PCA2'], centroid_y],
                   color=conf_color, alpha=0.15, linewidth=0.8, 
                   linestyle='-', zorder=1.5)
    
    # Plot non-reference models by confidence
    for conf_level in ['Weak', 'Moderate', 'Strong']:
        conf_data = non_ref[non_ref['confidence_level'] == conf_level]
        if len(conf_data) > 0:
            ax.scatter(conf_data['PCA1'], conf_data['PCA2'],
                      c=confidence_colors[conf_level], s=100, alpha=0.6,
                      edgecolors='black', linewidth=0.5,
                      label=f'{conf_level} Confidence', zorder=2)
    
    # Plot species centroids and prepare labels for adjustment
    texts = []
    for species in non_ref['species'].unique():
        centroid_x, centroid_y = species_centroids[species]
        conf_color = confidence_colors[species_conf_map[species]]
        
        ax.scatter(centroid_x, centroid_y,
                  c=conf_color, s=120, alpha=1.0,
                  edgecolors='black', linewidth=2, marker='D',
                  zorder=5)
        
        # Create text annotation but don't position it yet
        text = ax.text(centroid_x, centroid_y, species,
                      fontsize=10, fontweight='bold',
                      bbox=dict(boxstyle='round,pad=0.4',
                               facecolor=conf_color, alpha=0.8,
                               edgecolor='black', linewidth=1.5),
                      zorder=6, ha='center', va='center')
        texts.append(text)
    
    # Plot reference cluster with ellipse
    ref_coords = np.column_stack([ref_data['PCA1'], ref_data['PCA2']])
    ref_centroid = ref_coords.mean(axis=0)
    
    # Covariance ellipse (2 SD)
    cov_matrix = np.cov(ref_coords.T)
    eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)
    order = eigenvalues.argsort()[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    
    angle = np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]))
    width, height = 2 * 2 * np.sqrt(eigenvalues)
    
    ellipse = Ellipse(ref_centroid, width, height, angle=angle,
                     facecolor="#033B57", edgecolor="#0A5276",
                     linewidth=3, alpha=0.3, zorder=3)
    ax.add_patch(ellipse)
    
    # Plot reference points
    ax.scatter(ref_data['PCA1'], ref_data['PCA2'],
              c='#0A5276', s=120, alpha=0.8, marker='^',
              edgecolors='#033B57', linewidth=1.5,
              label=f'Reference ({ref_species})', zorder=4)
    
    # Plot reference centroid
    ax.scatter(ref_centroid[0], ref_centroid[1],
              c='black', s=120, alpha=1.0, marker='s',
              edgecolors='white', linewidth=2,
              label='Reference Centroid', zorder=7)
    
    # Adjust text positions to avoid overlaps
    # This iteratively moves labels to minimize overlaps
    adjust_text(texts, 
                arrowprops=dict(arrowstyle='-', color='gray', lw=0.5, alpha=0.5),
                expand_points=(1.5, 1.5),  # Expand bounding boxes for clearance
                expand_text=(1.2, 1.2),
                force_points=(0.3, 0.3),   # Force away from points
                force_text=(0.5, 0.5),     # Force labels away from each other
                ax=ax)
    
    # Configure plot
    ax.set_xlabel(f'PC1 ({explained_var[0]*100:.1f}%)', fontsize=14)
    ax.set_ylabel(f'PC2 ({explained_var[1]*100:.1f}%)', fontsize=14)
    
    # Add centroid marker to legend
    handles, labels = ax.get_legend_handles_labels()
    centroid_handle = Line2D([0], [0], marker='D', color='w',
                            markerfacecolor='gray', markersize=10,
                            markeredgecolor='black', markeredgewidth=2,
                            label='Species Centroid', linestyle='')
    handles.append(centroid_handle)
    
    ax.legend(handles=handles, fontsize=11, loc='best', 
             framealpha=0.95, edgecolor='black')
    ax.grid(True, alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    
    # Save figure
    output_path = os.path.join(output_dir, "pca_plot.png")
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\n✓ PCA plot saved: {output_path}")
    plt.close()
    
    return df


# ============================================================================
# STEP 7: SAVE RESULTS
# ============================================================================

def save_results(df, species_summary, outlier_info, thresholds, 
                use_permissive, permissive_factor, output_dir, cleaned_indices):
    """
    Save all analysis results to files.
    
    Outputs:
    1. per_model_results.csv - All models with confidence scores
    2. species_summary.csv - Per-species statistics
    3. outliers_removed.csv - Tracked outliers (if any)
    4. analysis_summary.txt - Human-readable summary
    
    Args:
        df: Full DataFrame with all results
        species_summary: DataFrame with species statistics
        outlier_info: List of outlier dictionaries
        thresholds: Dictionary of thresholds used
        use_permissive: Boolean for permissive thresholding
        permissive_factor: SD multiplier (if permissive)
        output_dir: Output directory
        cleaned_indices: Indices after outlier removal
    """
    print("\n" + "="*70)
    print("SAVING RESULTS")
    print("="*70)
    
    # 1. Per-model results
    model_file = os.path.join(output_dir, "per_model_results.csv")
    model_cols = ['binding_model', 'species', 'ensemble',
                  'binding_affinity_orig', 'ppsscore_orig',
                  'lig_rmsd_orig', 'plif_tanimoto_orig',
                  'confidence_level', 'n_metrics_pass',
                  'ba_pass', 'pps_pass', 'rmsd_pass', 'plif_pass',
                  'PCA1', 'PCA2']
    
    df[model_cols].to_csv(model_file, index=False, float_format='%.4f')
    print(f"\n✓ Per-model results: {model_file}")
    
    # 2. Species summary
    species_file = os.path.join(output_dir, "species_summary.csv")
    species_summary.to_csv(species_file, index=False)
    print(f"✓ Species summary: {species_file}")
    
    # 3. Outliers
    if len(outlier_info) > 0:
        outlier_file = os.path.join(output_dir, "outliers_removed.csv")
        pd.DataFrame(outlier_info).to_csv(outlier_file, index=False, float_format='%.4f')
        print(f"✓ Outliers tracked: {outlier_file} (n={len(outlier_info)})")
    else:
        print(f"✓ No outliers to save")
    
    # 4. Text summary
    summary_file = os.path.join(output_dir, "analysis_summary.txt")
    with open(summary_file, 'w', encoding='utf-8') as f:
        f.write("="*70 + "\n")
        f.write("SPECIES SUSCEPTIBILITY ASSESSMENT - ANALYSIS SUMMARY\n")
        f.write("="*70 + "\n\n")
        
        f.write("METHODOLOGY\n")
        f.write("-"*70 + "\n")
        f.write("1. Reference-only standardization\n")
        f.write("2. Species-specific outlier removal (Mahalanobis distance)\n")
        f.write("3. Ensemble-level confidence scoring\n")
        f.write("4. PCA visualization\n\n")
        
        f.write("DATA SUMMARY\n")
        f.write("-"*70 + "\n")
        f.write(f"Total models: {len(df)}\n")
        f.write(f"Outliers removed: {len(df) - len(cleaned_indices)}\n")
        f.write(f"Models analyzed: {len(cleaned_indices)}\n")
        f.write(f"Species: {df['species'].nunique()}\n\n")
        
        f.write("THRESHOLDS\n")
        f.write("-"*70 + "\n")
        if use_permissive:
            f.write(f"Type: permissive (±{permissive_factor} SD from reference)\n\n")
        else:
            f.write("Type: Hard\n\n")
        
        for metric, thresh in thresholds.items():
            f.write(f"  {metric}: {thresh:.3f}\n")
        
        f.write("\n\nCONFIDENCE CRITERIA\n")
        f.write("-"*70 + "\n")
        f.write("Strong: 3-4 metrics pass thresholds\n")
        f.write("Moderate: 2 metrics pass thresholds\n")
        f.write("Weak: 0-1 metrics pass thresholds\n\n")
        
        f.write("SPECIES RESULTS\n")
        f.write("="*70 + "\n\n")
        
        # Sort by confidence
        species_sorted = species_summary.sort_values(
            by='species_confidence',
            key=lambda x: x.map({'Strong': 0, 'Moderate': 1, 'Weak': 2})
        )
        
        for _, row in species_sorted.iterrows():
            f.write(f"{row['species']}: {row['species_confidence'].upper()}\n")
            f.write(f"  Ensembles: {row['n_ensembles']}\n")
            f.write(f"  Best: Ensemble {row['best_ensemble']} "
                   f"({row['best_n_metrics_pass']}/4 metrics)\n")
            f.write(f"  Distribution: Strong={row['strong_count']}, "
                   f"Moderate={row['moderate_count']}, "
                   f"Weak={row['weak_count']}\n\n")
    
    print(f"✓ Analysis summary: {summary_file}")


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """
    Main execution function for susceptibility assessment.
    """
    print("\n" + "="*70)
    print("CROSS-SPECIES SUSCEPTIBILITY ASSESSMENT")
    print("="*70)
    
    # Get input file
    file_path = input("\nEnter CSV summary file path: ").strip().replace('"', '')
    
    if not os.path.exists(file_path):
        print(f"\n✗ Error: File not found: {file_path}")
        return
    
    output_dir = os.path.dirname(file_path) if os.path.dirname(file_path) else '.'
    
    # STEP 1: Load data and identify reference
    df_temp = pd.read_csv(file_path)
    ref_species = get_reference_species(df_temp)
    
    df, scaler, metric_cols = load_and_preprocess_data(file_path, ref_species)
    
    # STEP 2: Remove outliers
    cleaned_indices, outlier_info = remove_outliers_by_species(df, metric_cols)
    
    # Get reference data for permissive thresholding
    ref_mask = (df['species'] == ref_species) & df.index.isin(cleaned_indices)
    ref_data_orig = df.loc[ref_mask]
    
    # STEP 3: Get threshold configuration
    thresholds, use_permissive, permissive_factor = get_threshold_configuration()
    
    # STEP 4-5: Evaluate confidence
    df = evaluate_ensemble_confidence(
        df, thresholds, use_permissive, permissive_factor, ref_data_orig, cleaned_indices
    )
    
    species_summary = calculate_species_summary(df, cleaned_indices)
    
    # STEP 6: Generate PCA visualization
    df = generate_pca_visualization(
        df, metric_cols, ref_species, species_summary, output_dir, cleaned_indices
    )
    
    # STEP 7: Save all results
    save_results(
        df, species_summary, outlier_info, thresholds,
        use_permissive, permissive_factor, output_dir, cleaned_indices
    )
    
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print(f"\nResults saved to: {output_dir}")
    print("\nGenerated files:")
    print("  - pca_plot.png")
    print("  - per_model_results.csv")
    print("  - species_summary.csv")
    if len(outlier_info) > 0:
        print("  - outliers_removed.csv")
    print("  - analysis_summary.txt")


if __name__ == "__main__":
    main()