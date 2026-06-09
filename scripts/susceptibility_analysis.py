"""
Cross-Species Susceptibility Assessment via Threshold-Based Ensemble Analysis
==============================================================================

This script performs susceptibility assessment for cross-species molecular docking
data using a multi-metric threshold approach with PCA visualization.

Input:  results/<ligand>_summary.csv  (produced by get_summary.py)
Output: results/susceptibility_analysis/  (pca_plot.png, per_model_results.csv,
                                          species_summary.csv, outliers_removed.csv,
                                          analysis_summary.txt)

Methodology:
1. Auto-detect project directory; load *_summary.csv from results/
2. Remove outliers from reference and test species using species-specific Mahalanobis distances
3. Configure confidence thresholds (hard or permissive)
4. Evaluate ensemble-level confidence scores
5. Generate PCA visualization
6. Save comprehensive results
"""

import os
os.environ["OMP_NUM_THREADS"] = '1'
import sys
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
from sklearn.covariance import LedoitWolf, MinCovDet
import warnings
# Suppress known benign numerical warnings from sparse covariance estimation
warnings.filterwarnings('ignore', category=UserWarning, module='sklearn')
warnings.filterwarnings('ignore', message='Degrees of freedom <= 0 for slice')

# ============================================================================
# CONFIGURATION
# ============================================================================

# Whether each metric passes when its value is ABOVE (True) or BELOW (False)
# the threshold.  This encodes both biological direction and the pass condition
# used in evaluate_ensemble_confidence.
#   True  → pass: value >= threshold  (higher is better: plif_tanimoto, ppsscore)
#   False → pass: value <= threshold  (lower is better:  binding_affinity, lig_rmsd)
METRIC_PASS_HIGHER = {
    'plif_tanimoto':    True,   # pass: plif >= threshold (higher = better similarity)
    'ppsscore':         True,   # pass: pps  >= threshold (higher = better score)
    'binding_affinity': False,  # pass: ba   <= threshold (more negative = better binding)
    'lig_rmsd':         False,  # pass: rmsd <= threshold (smaller = better fit)
}

# Default confidence thresholds (original, unstandardized scale)
DEFAULT_THRESHOLDS = {
    'plif_tanimoto': 0.5,
    'ppsscore': 0.5,
    'binding_affinity': -6.0,
    'lig_rmsd': 2.0
}

# Nominal chi-square percentile for outlier detection (Bonferroni-adjusted
# per-species inside remove_outliers_by_species)
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
        print(f"  - {species} (n={n})")

    species_map = {s.lower(): s for s in df['species'].unique()}

    while True:
        entry = input("\nEnter reference species name: ").strip()

        match = species_map.get(entry.lower())
        if match:
            n_ref = (df['species'] == match).sum()
            print(f"\n✓ Using {match} as reference (n={n_ref} models)")
            return match
        else:
            print(f"✗ Error: '{entry}' not found. Please try again.")


def load_and_preprocess_data(df_raw, ref_species):
    """
    Standardize docking data using REFERENCE-ONLY scaling.

    Accepts an already-loaded DataFrame so the CSV is not read twice.
    Missing values are removed BEFORE fitting the scaler to prevent NaNs
    in the reference subset from corrupting StandardScaler parameters.

    Args:
        df_raw: DataFrame already loaded from the CSV
        ref_species: Name of reference species

    Returns:
        df: Processed DataFrame
        scaler: Fitted StandardScaler
        metric_cols: List of metric column names
    """
    print("\n" + "="*70)
    print("DATA LOADING AND PREPROCESSING")
    print("="*70)

    df = df_raw.copy()

    # Validate required columns
    required_cols = ['binding_model', 'species', 'ensemble', 'binding_affinity',
                     'ppsscore', 'lig_rmsd', 'plif_tanimoto']

    missing_cols = set(required_cols) - set(df.columns)
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    print(f"\nLoaded {len(df)} models from {df['species'].nunique()} species")

    metric_cols = ['binding_affinity', 'ppsscore', 'lig_rmsd', 'plif_tanimoto']

    # Remove rows with missing metric values BEFORE fitting the scaler.
    # NaNs in ref_data would produce NaN mean/std, silently corrupting all
    # downstream standardized values.
    if df[metric_cols].isnull().any().any():
        n_missing = df[metric_cols].isnull().any(axis=1).sum()
        print(f"\n⚠ Warning: Removing {n_missing} rows with missing values")
        df = df.dropna(subset=metric_cols).reset_index(drop=True)

    # Store original metric values (used for threshold evaluation and permissive
    # threshold calculation, both of which operate on the unstandardized scale)
    for col in metric_cols:
        df[f'{col}_orig'] = df[col]

    # Reference-only standardization
    print(f"\nStandardizing using REFERENCE-ONLY scaling ({ref_species})")

    ref_mask = df['species'] == ref_species
    ref_data = df.loc[ref_mask, metric_cols]

    print(f"  Reference samples: {ref_mask.sum()}")
    print(f"\nReference statistics (before scaling):")
    print(ref_data.describe().T[['mean', 'std', 'min', 'max']].to_string())

    # Fit scaler on reference only, transform all species
    scaler = StandardScaler()
    scaler.fit(ref_data)
    df[metric_cols] = scaler.transform(df[metric_cols])

    print(f"\n✓ All data standardized using reference parameters")

    return df, scaler, metric_cols


# ============================================================================
# STEP 2: OUTLIER REMOVAL
# ============================================================================

def remove_outliers_by_species(df, metric_cols, outlier_threshold=OUTLIER_THRESHOLD_PERCENTILE):
    """
    Remove outliers from each species using species-specific Mahalanobis distances.

    Each species' centroid and covariance are computed independently, and outliers
    are identified relative to their own species distribution.

    Small-sample correction: the chi-square threshold is Bonferroni-adjusted
    per-species (alpha / n_species) to reduce the inflated false-positive rate
    that occurs when the chi-square approximation is applied to small n.

    Robust location: MinCovDet is used when n >= 2*(p+1) to provide a robust
    centroid estimate that resists masking (where clustered outliers pull the
    arithmetic mean toward themselves and reduce their own distances).
    LedoitWolf + sample mean is used as a fallback for smaller samples.

    Args:
        df: DataFrame with standardized metrics
        metric_cols: List of metric column names
        outlier_threshold: Nominal chi-square percentile before Bonferroni adjustment

    Returns:
        cleaned_indices: Index of non-outlier models
        outlier_info: List of dictionaries with outlier details
    """
    print("\n" + "="*70)
    print("OUTLIER REMOVAL (SPECIES-SPECIFIC MAHALANOBIS)")
    print("="*70)

    n_dims = len(metric_cols)
    alpha_nominal = 1.0 - outlier_threshold  # e.g. 0.025

    print(f"\nOutlier detection parameters:")
    print(f"  Dimensions: {n_dims}")
    print(f"  Nominal alpha: {alpha_nominal:.4f} (Bonferroni-adjusted per species)")
    print(f"  MinCovDet used when n >= {2 * (n_dims + 1)}, LedoitWolf otherwise")

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

        # Bonferroni-adjusted threshold: dividing alpha by n_species keeps the
        # per-species type-I error at alpha_nominal across all n_species tests.
        alpha_adj = alpha_nominal / n_species
        threshold = np.sqrt(chi2.ppf(1.0 - alpha_adj, df=n_dims))

        print(f"  Bonferroni threshold (√χ²): {threshold:.3f}")

        try:
            # Robust covariance and location estimation.
            # MinCovDet provides masking-resistant centroid and precision matrix.
            # LedoitWolf shrinkage stabilizes the precision matrix for small n
            # but uses the non-robust arithmetic mean as its location.
            if n_species >= 2 * (n_dims + 1):
                cov_obj = MinCovDet(random_state=42).fit(species_data)
                species_centroid = cov_obj.location_
                cov_inv = cov_obj.precision_
                estimator_name = "MinCovDet"
            else:
                cov_obj = LedoitWolf().fit(species_data)
                species_centroid = species_data.mean(axis=0).values
                cov_inv = cov_obj.precision_
                estimator_name = "LedoitWolf"

            print(f"  Covariance estimator: {estimator_name}")

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
    print("\nPermissive thresholding widens each cutoff to the reference distribution's")
    print("bad tail (mean ± factor × SD).  The effective threshold is the more lenient")
    print("of the hard cutoff and this data-driven value.")

    use_permissive = input("\nUse permissive thresholding? (y/n): ").strip().lower() in ['y', 'yes']

    if use_permissive:
        permissive_factor_input = input("  Tolerance (SD multiplier) [1.5]: ").strip()
        permissive_factor = float(permissive_factor_input) if permissive_factor_input else 1.5
        print(f"\n✓ Using permissive thresholds (factor={permissive_factor} SD)")
    else:
        permissive_factor = None
        print("\n✓ Using hard thresholds")

    return thresholds, use_permissive, permissive_factor


# ============================================================================
# STEP 4: PERMISSIVE THRESHOLD CALCULATION
# ============================================================================

def calculate_permissive_thresholds(thresholds, ref_data_orig, permissive_factor):
    """
    Calculate permissive thresholds based on reference variance.

    For each metric, the permissive threshold is placed at the reference
    distribution's bad tail: the direction that causes a model to FAIL.

    - Metrics that pass when value >= threshold (plif_tanimoto, ppsscore):
        bad tail = ref_mean - factor * ref_std  (low values are bad)
        effective = min(hard_thresh, permissive_thresh)  [lower = more lenient]

    - Metrics that pass when value <= threshold (binding_affinity, lig_rmsd):
        bad tail = ref_mean + factor * ref_std  (high values are bad)
        effective = max(hard_thresh, permissive_thresh)  [higher = more lenient]

    Args:
        thresholds: Dictionary of hard thresholds
        ref_data_orig: DataFrame of reference data (original scale)
        permissive_factor: Standard deviation multiplier

    Returns:
        effective_thresholds: Dictionary of adjusted thresholds, guaranteed
                              to be at least as lenient as the hard thresholds
    """
    print("\n" + "="*70)
    print("CALCULATING PERMISSIVE THRESHOLDS")
    print("="*70)
    print(f"\nTolerance: {permissive_factor} standard deviations")

    effective_thresholds = {}

    for metric in thresholds.keys():
        ref_mean = ref_data_orig[f'{metric}_orig'].mean()
        ref_std  = ref_data_orig[f'{metric}_orig'].std()
        hard_thresh  = thresholds[metric]
        pass_higher  = METRIC_PASS_HIGHER[metric]

        if pass_higher:
            # pass: value >= threshold  →  more lenient = lower threshold
            # bad tail is the lower end of the reference distribution
            permissive_thresh = ref_mean - permissive_factor * ref_std
            effective_thresholds[metric] = min(hard_thresh, permissive_thresh)
            direction_str = "≥"
        else:
            # pass: value <= threshold  →  more lenient = higher threshold
            # bad tail is the upper end of the reference distribution
            permissive_thresh = ref_mean + permissive_factor * ref_std
            effective_thresholds[metric] = max(hard_thresh, permissive_thresh)
            direction_str = "≤"

        print(f"\n{metric} (pass: {direction_str} threshold):")
        print(f"  Reference: {ref_mean:.3f} ± {ref_std:.3f}")
        print(f"  Hard threshold: {hard_thresh:.3f}")
        print(f"  Reference bad-tail threshold: {permissive_thresh:.3f}")
        print(f"  Effective (more lenient of the two): {effective_thresholds[metric]:.3f}")

    return effective_thresholds


# ============================================================================
# STEP 5: CONFIDENCE SCORING
# ============================================================================

def evaluate_ensemble_confidence(df, thresholds, use_permissive, permissive_factor,
                                 ref_data_orig, cleaned_indices):
    """
    Evaluate confidence level for each ensemble based on threshold criteria.

    Confidence levels:
    - Strong:   3-4 metrics pass thresholds
    - Moderate: 2 metrics pass thresholds
    - Weak:     0-1 metrics pass thresholds

    All four metrics are equally weighted.  The reference species is scored
    against the same criteria as test species.  Outlier rows (not in
    cleaned_indices) retain None for all confidence columns.

    Args:
        df: Full DataFrame
        thresholds: Dictionary of thresholds
        use_permissive: Boolean for permissive thresholding
        permissive_factor: SD multiplier (if permissive thresholding)
        ref_data_orig: Reference data in original scale (cleaned only)
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

    # Initialize confidence columns; outlier rows stay None/NaN
    df['confidence_level'] = None
    df['n_metrics_pass']   = None
    df['plif_pass']        = None
    df['pps_pass']         = None
    df['ba_pass']          = None
    df['rmsd_pass']        = None

    print("\n" + "-"*70)
    print("EVALUATING MODELS")
    print("-"*70)

    for idx in cleaned_indices:
        row = df.loc[idx]

        plif_pass = row['plif_tanimoto_orig'] >= effective_thresholds['plif_tanimoto']
        pps_pass  = row['ppsscore_orig']       >= effective_thresholds['ppsscore']
        ba_pass   = row['binding_affinity_orig'] <= effective_thresholds['binding_affinity']
        rmsd_pass = row['lig_rmsd_orig']        <= effective_thresholds['lig_rmsd']

        n_pass = sum([plif_pass, pps_pass, ba_pass, rmsd_pass])

        if n_pass >= 3:
            confidence = 'Strong'
        elif n_pass == 2:
            confidence = 'Moderate'
        else:
            confidence = 'Weak'

        df.loc[idx, 'confidence_level'] = confidence
        df.loc[idx, 'n_metrics_pass']   = n_pass
        df.loc[idx, 'plif_pass']        = plif_pass
        df.loc[idx, 'pps_pass']         = pps_pass
        df.loc[idx, 'ba_pass']          = ba_pass
        df.loc[idx, 'rmsd_pass']        = rmsd_pass

    # Print summary by species.
    # Read confidence values from df AFTER the loop above has populated them.
    # (A copy made before the loop would be stale and show all-None counts.)
    print("\nConfidence distribution by species:")
    df_clean_updated = df.loc[cleaned_indices]
    for species in sorted(df_clean_updated['species'].unique()):
        species_data = df_clean_updated[df_clean_updated['species'] == species]
        conf_counts  = species_data['confidence_level'].value_counts()
        n_sp = len(species_data)

        print(f"\n{species} (n={n_sp}):")
        for level in ['Strong', 'Moderate', 'Weak']:
            count = conf_counts.get(level, 0)
            pct   = 100 * count / n_sp if n_sp > 0 else 0
            print(f"  {level}: {count} ({pct:.1f}%)")

    return df


def calculate_species_summary(df, cleaned_indices):
    """
    Calculate species-level confidence summary.

    Species confidence uses majority vote: a level is assigned only when
    more than 50 % of that species' ensembles reach it (or better).  This
    prevents a single high-performing ensemble from inflating the species
    label when the remaining ensembles are weaker.

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
        n_total = len(species_data)

        # Best single ensemble (for reporting)
        best_idx      = species_data['n_metrics_pass'].idxmax()
        best_ensemble = species_data.loc[best_idx]

        # Count confidence levels
        conf_counts = species_data['confidence_level'].value_counts()
        n_strong   = int(conf_counts.get('Strong',   0))
        n_moderate = int(conf_counts.get('Moderate', 0))
        n_weak     = int(conf_counts.get('Weak',     0))

        # Majority-vote: >50 % of ensembles must reach the claimed level (or better)
        if n_strong / n_total > 0.5:
            species_confidence = 'Strong'
        elif (n_strong + n_moderate) / n_total > 0.5:
            species_confidence = 'Moderate'
        else:
            species_confidence = 'Weak'

        species_summary.append({
            'species':             species,
            'n_ensembles':         n_total,
            'best_ensemble':       int(best_ensemble['ensemble']),
            'best_n_metrics_pass': int(best_ensemble['n_metrics_pass']),
            'species_confidence':  species_confidence,
            'strong_count':        n_strong,
            'moderate_count':      n_moderate,
            'weak_count':          n_weak
        })

        print(f"\n{species}: {species_confidence.upper()}")
        print(f"  Best: Ensemble {best_ensemble['ensemble']} "
              f"({best_ensemble['n_metrics_pass']}/4 metrics)")
        print(f"  Distribution: Strong={n_strong}, "
              f"Moderate={n_moderate}, "
              f"Weak={n_weak}")

    return pd.DataFrame(species_summary)


# ============================================================================
# STEP 6: PCA VISUALIZATION
# ============================================================================

def generate_pca_visualization(df, metric_cols, ref_species, species_summary,
                               output_dir, cleaned_indices):
    """
    Generate PCA visualization colored by species-level confidence.

    PCA is fitted on cleaned data only; all data (including outliers) are
    then projected for visualization.

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

    # Transform all data (including outliers — extrapolated into the fitted space)
    pca_coords  = pca.transform(df[metric_cols])
    df['PCA1']  = pca_coords[:, 0]
    df['PCA2']  = pca_coords[:, 1]

    explained_var = pca.explained_variance_ratio_
    print(f"\nExplained variance:")
    print(f"  PC1: {explained_var[0]*100:.1f}%")
    print(f"  PC2: {explained_var[1]*100:.1f}%")
    print(f"  Total: {sum(explained_var)*100:.1f}%")

    # Prepare data subsets
    df_clean    = df.loc[cleaned_indices]
    df_outliers = df.loc[~df.index.isin(cleaned_indices)]
    ref_data    = df_clean[df_clean['species'] == ref_species]
    non_ref     = df_clean[df_clean['species'] != ref_species]

    # Confidence colour mapping
    confidence_colors = {
        'Strong':   '#006400',  # Dark green
        'Moderate': '#FFA500',  # Orange
        'Weak':     '#8B0000'   # Dark red
    }

    species_conf_map = dict(zip(
        species_summary['species'],
        species_summary['species_confidence']
    ))

    # ---- Figure ----
    fig, ax = plt.subplots(figsize=(12, 8))

    # Outliers
    if len(df_outliers) > 0:
        ax.scatter(df_outliers['PCA1'], df_outliers['PCA2'],
                   c='lightgray', s=50, alpha=0.4, marker='x',
                   label='Outliers (removed)', zorder=1)

    # Species centroids (needed before drawing connecting lines)
    species_centroids = {}
    for species in non_ref['species'].unique():
        sd = non_ref[non_ref['species'] == species]
        species_centroids[species] = (sd['PCA1'].mean(), sd['PCA2'].mean())

    # Connecting lines from each point to its species centroid
    for species in non_ref['species'].unique():
        sd = non_ref[non_ref['species'] == species]
        cx, cy = species_centroids[species]
        conf_color = confidence_colors[species_conf_map[species]]
        for _, row in sd.iterrows():
            ax.plot([row['PCA1'], cx], [row['PCA2'], cy],
                    color=conf_color, alpha=0.15, linewidth=0.8,
                    linestyle='-', zorder=1.5)

    # Non-reference models coloured by ensemble confidence
    for conf_level in ['Weak', 'Moderate', 'Strong']:
        conf_data = non_ref[non_ref['confidence_level'] == conf_level]
        if len(conf_data) > 0:
            ax.scatter(conf_data['PCA1'], conf_data['PCA2'],
                       c=confidence_colors[conf_level], s=100, alpha=0.6,
                       edgecolors='black', linewidth=0.5,
                       label=f'{conf_level} Confidence', zorder=2)

    # Species centroid markers and labels
    texts = []
    for species in non_ref['species'].unique():
        cx, cy    = species_centroids[species]
        conf_color = confidence_colors[species_conf_map[species]]

        ax.scatter(cx, cy, c=conf_color, s=120, alpha=1.0,
                   edgecolors='black', linewidth=2, marker='D', zorder=5)

        text = ax.text(cx, cy, species, fontsize=10, fontweight='bold',
                       bbox=dict(boxstyle='round,pad=0.4',
                                 facecolor=conf_color, alpha=0.8,
                                 edgecolor='black', linewidth=1.5),
                       zorder=6, ha='center', va='center')
        texts.append(text)

    # Reference cluster
    ref_coords  = np.column_stack([ref_data['PCA1'], ref_data['PCA2']])
    ref_centroid = ref_coords.mean(axis=0)

    # 2-SD covariance ellipse around reference cluster.
    # Requires ≥3 points and a full-rank covariance matrix.
    # eigh (symmetric solver) is used instead of eig for numerical stability;
    # eig can return complex eigenvalues for near-singular symmetric matrices.
    if len(ref_coords) >= 3:
        cov_matrix = np.cov(ref_coords.T)
        eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)
        order = eigenvalues.argsort()[::-1]
        eigenvalues  = eigenvalues[order]
        eigenvectors = eigenvectors[:, order]

        if np.all(eigenvalues > 0):
            angle        = np.degrees(np.arctan2(eigenvectors[1, 0], eigenvectors[0, 0]))
            width, height = 2 * 2 * np.sqrt(eigenvalues)
            ellipse = Ellipse(ref_centroid, width, height, angle=angle,
                              facecolor="#033B57", edgecolor="#0A5276",
                              linewidth=3, alpha=0.3, zorder=3)
            ax.add_patch(ellipse)
        else:
            print("  ⚠ Reference covariance is degenerate; 2-SD ellipse not drawn")
    else:
        print(f"  ⚠ Only {len(ref_coords)} reference points after cleaning — "
              f"ellipse requires ≥3")

    # Reference points
    ax.scatter(ref_data['PCA1'], ref_data['PCA2'],
               c='#0A5276', s=120, alpha=0.8, marker='^',
               edgecolors='#033B57', linewidth=1.5,
               label=f'Reference ({ref_species})', zorder=4)

    # Reference centroid marker
    ax.scatter(ref_centroid[0], ref_centroid[1],
               c='black', s=120, alpha=1.0, marker='s',
               edgecolors='white', linewidth=2,
               label='Reference Centroid', zorder=7)

    # Adjust text positions to avoid overlaps
    adjust_text(texts,
                arrowprops=dict(arrowstyle='-', color='gray', lw=0.5, alpha=0.5),
                expand_points=(1.5, 1.5),
                expand_text=(1.2, 1.2),
                force_points=(0.3, 0.3),
                force_text=(0.5, 0.5),
                ax=ax)

    ax.set_xlabel(f'PC1 ({explained_var[0]*100:.1f}%)', fontsize=14)
    ax.set_ylabel(f'PC2 ({explained_var[1]*100:.1f}%)', fontsize=14)

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
    1. per_model_results.csv - All models with confidence scores and is_outlier flag
    2. species_summary.csv   - Per-species statistics
    3. outliers_removed.csv  - Tracked outliers (if any)
    4. analysis_summary.txt  - Human-readable summary

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

    # Add outlier flag so downstream users can distinguish excluded rows
    # from genuinely un-scored ones
    df['is_outlier'] = ~df.index.isin(cleaned_indices)

    # 1. Per-model results (all rows; outliers have is_outlier=True and NaN confidence)
    model_file = os.path.join(output_dir, "per_model_results.csv")
    model_cols = ['binding_model', 'species', 'ensemble',
                  'binding_affinity_orig', 'ppsscore_orig',
                  'lig_rmsd_orig', 'plif_tanimoto_orig',
                  'is_outlier', 'confidence_level', 'n_metrics_pass',
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
        f.write("2. Species-specific outlier removal (Bonferroni-adjusted Mahalanobis;\n")
        f.write("   MinCovDet location when n >= 2*(p+1), else LedoitWolf + sample mean)\n")
        f.write("3. Ensemble-level confidence scoring (equal metric weighting, 3/4 = Strong)\n")
        f.write("4. Species-level confidence via majority vote (>50% of ensembles)\n")
        f.write("5. PCA visualization\n\n")
        f.write("Note: all four metrics are equally weighted.  The reference species\n")
        f.write("is scored against the same criteria as test species.\n\n")

        f.write("DATA SUMMARY\n")
        f.write("-"*70 + "\n")
        f.write(f"Total models:     {len(df)}\n")
        f.write(f"Outliers removed: {len(df) - len(cleaned_indices)}\n")
        f.write(f"Models analyzed:  {len(cleaned_indices)}\n")
        f.write(f"Species:          {df['species'].nunique()}\n\n")

        f.write("THRESHOLDS\n")
        f.write("-"*70 + "\n")
        if use_permissive:
            f.write(f"Type: Permissive (reference bad-tail; factor={permissive_factor} SD)\n\n")
        else:
            f.write("Type: Hard\n\n")

        for metric, thresh in thresholds.items():
            f.write(f"  {metric}: {thresh:.3f}\n")

        f.write("\n\nCONFIDENCE CRITERIA\n")
        f.write("-"*70 + "\n")
        f.write("Strong:   3-4 metrics pass  |  species requires >50% ensembles Strong\n")
        f.write("Moderate: 2 metrics pass     |  species requires >50% ensembles Strong+Moderate\n")
        f.write("Weak:     0-1 metrics pass   |  default\n\n")

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

    # Locate project and resolve results/ directory
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from utils import resolve_project_dir, get_project_paths

    project_dir = resolve_project_dir()
    details    = get_project_paths(project_dir)["results"]
    output_dir = os.path.join(details, "susceptibility_analysis")
    os.makedirs(output_dir, exist_ok=True)

    # Search results/ for *_summary.csv files produced by get_summary.py
    summary_csvs = sorted(
        f for f in os.listdir(details) if f.endswith("_summary.csv")
    ) if os.path.isdir(details) else []

    if not summary_csvs:
        print(f"\n✗ No *_summary.csv found in:\n  {details}")
        print("Run get_summary.py first to generate the summary file.")
        return
    elif len(summary_csvs) == 1:
        file_path = os.path.join(details, summary_csvs[0])
        print(f"\nUsing summary file: {summary_csvs[0]}")
    else:
        print(f"\nMultiple summary files found in {details}:")
        for i, name in enumerate(summary_csvs, 1):
            print(f"  [{i}] {name}")
        while True:
            choice = input("Select file number: ").strip()
            if choice.isdigit() and 1 <= int(choice) <= len(summary_csvs):
                file_path = os.path.join(details, summary_csvs[int(choice) - 1])
                break
            print("  Invalid selection, please try again.")

    # STEP 1: Load raw data once and reuse for both species selection and
    # preprocessing (avoids a redundant second pd.read_csv call)
    df_raw = pd.read_csv(file_path)
    ref_species = get_reference_species(df_raw)
    df, scaler, metric_cols = load_and_preprocess_data(df_raw, ref_species)

    # STEP 2: Remove outliers
    cleaned_indices, outlier_info = remove_outliers_by_species(df, metric_cols)

    # Reference data (cleaned only) used for permissive threshold calculation
    ref_mask      = (df['species'] == ref_species) & df.index.isin(cleaned_indices)
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
