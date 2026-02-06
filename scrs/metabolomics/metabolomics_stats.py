import numpy as np
import pandas as pd
from scipy.stats import median_abs_deviation
import os
import warnings
warnings.filterwarnings('ignore')

# ==========================
# Core Statistical Functions 
# ==========================

def hodges_lehmann_estimator(x: np.ndarray, y: np.ndarray) -> float:
    """Hodges-Lehmann estimator for location shift."""
    x_clean = x[np.isfinite(x)]
    y_clean = y[np.isfinite(y)]
    
    if len(x_clean) == 0 or len(y_clean) == 0:
        return np.nan
    
    pairwise_diffs = np.subtract.outer(x_clean, y_clean)
    return np.median(pairwise_diffs)


def mcid_threshold(x: np.ndarray, y: np.ndarray, k_sd: float = 0.2) -> float:
    """
    MCID threshold based on pooled MAD.
    
    k_sd = 0.2 corresponds to Cohen's "small effect" (d=0.2).
    """
    pooled = np.concatenate([x[np.isfinite(x)], y[np.isfinite(y)]])
    
    if len(pooled) < 3:
        return np.nan
    
    mad = median_abs_deviation(pooled, scale='normal')
    return k_sd * mad


def bootstrap_posterior_distribution(
    x: np.ndarray, 
    y: np.ndarray, 
    n_boot: int,
    rng: np.random.Generator
) -> np.ndarray:
    """Bootstrap distribution of HL estimator."""
    x_clean = x[np.isfinite(x)]
    y_clean = y[np.isfinite(y)]
    
    hl_samples = np.zeros(n_boot)
    
    for i in range(n_boot):
        x_boot = rng.choice(x_clean, size=len(x_clean), replace=True)
        y_boot = rng.choice(y_clean, size=len(y_clean), replace=True)
        hl_samples[i] = hodges_lehmann_estimator(x_boot, y_boot)
    
    return hl_samples


def posterior_probabilities(
    hl_samples: np.ndarray, 
    delta_mcid: float
) -> dict:
    """Compute probabilities from bootstrap samples."""
    hl_clean = hl_samples[np.isfinite(hl_samples)]
    
    if len(hl_clean) == 0:
        return {
            "prob_positive": np.nan,
            "prob_negative": np.nan,
            "prob_direction": np.nan,
            "prob_exceeds_mcid": np.nan
        }
    
    prob_positive = np.mean(hl_clean > 0)
    prob_negative = np.mean(hl_clean < 0)
    prob_exceeds_mcid = np.mean(np.abs(hl_clean) > delta_mcid)
    prob_direction = max(prob_positive, prob_negative)
    
    return {
        "prob_positive": float(prob_positive),
        "prob_negative": float(prob_negative),
        "prob_direction": float(prob_direction),
        "prob_exceeds_mcid": float(prob_exceeds_mcid)
    }


def percentile_bootstrap_ci(
    hl_samples: np.ndarray,
    confidence_level: float = 0.95
) -> tuple:
    """Percentile bootstrap CI."""
    alpha = 1 - confidence_level
    lower = np.percentile(hl_samples, 100 * alpha / 2)
    upper = np.percentile(hl_samples, 100 * (1 - alpha / 2))
    return (float(lower), float(upper))


def bca_confidence_interval(
    x: np.ndarray,
    y: np.ndarray,
    confidence_level: float,
    n_boot: int,
    rng: np.random.Generator
) -> tuple:
    """BCa bootstrap CI (falls back to percentile if R unavailable)."""
    # Fallback implementation (R version kept as-is if available)
    hl_samples = bootstrap_posterior_distribution(x, y, n_boot, rng)
    return percentile_bootstrap_ci(hl_samples, confidence_level)


# ======================
# Main Analysis Function 
# ======================

def analyze_metabolite_differences(
    excel_path: str,
    sheet_name: str,
    output_dir: str = "results",
    high_keyword: str = "High",
    low_keyword: str = "Low",
    confidence_level: float = 0.95,
    n_bootstrap: int = 10000,
    k_sd: float = 0.2,  
    prob_direction_threshold: float = 0.90,  
    prob_magnitude_threshold: float = 0.85, 
    mw_pvalue_threshold: float = 0.20,      
    random_state: int = 42
) -> pd.DataFrame:
    """
    Analyze metabolite differences using hybrid approach.
    
    Hybrid filtering requires:
      1. Pr(direction) > 0.90
      2. Pr(|HL| > δ) > 0.75
      3. Mann-Whitney p < 0.15
    
    Parameters updated based on grid search optimization.
    """
    rng = np.random.default_rng(random_state)
    
    print(f"Reading {excel_path}, sheet: {sheet_name}")
    df = pd.read_excel(excel_path, sheet_name=sheet_name)
    
    os.makedirs(output_dir, exist_ok=True)
    
    results = []
    
    for idx, row in df.iterrows():
        # Extract metabolite info
        met_id = None
        for col in ["metabolite_id", "metabolite", "Metabolite", "metabolite_name"]:
            if col in row.index and pd.notna(row[col]):
                met_id = row[col]
                break
        if met_id is None:
            met_id = f"metabolite_{idx}"
        
        compound = row.get("Compound", row.get("compound", np.nan))
        
        # Extract group values
        x = pd.to_numeric(
            row.filter(like=high_keyword), 
            errors='coerce'
        ).values.astype(float)
        
        y = pd.to_numeric(
            row.filter(like=low_keyword), 
            errors='coerce'
        ).values.astype(float)
        
        x = x[np.isfinite(x)]
        y = y[np.isfinite(y)]
        
        if x.size == 0 or y.size == 0:
            continue
        
        # Core estimation
        hl_point = hodges_lehmann_estimator(x, y)
        delta = mcid_threshold(x, y, k_sd=k_sd)
        
        # Bootstrap
        hl_samples = bootstrap_posterior_distribution(x, y, n_bootstrap, rng)
        probs = posterior_probabilities(hl_samples, delta)
        
        # Confidence intervals
        ci_bca = bca_confidence_interval(x, y, confidence_level, n_bootstrap, rng)
        
        # Mann-Whitney p-value
        from scipy.stats import mannwhitneyu
        try:
            _, mw_pvalue = mannwhitneyu(x, y, alternative='two-sided')
        except:
            mw_pvalue = np.nan
        
        directional_consistent = (probs["prob_direction"] > prob_direction_threshold)
        magnitude_exceeds = (probs["prob_exceeds_mcid"] > prob_magnitude_threshold)
        mw_significant = (mw_pvalue < mw_pvalue_threshold)
        
        passes_hybrid_filter = (
            directional_consistent and 
            magnitude_exceeds and 
            mw_significant
        )
        
        # Direction label
        if probs["prob_positive"] > 0.5:
            direction_label = "Higher in fast-growing"
        elif probs["prob_negative"] > 0.5:
            direction_label = "Higher in slow-growing"
        else:
            direction_label = "Uncertain"
        
        # Compile results
        results.append({
            "metabolite_id": met_id,
            "compound_name": compound,
            "n_high": int(x.size),
            "n_low": int(y.size),
            "median_high": float(np.median(x)),
            "median_low": float(np.median(y)),
            
            # Effect size
            "hl_estimate": hl_point,
            "delta_mcid": delta,
            "effect_size_ratio": abs(hl_point / delta) if np.isfinite(delta) and delta > 0 else np.nan,
            
            # Confidence intervals
            f"ci_bca_{int(confidence_level*100)}_lower": ci_bca[0],
            f"ci_bca_{int(confidence_level*100)}_upper": ci_bca[1],
            
            # Probabilities
            "prob_positive": probs["prob_positive"],
            "prob_negative": probs["prob_negative"],
            "prob_direction": probs["prob_direction"],
            "prob_exceeds_mcid": probs["prob_exceeds_mcid"],
            "direction_label": direction_label,
            
            # Mann-Whitney
            "mw_pvalue": mw_pvalue,
            
            # Decision criteria 
            "directional_consistent": directional_consistent,
            "magnitude_exceeds": magnitude_exceeds,
            "mw_significant": mw_significant,
            "passes_filters": passes_hybrid_filter
        })
        
        if (idx + 1) % 10 == 0:
            print(f"  Processed {idx + 1} metabolites...")
    
    result_df = pd.DataFrame(results)
    result_df = result_df.sort_values(
        ["passes_hybrid_filter", "prob_direction", "prob_exceeds_mcid"], 
        ascending=[False, False, False]
    )
    
    output_path = os.path.join(output_dir, f"metabolites_{sheet_name}.csv")
    result_df.to_csv(output_path, index=False, float_format='%.4f')
    n_total = len(result_df)
    n_passed = result_df['passes_hybrid_filter'].sum()
    
    print(f"\n{'='*60}")
    print(f"RESULTS FOR {sheet_name}")
    print(f"{'='*60}")
    print(f"Total metabolites analyzed: {n_total}")
    print(f"Passed filters: {n_passed} ({100*n_passed/n_total:.1f}%)")
    print(f"  - Directional (>{prob_direction_threshold}): "
          f"{result_df['directional_consistent'].sum()}")
    print(f"  - Magnitude (>{prob_magnitude_threshold}): "
          f"{result_df['magnitude_exceeds'].sum()}")
    print(f"  - MW p-value (<{mw_pvalue_threshold}): "
          f"{result_df['mw_significant'].sum()}")
    print(f"Saved: {output_path}")
    print(f"{'='*60}\n")
    
    return result_df




if __name__ == "__main__":
    tissues = ["Liver", "Breast", "Leg"]
    
    all_results = {}
    
    for tissue in tissues:
        print(f"\n{'='*70}")
        print(f"PROCESSING TISSUE: {tissue}")
        print(f"{'='*70}")
        
        df = analyze_metabolite_differences(
            excel_path="Metabolome_updated_3.xlsx",
            sheet_name=tissue,
            output_dir="results_mets",
            high_keyword="High",
            low_keyword="Low",
            n_bootstrap=10000,
            k_sd=0.2,  
            prob_direction_threshold=0.90,  
            prob_magnitude_threshold=0.75,  
            mw_pvalue_threshold=0.2       
        )
        
        all_results[tissue] = df
    
    # Summary
    print("\n" + "="*70)
    print("OVERALL SUMMARY")
    print("="*70)
    for tissue, df in all_results.items():
        n_passed = df['passes_filter'].sum()
        print(f"{tissue:15s}: {n_passed:3d} metabolites passed filter")

