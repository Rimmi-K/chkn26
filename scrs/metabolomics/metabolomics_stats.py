import os
import warnings

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, median_abs_deviation
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings('ignore')


def hodges_lehmann_estimator(x: np.ndarray, y: np.ndarray) -> float:
    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
    if len(x) == 0 or len(y) == 0:
        return np.nan
    return np.median(np.subtract.outer(x, y))


def mcid_threshold(x: np.ndarray, y: np.ndarray, k_sd: float = 0.2) -> float:
    """MCID based on pooled MAD. k_sd=0.2 corresponds to Cohen's small effect."""
    pooled = np.concatenate([x[np.isfinite(x)], y[np.isfinite(y)]])
    if len(pooled) < 3:
        return np.nan
    return k_sd * median_abs_deviation(pooled, scale='normal')


def bootstrap_hl_distribution(
    x: np.ndarray,
    y: np.ndarray,
    n_boot: int,
    rng: np.random.Generator,
) -> np.ndarray:
    x, y = x[np.isfinite(x)], y[np.isfinite(y)]
    return np.array([
        hodges_lehmann_estimator(
            rng.choice(x, size=len(x), replace=True),
            rng.choice(y, size=len(y), replace=True),
        )
        for _ in range(n_boot)
    ])


def posterior_probabilities(hl_samples: np.ndarray, delta_mcid: float) -> dict:
    samples = hl_samples[np.isfinite(hl_samples)]
    if len(samples) == 0:
        return {k: np.nan for k in ("prob_positive", "prob_negative", "prob_direction", "prob_exceeds_mcid")}

    prob_pos = float(np.mean(samples > 0))
    prob_neg = float(np.mean(samples < 0))
    return {
        "prob_positive": prob_pos,
        "prob_negative": prob_neg,
        "prob_direction": max(prob_pos, prob_neg),
        "prob_exceeds_mcid": float(np.mean(np.abs(samples) > delta_mcid)),
    }


def percentile_ci(hl_samples: np.ndarray, confidence_level: float = 0.95) -> tuple:
    alpha = 1 - confidence_level
    return (
        float(np.percentile(hl_samples, 100 * alpha / 2)),
        float(np.percentile(hl_samples, 100 * (1 - alpha / 2))),
    )


def _get_metabolite_id(row, idx: int) -> str:
    for col in ("metabolite_id", "metabolite", "Metabolite", "metabolite_name"):
        if col in row.index and pd.notna(row[col]):
            return row[col]
    return f"metabolite_{idx}"


def combine_csv_to_excel(
    results_dir: str,
    output_filename: str = "metabolomics_combined_results.xlsx",
) -> None:
    """
    Combine all CSV files from results directory into a single Excel file
    with each CSV as a separate sheet.

    Args:
        results_dir: Directory containing CSV files
        output_filename: Name of the output Excel file
    """
    import glob

    csv_files = glob.glob(os.path.join(results_dir, "metabolites_*.csv"))

    if not csv_files:
        print(f"No CSV files found in {results_dir}")
        return

    output_path = os.path.join(results_dir, output_filename)

    with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
        for csv_file in sorted(csv_files):
            # Extract sheet name from filename (e.g., "metabolites_Liver.csv" -> "Liver")
            sheet_name = os.path.basename(csv_file).replace("metabolites_", "").replace(".csv", "")
            df = pd.read_csv(csv_file)
            df.to_excel(writer, sheet_name=sheet_name, index=False)
            print(f"  Added sheet: {sheet_name} ({len(df)} metabolites)")

    print(f"\nCombined Excel file saved: {output_path}")


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
    random_state: int = 42,
) -> pd.DataFrame:
    """
    Analyze metabolite differences using a hybrid approach combining:
      - Hodges-Lehmann estimator for robust effect size
      - Bootstrap-based probability estimates
      - Mann-Whitney U test with FDR correction

    Results include both raw and FDR-corrected p-values.
    Users can apply custom thresholds based on:
      1. prob_direction (directional consistency)
      2. prob_exceeds_mcid (magnitude of effect)
      3. mw_pvalue_fdr (statistical significance with FDR correction)

    FDR correction uses the Benjamini-Hochberg method to control
    false discovery rate across all metabolites.
    """
    rng = np.random.default_rng(random_state)
    print(f"Reading {excel_path}, sheet: {sheet_name}")
    df = pd.read_excel(excel_path, sheet_name=sheet_name)
    os.makedirs(output_dir, exist_ok=True)

    ci_key = int(confidence_level * 100)
    records = []

    for idx, row in df.iterrows():
        met_id = _get_metabolite_id(row, idx)
        compound = row.get("Compound", row.get("compound", np.nan))

        x = pd.to_numeric(row.filter(like=high_keyword), errors='coerce').values.astype(float)
        y = pd.to_numeric(row.filter(like=low_keyword), errors='coerce').values.astype(float)
        x, y = x[np.isfinite(x)], y[np.isfinite(y)]

        if x.size == 0 or y.size == 0:
            continue

        hl_point = hodges_lehmann_estimator(x, y)
        delta = mcid_threshold(x, y, k_sd=k_sd)
        hl_samples = bootstrap_hl_distribution(x, y, n_bootstrap, rng)
        probs = posterior_probabilities(hl_samples, delta)
        ci = percentile_ci(hl_samples, confidence_level)

        try:
            _, mw_pvalue = mannwhitneyu(x, y, alternative='two-sided')
        except Exception:
            mw_pvalue = np.nan

        if probs["prob_positive"] > 0.5:
            direction_label = "Higher in fast-growing"
        elif probs["prob_negative"] > 0.5:
            direction_label = "Higher in slow-growing"
        else:
            direction_label = "Uncertain"

        records.append({
            "metabolite_id": met_id,
            "compound_name": compound,
            "n_high": int(x.size),
            "n_low": int(y.size),
            "median_high": float(np.median(x)),
            "median_low": float(np.median(y)),
            "hl_estimate": hl_point,
            "delta_mcid": delta,
            "effect_size_ratio": abs(hl_point / delta) if np.isfinite(delta) and delta > 0 else np.nan,
            f"ci_{ci_key}_lower": ci[0],
            f"ci_{ci_key}_upper": ci[1],
            "prob_positive": probs["prob_positive"],
            "prob_negative": probs["prob_negative"],
            "prob_direction": probs["prob_direction"],
            "prob_exceeds_mcid": probs["prob_exceeds_mcid"],
            "direction_label": direction_label,
            "mw_pvalue": mw_pvalue,
        })

        if (idx + 1) % 10 == 0:
            print(f"  Processed {idx + 1} metabolites...")

    result_df = pd.DataFrame(records)

    # Apply FDR correction to Mann-Whitney p-values
    valid_pvals = result_df['mw_pvalue'].notna()
    if valid_pvals.sum() > 0:
        # Apply Benjamini-Hochberg FDR correction
        pvals = result_df.loc[valid_pvals, 'mw_pvalue'].values
        _, pvals_corrected, _, _ = multipletests(pvals, alpha=mw_pvalue_threshold, method='fdr_bh')

        # Add corrected p-values to dataframe
        result_df['mw_pvalue_fdr'] = np.nan
        result_df.loc[valid_pvals, 'mw_pvalue_fdr'] = pvals_corrected
    else:
        result_df['mw_pvalue_fdr'] = np.nan

    # Sort by probability metrics
    result_df = result_df.sort_values(
        ["prob_direction", "prob_exceeds_mcid", "mw_pvalue_fdr"],
        ascending=[False, False, True],
    )

    output_path = os.path.join(output_dir, f"metabolites_{sheet_name}.csv")
    result_df.to_csv(output_path, index=False, float_format='%.4f')

    n_total = len(result_df)
    n_directional = (result_df['prob_direction'] > prob_direction_threshold).sum()
    n_magnitude = (result_df['prob_exceeds_mcid'] > prob_magnitude_threshold).sum()
    n_mw_uncorrected = (result_df['mw_pvalue'] < mw_pvalue_threshold).sum()
    n_mw_fdr = (result_df['mw_pvalue_fdr'] < mw_pvalue_threshold).sum()
    n_all_criteria = (
        (result_df['prob_direction'] > prob_direction_threshold) &
        (result_df['prob_exceeds_mcid'] > prob_magnitude_threshold) &
        (result_df['mw_pvalue_fdr'] < mw_pvalue_threshold)
    ).sum()

    print(f"\n{'='*60}")
    print(f"RESULTS FOR {sheet_name}")
    print(f"{'='*60}")
    print(f"Total metabolites analyzed: {n_total}")
    print(f"Meeting all criteria (FDR-corrected): {n_all_criteria} ({100*n_all_criteria/n_total:.1f}%)")
    print(f"  - Directional (>{prob_direction_threshold}): {n_directional}")
    print(f"  - Magnitude (>{prob_magnitude_threshold}): {n_magnitude}")
    print(f"  - MW p-value uncorrected (<{mw_pvalue_threshold}): {n_mw_uncorrected}")
    print(f"  - MW p-value FDR-corrected (<{mw_pvalue_threshold}): {n_mw_fdr}")
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
            excel_path="data/metabolomics/Metabolome.xlsx",
            sheet_name=tissue,
            output_dir="results/metabolomics/",
            n_bootstrap=10000,
            k_sd=0.2,
            prob_direction_threshold=0.90,
            prob_magnitude_threshold=0.75,
            mw_pvalue_threshold=0.2,
        )
        all_results[tissue] = df

    print("\n" + "="*70)
    print("OVERALL SUMMARY")
    print("="*70)
    for tissue, df in all_results.items():
        n_total = len(df)
        n_passed = (
            (df['prob_direction'] > 0.90) &
            (df['prob_exceeds_mcid'] > 0.75) &
            (df['mw_pvalue_fdr'] < 0.2)
        ).sum()
        print(f"{tissue:15s}: {n_passed:3d}/{n_total} metabolites meet all criteria (FDR-corrected)")

    # Combine all CSV files into a single Excel file
    print("\n" + "="*70)
    print("COMBINING RESULTS INTO EXCEL")
    print("="*70)
    combine_csv_to_excel(
        results_dir="results/metabolomics/",
        output_filename="metabolomics_all_tissues.xlsx"
    )
    print("="*70)
