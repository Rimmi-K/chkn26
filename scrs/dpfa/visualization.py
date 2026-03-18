import os
import logging
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import math
import matplotlib.transforms as mtransforms

def _explode_pathways(df: pd.DataFrame, pathway_column: str = "Pathways") -> pd.DataFrame:
    """
    Explode pathways column and remove duplicate pathway entries per reaction.

    After pathway merging, a reaction may have duplicate pathway names.
    We deduplicate to avoid counting reactions multiple times in DRF bar plots.

    Parameters:
    -----------
    df : pd.DataFrame
        Input dataframe
    pathway_column : str
        Name of pathway column to explode (default: "Pathways")
    """
    if pathway_column not in df.columns:
        return df

    out = df.copy()
    out[pathway_column] = out[pathway_column].fillna("").astype(str)
    out[pathway_column] = (out[pathway_column]
                           .str.replace("|", ";", regex=False)
                           .str.replace(" / ", ";", regex=False))
    out[pathway_column] = out[pathway_column].str.split(";")
    out = out.explode(pathway_column)
    out[pathway_column] = out[pathway_column].str.strip()
    out = out[out[pathway_column] != ""]

    if "reaction_id" in out.columns:
        out = out.drop_duplicates(subset=["reaction_id", pathway_column], keep="first")

    return out

def plot_regulation_counts(df: pd.DataFrame, output_dir: str, tissue: str,
                           pathways_filter=None, pathway_flux_diff=None,
                           flux_diff_threshold: float = 0.0,
                           xlim: tuple = None,
                           figsize=None):
    """
    Plot DRF histogram by pathways with variable bar width.
    Bar width depends on |total_flux_shift| (normalized by maximum).
    """
    import matplotlib.transforms as mtransforms
    from matplotlib.patches import Rectangle, FancyBboxPatch

    os.makedirs(os.path.join(output_dir, 'plots'), exist_ok=True)
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.sans-serif'] = ['Arial']
    sns.set(style="whitegrid", font_scale=0.8)

    color_dict = {
        'Increased': 'skyblue',
        'Decreased': 'lightcoral',
        'Reversed': 'plum',
        'Off': 'black',
        'On': 'limegreen',
        'Unchanged': 'gray'
    }

    # Use Subsystems for visualization if available, otherwise use Pathways
    pathway_col = "Subsystems" if "Subsystems" in df.columns else "Pathways"

    if pathway_col not in df.columns or 'DRF_category' not in df.columns:
        logging.warning(f"No {pathway_col} or DRF_category for plotting.")
        return

    df = _explode_pathways(df, pathway_column=pathway_col)

    # Apply pathway filtering if pathways_filter is provided
    if pathways_filter is not None:
        df = df[df[pathway_col].isin(pathways_filter)]

    df = df[df['DRF_category'].isin(color_dict.keys())]

    df['Pathways label'] = df[pathway_col]
    counts = df.groupby(['Pathways label', 'DRF_category']).size().unstack(fill_value=0)
    pathway_order = [p for p in pathways_filter if p in counts.index]
    counts = counts.reindex(pathway_order)

    if flux_diff_threshold > 0 and pathway_flux_diff is not None:
        flux_col_diff = "Subsystems" if "Subsystems" in pathway_flux_diff.columns else "Pathways"
        flux_diff_dict = dict(zip(
            pathway_flux_diff[flux_col_diff],
            abs(pathway_flux_diff['total_flux_shift'])
        ))

        pathway_order_before = len(pathway_order)
        pathway_order = [
            p for p in pathway_order
            if flux_diff_dict.get(p, 0) >= flux_diff_threshold
        ]
        counts = counts.reindex(pathway_order)

        logging.info(
            f"[{tissue}] Flux diff threshold={flux_diff_threshold}: "
            f"{pathway_order_before} → {len(pathway_order)} pathways"
        )

    bar_heights_pt = {}
    min_pt, max_pt = 5, 15

    if pathway_flux_diff is not None and not pathway_flux_diff.empty:
        flux_col_diff = "Subsystems" if "Subsystems" in pathway_flux_diff.columns else "Pathways"
        flux_diff_dict = dict(zip(
            pathway_flux_diff[flux_col_diff],
            abs(pathway_flux_diff['total_flux_shift'])
        ))

        relevant_diffs = [flux_diff_dict.get(p, 0) for p in pathway_order]
        max_diff = max(relevant_diffs) if relevant_diffs else 1.0

        for pathway in pathway_order:
            diff = flux_diff_dict.get(pathway, 0)
            if max_diff > 0:
                normalized = (diff / max_diff) * (max_pt - min_pt) + min_pt
            else:
                normalized = min_pt
            bar_heights_pt[pathway] = normalized
    else:
        bar_heights_pt = {p: 20 for p in pathway_order}

    # Use tissue-specific figure size or default
    if figsize and len(figsize) == 2:
        fig_size = tuple(figsize)
    else:
        fig_size = (12, 8)  # Default size for DRF histogram

    fig, ax = plt.subplots(figsize=fig_size)

    y_positions = np.arange(len(pathway_order))
    categories = [c for c in color_dict.keys() if c in counts.columns]

    fig_height_inches = fig.get_size_inches()[1]
    n_pathways = len(pathway_order)
    data_height_inches = fig_height_inches / (n_pathways + 1)
    pts_to_data = (1/72.) / data_height_inches

    left = np.zeros(len(pathway_order))

    for category in categories:
        if category not in counts.columns:
            continue

        values = counts[category].values
        heights = [bar_heights_pt[p] * pts_to_data for p in pathway_order]

        ax.barh(
            y_positions,
            values,
            height=heights,
            left=left,
            color=color_dict[category],
            label=category,
            edgecolor='white',
            linewidth=0.5
        )
        left += values

    ax.set_yticks(y_positions)
    ax.set_yticklabels(pathway_order, fontsize=8, fontfamily='Arial')
    ax.set_xlabel('Number of Reactions', fontsize=9, fontfamily='Arial')
    ax.set_ylabel('Metabolic Pathways', fontsize=9, fontfamily='Arial')
    ax.tick_params(axis='x', labelsize=8)

    if xlim is not None:
        ax.set_xlim(xlim)

    for label in ax.get_xticklabels():
        label.set_fontfamily('Arial')

    ax.invert_yaxis()

    if pathway_flux_diff is not None and not pathway_flux_diff.empty:
        flux_col_final = "Subsystems" if "Subsystems" in pathway_flux_diff.columns else "Pathways"
        flux_diff_dict = dict(zip(
            pathway_flux_diff[flux_col_final],
            pathway_flux_diff['total_flux_shift']
        ))

        bar_totals = counts.sum(axis=1)

        for i, pathway in enumerate(pathway_order):
            if pathway in flux_diff_dict:
                delta_flux = flux_diff_dict[pathway]
                bar_end_x = bar_totals[pathway]
                text = f'{delta_flux:.2f}'

                for i, pathway in enumerate(pathway_order):
                    if pathway in flux_diff_dict:
                        delta_flux = flux_diff_dict[pathway]
                        bar_end_x = bar_totals[pathway]

                        if delta_flux >= 0:
                            text = f'+{delta_flux:.2f}'
                            color = 'darkblue'
                        else:
                            text = f'{delta_flux:.2f}'
                            color = 'darkred'

                        ax.text(
                            bar_end_x + 0.8,
                            y_positions[i],
                            text,
                            va='center',
                            ha='left',
                            fontsize=8,
                            fontfamily='Arial',
                            color=color,
                            fontweight='normal',
                            bbox=dict(
                                boxstyle='round,pad=0.3',
                                facecolor='white',
                                edgecolor='none',
                                alpha=0.85
                            )
                        )

    legend = ax.legend(fontsize=8, frameon=True, edgecolor='gray',  loc='upper right', bbox_to_anchor=(1.2, 1))
    legend.get_frame().set_alpha(0.95)
    for text in legend.get_texts():
        text.set_fontfamily('Arial')

    plt.tight_layout()

    plot_path = os.path.join(output_dir, 'plots', f'drf_{tissue}.png')
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    logging.info(f"Plot saved to {plot_path}")

def create_rdpfa_summary(merged_df: pd.DataFrame, output_dir: str, tissue: str) -> pd.DataFrame:
    """
    Create unified rDPFA summary file combining DRF categories and flux analysis.

    Returns a hierarchical table with:
    - Subsystems (first column)
    - Pathways (second column)
    - DRF category counts (Increased, Decreased, Reversed, On, Off, Unchanged)
    - DRF reaction lists (which reactions belong to each category)
    - Flux metrics (flux_slow, flux_fast, total_flux_shift)

    Parameters:
    -----------
    merged_df : pd.DataFrame
        DataFrame with reaction-level DRF data
    output_dir : str
        Output directory for the summary file
    tissue : str
        Tissue name for the output filename

    Returns:
    --------
    pd.DataFrame
        Summary dataframe at pathway level
    """
    if not {'flux_slow', 'flux_fast', 'DRF_category'}.issubset(merged_df.columns):
        logging.warning(f"[{tissue}] rDPFA summary skipped (missing required columns).")
        return None

    if "Pathways" not in merged_df.columns or "Subsystems" not in merged_df.columns:
        logging.warning(f"[{tissue}] rDPFA summary skipped (missing pathway columns).")
        return None

    # Explode pathways to have one row per pathway per reaction
    df_exploded = _explode_pathways(merged_df.copy(), pathway_column="Pathways")

    # Build pathway to subsystem mapping
    pathway_to_subsystem = {}
    for _, row in merged_df.iterrows():
        if pd.notna(row.get('Pathways')) and pd.notna(row.get('Subsystems')):
            pathways = str(row['Pathways']).split(';')
            subsystems = str(row['Subsystems']).split(';')
            for pw in pathways:
                pw_clean = pw.strip()
                for sub in subsystems:
                    sub_clean = sub.strip()
                    if sub_clean not in pathway_to_subsystem.get(pw_clean, []):
                        if pw_clean not in pathway_to_subsystem:
                            pathway_to_subsystem[pw_clean] = []
                        pathway_to_subsystem[pw_clean].append(sub_clean)

    # Calculate DRF category counts per pathway
    drf_counts = df_exploded.groupby(['Pathways', 'DRF_category']).size().unstack(fill_value=0)

    # Collect reaction IDs with flux_ratio for each DRF category per pathway
    drf_categories = ['Increased', 'Decreased', 'Reversed', 'On', 'Off', 'Unchanged']

    drf_reactions = {}
    for category in drf_categories:
        category_df = df_exploded[df_exploded['DRF_category'] == category].copy()
        if len(category_df) > 0:
            # Create reaction strings with flux_ratio
            category_df['rxn_with_ratio'] = category_df.apply(
                lambda row: f"{row['reaction_id']} ({row['flux_ratio']:.2f})", axis=1
            )
            # Group by pathway and join reaction strings
            reactions_by_pathway = (
                category_df.groupby('Pathways', as_index=True)['rxn_with_ratio']
                .apply(lambda x: '; '.join(sorted(set(x))))
            )
            drf_reactions[f'{category}_reactions'] = reactions_by_pathway
        else:
            # Empty series for this category
            drf_reactions[f'{category}_reactions'] = pd.Series(dtype=str)

    drf_reactions_df = pd.DataFrame(drf_reactions)

    # Calculate flux sums and reaction counts per pathway
    flux_sums = df_exploded.groupby('Pathways')[['flux_slow', 'flux_fast']].apply(
        lambda x: x.abs().sum()
    ).reset_index()

    n_reactions = df_exploded.groupby('Pathways').size().reset_index(name='n_reactions')

    # Clean up near-zero values
    flux_sums['flux_slow'] = flux_sums['flux_slow'].where(flux_sums['flux_slow'] >= 1e-6, 0)
    flux_sums['flux_fast'] = flux_sums['flux_fast'].where(flux_sums['flux_fast'] >= 1e-6, 0)
    flux_sums['total_flux_shift'] = flux_sums['flux_fast'] - flux_sums['flux_slow']

    # Merge reaction counts
    flux_sums = flux_sums.merge(n_reactions, on='Pathways', how='left')

    # Calculate mean flux shift (normalized by number of reactions)
    flux_sums['mean_flux_shift'] = flux_sums['total_flux_shift'] / flux_sums['n_reactions']

    # Merge all data together
    summary_df = drf_counts.reset_index()
    summary_df = summary_df.merge(drf_reactions_df, left_on='Pathways', right_index=True, how='left')
    summary_df = summary_df.merge(flux_sums, on='Pathways', how='left')

    # Add subsystem information
    summary_df['Subsystems'] = summary_df['Pathways'].apply(
        lambda pw: '; '.join(pathway_to_subsystem.get(pw, ['Unknown']))
    )

    # Reorder columns: Subsystems first, then Pathways, then DRF counts, then DRF reaction lists, then flux metrics
    available_drf_cats = [cat for cat in drf_categories if cat in summary_df.columns]
    available_drf_rxn_cols = [f'{cat}_reactions' for cat in drf_categories if f'{cat}_reactions' in summary_df.columns]

    column_order = (
        ['Subsystems', 'Pathways'] +
        available_drf_cats +
        available_drf_rxn_cols +
        ['n_reactions', 'flux_fast', 'flux_slow', 'total_flux_shift', 'mean_flux_shift']
    )

    summary_df = summary_df[[col for col in column_order if col in summary_df.columns]]

    # Fill NaN in reaction lists with empty strings
    for col in available_drf_rxn_cols:
        if col in summary_df.columns:
            summary_df[col] = summary_df[col].fillna('')

    # Sort by Subsystems, then by Pathways
    summary_df = summary_df.sort_values(['Subsystems', 'Pathways'])

    # Save the pathway-level summary file
    out_path_pathways = os.path.join(output_dir, f'rdpfa_summary_pathways_{tissue}.csv')
    summary_df.to_csv(out_path_pathways, index=False)
    logging.info(f"[{tissue}] rDPFA pathway-level summary saved: {out_path_pathways}")

    # Create pathway groups level summary directly from merged_df to avoid double-counting
    groups_summary = _create_pathway_groups_summary_from_reactions(
        merged_df, drf_categories, tissue
    )

    # Save the pathway groups summary file
    out_path_subsystems = os.path.join(output_dir, f'rdpfa_summary_subsystems_{tissue}.csv')
    groups_summary.to_csv(out_path_subsystems, index=False)
    logging.info(f"[{tissue}] rDPFA subsystems-level summary saved: {out_path_subsystems}")

    return summary_df


def _create_pathway_groups_summary_from_reactions(merged_df: pd.DataFrame,
                                                   drf_categories: list,
                                                   tissue: str) -> pd.DataFrame:
    """
    Aggregate reaction-level data to pathway groups level.
    This avoids double-counting reactions that belong to multiple pathways within the same group.

    Parameters:
    -----------
    merged_df : pd.DataFrame
        DataFrame with reaction-level DRF data
    drf_categories : list
        List of DRF category names
    tissue : str
        Tissue name for logging

    Returns:
    --------
    pd.DataFrame
        Aggregated summary at pathway groups level
    """
    if "Subsystems" not in merged_df.columns or "reaction_id" not in merged_df.columns:
        logging.warning(f"[{tissue}] Cannot create pathway groups summary: missing required columns")
        return None

    # Explode by Subsystems and remove duplicates at reaction level
    df_groups = _explode_pathways(merged_df.copy(), pathway_column="Subsystems")

    # Now df_groups has unique (reaction_id, Subsystems) pairs

    # Calculate DRF category counts per pathway group
    drf_counts = df_groups.groupby(['Subsystems', 'DRF_category']).size().unstack(fill_value=0)

    # Collect reaction IDs with flux_ratio for each DRF category per pathway group
    drf_reactions = {}
    for category in drf_categories:
        category_df = df_groups[df_groups['DRF_category'] == category].copy()
        if len(category_df) > 0:
            # Create reaction strings with flux_ratio
            category_df['rxn_with_ratio'] = category_df.apply(
                lambda row: f"{row['reaction_id']} ({row['flux_ratio']:.2f})", axis=1
            )
            # Group by pathway group and join reaction strings
            reactions_by_group = (
                category_df.groupby('Subsystems', as_index=True)['rxn_with_ratio']
                .apply(lambda x: '; '.join(sorted(set(x))))
            )
            drf_reactions[f'{category}_reactions'] = reactions_by_group
        else:
            drf_reactions[f'{category}_reactions'] = pd.Series(dtype=str)

    drf_reactions_df = pd.DataFrame(drf_reactions)

    # Calculate flux sums and reaction counts per pathway group
    flux_sums = df_groups.groupby('Subsystems')[['flux_slow', 'flux_fast']].apply(
        lambda x: x.abs().sum()
    ).reset_index()

    n_reactions = df_groups.groupby('Subsystems').size().reset_index(name='n_reactions')

    # Clean up near-zero values
    flux_sums['flux_slow'] = flux_sums['flux_slow'].where(flux_sums['flux_slow'] >= 1e-6, 0)
    flux_sums['flux_fast'] = flux_sums['flux_fast'].where(flux_sums['flux_fast'] >= 1e-6, 0)
    flux_sums['total_flux_shift'] = flux_sums['flux_fast'] - flux_sums['flux_slow']

    # Merge reaction counts
    flux_sums = flux_sums.merge(n_reactions, on='Subsystems', how='left')

    # Calculate mean flux shift (normalized by number of reactions)
    flux_sums['mean_flux_shift'] = flux_sums['total_flux_shift'] / flux_sums['n_reactions']

    # Merge all data together
    summary_df = drf_counts.reset_index()
    summary_df = summary_df.merge(drf_reactions_df, left_on='Subsystems', right_index=True, how='left')
    summary_df = summary_df.merge(flux_sums, on='Subsystems', how='left')

    # Reorder columns
    available_drf_cats = [cat for cat in drf_categories if cat in summary_df.columns]
    available_drf_rxn_cols = [f'{cat}_reactions' for cat in drf_categories if f'{cat}_reactions' in summary_df.columns]

    column_order = (
        ['Subsystems'] +
        available_drf_cats +
        available_drf_rxn_cols +
        ['n_reactions', 'flux_fast', 'flux_slow', 'total_flux_shift', 'mean_flux_shift']
    )

    summary_df = summary_df[[col for col in column_order if col in summary_df.columns]]

    # Fill NaN in reaction lists with empty strings
    for col in available_drf_rxn_cols:
        if col in summary_df.columns:
            summary_df[col] = summary_df[col].fillna('')

    # Sort by Subsystems
    summary_df = summary_df.sort_values('Subsystems')

    logging.info(f"[{tissue}] Created pathway groups summary with {len(summary_df)} groups")

    return summary_df


def analyze_pathway_flux_difference(merged_df: pd.DataFrame, output_dir: str,
                                   tissue: str, threshold: float = 1.0) -> pd.DataFrame:
    """
    Analyze flux difference by pathways (sum of absolute fluxes).

    This function is kept for backward compatibility with plotting functions.
    It returns filtered pathway groups data for visualization.
    """
    if not {'flux_slow', 'flux_fast'}.issubset(merged_df.columns):
        logging.warning("Pathway flux analysis skipped (missing flux data).")
        return None

    subsystems_available = "Subsystems" in merged_df.columns

    if subsystems_available:
        df_groups = _explode_pathways(merged_df.copy(), pathway_column="Subsystems")
        n_reactions = df_groups.groupby('Subsystems').size().rename('n_reactions')
        path_sums_groups = df_groups.groupby('Subsystems')[['flux_slow', 'flux_fast']].apply(
            lambda x: x.abs().sum()
        ).reset_index()

        path_sums_groups['flux_slow'] = path_sums_groups['flux_slow'].where(
            path_sums_groups['flux_slow'] >= 1e-6, 0
        )
        path_sums_groups['flux_fast'] = path_sums_groups['flux_fast'].where(
            path_sums_groups['flux_fast'] >= 1e-6, 0
        )
        path_sums_groups = path_sums_groups.join(n_reactions, on='Subsystems')
        path_sums_groups['total_flux_shift'] = (
            path_sums_groups['flux_fast'] - path_sums_groups['flux_slow']
        )
        path_sums_groups['mean_flux_shift'] = (
            path_sums_groups['total_flux_shift'] / path_sums_groups['n_reactions']
        )

        filtered = path_sums_groups[abs(path_sums_groups['total_flux_shift']) > threshold]
        filtered = filtered.sort_values(by='total_flux_shift', ascending=False)
        return filtered
    else:
        if "Pathways" in merged_df.columns:
            df_original = _explode_pathways(merged_df.copy(), pathway_column="Pathways")
            path_sums_original = df_original.groupby('Pathways')[['flux_slow', 'flux_fast']].apply(
                lambda x: x.abs().sum()
            ).reset_index()

            path_sums_original['flux_slow'] = path_sums_original['flux_slow'].where(
                path_sums_original['flux_slow'] >= 1e-6, 0
            )
            path_sums_original['flux_fast'] = path_sums_original['flux_fast'].where(
                path_sums_original['flux_fast'] >= 1e-6, 0
            )
            path_sums_original['total_flux_shift'] = (
                path_sums_original['flux_fast'] - path_sums_original['flux_slow']
            )

            filtered = path_sums_original[abs(path_sums_original['total_flux_shift']) > threshold]
            filtered = filtered.sort_values(by='total_flux_shift', ascending=False)
            return filtered
        return None


def make_band_squash_transform(B: float = 1.0, k: float = 10.0):
    """Transform for compressing large values."""
    C = (1 - 1/k) * B
    def f(x):
        x = np.asarray(x)
        out = np.empty_like(x, dtype=float)
        m = (x >= -B) & (x <= B)
        out[m]  = x[m] / k
        out[~m] = np.sign(x[~m]) * (np.abs(x[~m]) - C)
        return out
    def finv(y):
        y = np.asarray(y)
        out = np.empty_like(y, dtype=float)
        m = (y >= -B/k) & (y <= B/k)
        out[m]  = y[m] * k
        out[~m] = np.sign(y[~m]) * (np.abs(y[~m]) + C)
        return out
    return f, finv


def filter_metabolites_by_content(heatmap_df: pd.DataFrame,
                                  min_non_empty: int = 2) -> pd.DataFrame:
    """Filter metabolites keeping only those with at least min_non_empty non-empty values."""
    non_empty_counts = heatmap_df.notna().sum(axis=0)
    keep_metabolites = non_empty_counts[non_empty_counts >= min_non_empty].index
    return heatmap_df[keep_metabolites]


def merge_identical_metabolites(heatmap_df: pd.DataFrame,
                                tolerance: int = 0) -> pd.DataFrame:
    """
    Merge metabolites with identical log2FC values using comma separation.
    """
    df_t = heatmap_df.T
    df_rounded = df_t.round(tolerance)

    def make_group_key(row):
        return tuple(x if pd.notna(x) else 'NA' for x in row)

    df_rounded['_group_key'] = df_rounded.apply(make_group_key, axis=1)

    merged_data = {}
    for group_key, group in df_rounded.groupby('_group_key'):
        if all(x == 'NA' for x in group_key):
            continue

        metabolite_names = group.index.tolist()
        combined_name = ', '.join(metabolite_names)
        merged_data[combined_name] = df_t.loc[metabolite_names[0]]

    result = pd.DataFrame(merged_data).T
    return result.T



def plot_fluxsum_log2fc_heatmap(df: pd.DataFrame, tissue: str, output_dir: str,
                                pathways_filter=None, min_metabolite_values=2,
                                merge_identical: bool = True,
                                metabolite_filter: dict = None,
                                log2fc_threshold: float = 0.0,
                                min_pathways_with_change: int = 2,
                                figsize=None):
    """
    Enhanced heatmap for mcPFA with unified pathways, metabolite filtering,
    identical metabolite merging, and Times New Roman font for compact display.

    Parameters:
    -----------
    log2fc_threshold : float
        Threshold for metabolite filtering. Default: 0.0 (no filtering by log2fc)
    min_pathways_with_change : int
        Minimum number of pathways where metabolite must have |log2FC| > threshold
        to be included. This ensures changes are reproducible. Default: 2
    """
    if df is None or df.empty:
        print(f"[plot_fluxsum_log2fc_heatmap] Empty DataFrame for {tissue} — skip.")
        return

    os.makedirs(os.path.join(output_dir, 'plots'), exist_ok=True)

    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times New Roman', 'Times', 'DejaVu Serif']

    label_fontsize = 7.5
    axis_title_fontsize = 9
    annot_fontsize = 8
    cbar_label_fontsize = 8
    cbar_tick_fontsize = 7

    pathway_column = "Subsystem" if "Subsystem" in df.columns else "Metabolic Pathways"
    metabolite_column = "Metabolite Shortname" if "Metabolite Shortname" in df.columns else "Metabolites"

    heatmap_df = df.pivot_table(
        index=pathway_column,
        columns=metabolite_column,
        values="log2FC_fluxsum",
        aggfunc="mean"
    )

    if heatmap_df.empty:
        print(f"[plot_fluxsum_log2fc_heatmap] Empty pivot for {tissue} — skip.")
        return

    # Apply pathway filtering if pathways_filter is provided
    if pathways_filter is not None:
        available_pathways = [p for p in pathways_filter if p in heatmap_df.index]
        if available_pathways:
            heatmap_df = heatmap_df.reindex(available_pathways)
            logging.info(f"[{tissue}] mcPFA pathway filtering: {len(heatmap_df.index)} pathways selected")
        else:
            logging.warning(f"[{tissue}] No pathways matched filter; using all pathways.")

    # Apply log2fc_threshold filtering: keep only metabolites with significant changes
    # in at least min_pathways_with_change pathways
    # This ensures metabolite changes are reproducible across multiple pathways
    if log2fc_threshold > 0:
        initial_metabolites = len(heatmap_df.columns)

        metabolites_to_keep = []
        for met in heatmap_df.columns:
            met_values = heatmap_df[met].dropna()
            # Count how many pathways have |log2FC| > threshold for this metabolite
            n_significant = (met_values.abs() > log2fc_threshold).sum()
            if n_significant >= min_pathways_with_change:
                metabolites_to_keep.append(met)

        heatmap_df = heatmap_df[metabolites_to_keep]

        logging.info(
            f"[{tissue}] mcPFA log2fc filtering (threshold={log2fc_threshold}, min_pathways={min_pathways_with_change}): "
            f"metabolites {initial_metabolites} → {len(metabolites_to_keep)}"
        )

    heatmap_df = filter_metabolites_by_content(heatmap_df, min_non_empty=min_metabolite_values)

    if heatmap_df.empty:
        print(f"[plot_fluxsum_log2fc_heatmap] No data after filtering for {tissue} — skip.")
        return

    # Filter out pathways that became empty after metabolite filtering
    pathway_content = heatmap_df.notna().sum(axis=1)
    heatmap_df = heatmap_df[pathway_content > 0]

    if heatmap_df.empty:
        print(f"[plot_fluxsum_log2fc_heatmap] No pathways with data after filtering for {tissue} — skip.")
        return

    if metabolite_filter:
        met_filter_mode = metabolite_filter.get("mode", "none")
        if met_filter_mode == "whitelist":
            whitelist = metabolite_filter.get("whitelist", [])
            if whitelist:
                available = [m for m in whitelist if m in heatmap_df.columns]
                if available:
                    heatmap_df = heatmap_df[available]
                    logging.info(f"[{tissue}] Metabolite whitelist: {len(available)} of {len(whitelist)} matched")
                else:
                    logging.warning(f"[{tissue}] Metabolite whitelist: no matches found, showing all")
        elif met_filter_mode == "blacklist":
            blacklist = metabolite_filter.get("blacklist", [])
            if blacklist:
                cols_to_keep = [c for c in heatmap_df.columns if c not in blacklist]
                heatmap_df = heatmap_df[cols_to_keep]
                logging.info(f"[{tissue}] Metabolite blacklist: removed {len(blacklist)} metabolites, {len(cols_to_keep)} remaining")

    if heatmap_df.empty:
        print(f"[plot_fluxsum_log2fc_heatmap] No data after metabolite filtering for {tissue} — skip.")
        return

    # Filter out pathways that became empty after all metabolite filtering
    pathway_content = heatmap_df.notna().sum(axis=1)
    heatmap_df = heatmap_df[pathway_content > 0]

    if heatmap_df.empty:
        print(f"[plot_fluxsum_log2fc_heatmap] No pathways with data after all filtering for {tissue} — skip.")
        return

    if merge_identical:
        logging.info(f"[{tissue}] Before merging: {len(heatmap_df.columns)} metabolites")
        heatmap_df = merge_identical_metabolites(heatmap_df)
        logging.info(f"[{tissue}] After merging: {len(heatmap_df.columns)} metabolite groups")

    # Use tissue-specific figure size or default
    if figsize and len(figsize) == 2:
        fig_width, fig_height = figsize
    else:
        # Default size
        fig_width = 14
        fig_height = 10

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    if tissue == "liver" or tissue == "leg":
        sns.heatmap(
            heatmap_df,
            cmap="coolwarm",
            center=0,
            linewidths=0.3,
            linecolor="gray",
            cbar_kws={"label": "log2FC", "pad": 0.02, "shrink": 0.7},
            annot=True,
            fmt=".0f",
            annot_kws={"size": annot_fontsize, "family": "Times New Roman"},
            ax=ax
        )

    else:
        sns.heatmap(
            heatmap_df,
            cmap="coolwarm",
            center=0,
            linewidths=0.3,
            linecolor="gray",
            cbar_kws={"label": "log2FC", "pad": 0.02, "shrink": 0.7},
            annot=True,
            fmt=".0f",
            annot_kws={"size": annot_fontsize, "family": "Times New Roman"},
            ax=ax
        )

    tick_rotation = -130
    leader_angle = tick_rotation
    leader_len_pts = 5
    leader_lw = 0.9
    leader_color = "gray"

    x_centers = np.arange(heatmap_df.shape[1]) + 0.5

    dx = leader_len_pts * math.cos(math.radians(leader_angle))
    dy = leader_len_pts * math.sin(math.radians(leader_angle))

    for x in x_centers:
        ax.annotate(
            '', xy=(x, 0), xycoords=ax.get_xaxis_transform(),
            xytext=(dx, dy), textcoords='offset points',
            arrowprops=dict(arrowstyle='-', lw=leader_lw, color=leader_color,
                            shrinkA=0, shrinkB=0),
            annotation_clip=False
        )

    ax.set_xlabel("Metabolites", fontsize=axis_title_fontsize, fontfamily='Times New Roman')
    ax.set_ylabel("Metabolic pathways", fontsize=axis_title_fontsize, fontfamily='Times New Roman')

    ax.tick_params(axis='x', rotation=45, labelsize=label_fontsize, pad=0)
    ax.tick_params(axis='y', rotation=0, labelsize=label_fontsize, pad=0)

    for label in ax.get_xticklabels():
        label.set_fontfamily('Times New Roman')
        label.set_fontsize(label_fontsize)
        label.set_horizontalalignment('right')

    for label in ax.get_yticklabels():
        label.set_fontfamily('Times New Roman')
        label.set_fontsize(label_fontsize)

    cbar = ax.collections[0].colorbar
    cbar.ax.set_ylabel('log2FC', fontsize=cbar_label_fontsize, fontfamily='Times New Roman')
    cbar.ax.tick_params(labelsize=cbar_tick_fontsize)
    for label in cbar.ax.get_yticklabels():
        label.set_fontfamily('Times New Roman')
        label.set_fontsize(cbar_tick_fontsize)

    plt.tight_layout()
    out_png = os.path.join(output_dir, 'plots', f"mcPFA_heatmap_{tissue}.png")
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()

    logging.info(f"mcPFA heatmap saved: {out_png}")

    # Save filtered version of original dataframe (not pivot table)
    # Get list of metabolites that passed filtering
    filtered_metabolites = list(heatmap_df.columns)
    filtered_pathways = list(heatmap_df.index)

    # Filter original dataframe by metabolites and pathways that appear on heatmap
    df_filtered = df[
        (df[metabolite_column].isin(filtered_metabolites)) &
        (df[pathway_column].isin(filtered_pathways))
    ].copy()

    csv_path = os.path.join(output_dir, f"mcPFA_filtered_{tissue}.csv")
    df_filtered.to_csv(csv_path, index=False)
    logging.info(f"Filtered mcPFA data saved: {csv_path} ({len(df_filtered)} rows, {df_filtered[metabolite_column].nunique()} metabolites, {df_filtered[pathway_column].nunique()} pathways)")



def plot_fluxsum_log2fc_heatmap_mets(df: pd.DataFrame, tissue: str, output_dir: str,
                                metabolites_filter=None, min_pathway_values=1,
                                merge_identical: bool = True):
    """
    Enhanced heatmap for mcPFA filtered by key metabolites instead of pathways,
    with empty pathway filtering and identical metabolite merging.
    """
    if df is None or df.empty:
        print(f"[plot_fluxsum_log2fc_heatmap] Empty DataFrame for {tissue} — skip.")
        return

    os.makedirs(output_dir, exist_ok=True)

    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.sans-serif'] = ['Arial']

    pathway_column = "Subsystem" if "Subsystem" in df.columns else "Metabolic Pathways"
    metabolite_column = "Metabolite Shortname" if "Metabolite Shortname" in df.columns else "Metabolites"

    heatmap_df = df.pivot_table(
        index=pathway_column,
        columns=metabolite_column,
        values="log2FC_fluxsum",
        aggfunc="mean"
    )

    if heatmap_df.empty:
        print(f"[plot_fluxsum_log2fc_heatmap] Empty pivot for {tissue} — skip.")
        return

    if metabolites_filter is None:
        metabolites_filter = [
            'ATP', 'ADP', 'NADH', 'NADPH', 'Flavin adenine dinucleotide reduced',
            'pyruvate', '(S)-lactate', 'L-Lactate', 'D-Glucose 6-phosphate', 'D-ribulose 5-phosphate(2-)',
            'Acetyl-CoA', 'Citrate', '(S)-malate(2-)', 'Succinate', '2-Oxoglutarate',
            'L-glutamate(-1)', 'L-glutamine', 'L-aspartate(1-)'
        ]

    available_metabolites = [m for m in metabolites_filter if m in heatmap_df.columns]
    if available_metabolites:
        heatmap_df = heatmap_df[available_metabolites]
        logging.info(f"[{tissue}] Filtered to {len(available_metabolites)} key metabolites: {available_metabolites}")
    else:
        logging.warning(f"[{tissue}] No metabolites matched filter; using all metabolites.")

    pathway_content = heatmap_df.notna().sum(axis=1)
    heatmap_df = heatmap_df[pathway_content >= min_pathway_values]

    if heatmap_df.empty:
        print(f"[plot_fluxsum_log2fc_heatmap] No data after filtering for {tissue} — skip.")
        return

    if merge_identical:
        logging.info(f"[{tissue}] Before merging: {len(heatmap_df.columns)} metabolites")
        heatmap_df = merge_identical_metabolites(heatmap_df)
        logging.info(f"[{tissue}] After merging: {len(heatmap_df.columns)} metabolite groups")

    n_pathways = len(heatmap_df.index)
    fig_width = 10
    fig_height = max(3, n_pathways * 0.3)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    sns.heatmap(
        heatmap_df,
        cmap="coolwarm",
        center=0,
        linewidths=0.3,
        linecolor="gray",
        cbar_kws={"label": "log2FC", "pad": 0.02, "shrink": 0.7},
        annot=True,
        fmt=".1f",
        annot_kws={"size": 7, "family": "Arial"},
        ax=ax
    )

    tick_rotation = -45
    leader_angle = tick_rotation
    leader_len_pts = 5
    leader_lw = 0.9
    leader_color = "gray"

    ax.tick_params(axis='x', rotation=tick_rotation, labelsize=8, pad=-3)

    offset = mtransforms.ScaledTranslation(-5/72, 0, fig.dpi_scale_trans)
    for label in ax.get_xticklabels():
        label.set_fontfamily('Arial')
        label.set_horizontalalignment('right')
        label.set_rotation_mode('anchor')
        label.set_transform(label.get_transform() + offset)

    x_centers = np.arange(heatmap_df.shape[1]) + 0.5

    dx = leader_len_pts * math.cos(math.radians(leader_angle))
    dy = leader_len_pts * math.sin(math.radians(leader_angle))

    for x in x_centers:
        ax.annotate(
            '', xy=(x, 0), xycoords=ax.get_xaxis_transform(),
            xytext=(dx, dy), textcoords='offset points',
            arrowprops=dict(arrowstyle='-', lw=leader_lw, color=leader_color,
                            shrinkA=0, shrinkB=0),
            annotation_clip=False
        )

    ax.set_xlabel("Key Metabolites", fontsize=10, fontfamily='Arial', fontweight='bold')
    ax.set_ylabel("Metabolic Pathways", fontsize=10, fontfamily='Arial', fontweight='bold')

    ax.tick_params(axis='y', rotation=0, labelsize=9)

    for label in ax.get_yticklabels():
        label.set_fontfamily('Arial')

    cbar = ax.collections[0].colorbar
    cbar.ax.set_ylabel('log2FC', fontsize=9, fontfamily='Arial')
    cbar.ax.tick_params(labelsize=8)
    for label in cbar.ax.get_yticklabels():
        label.set_fontfamily('Arial')

    plt.tight_layout()

    out_png = os.path.join(output_dir, f"mcPFA_heatmap_keymetabolites_{tissue}.png")
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()

    logging.info(f"mcPFA heatmap (key metabolites) saved: {out_png}")

    csv_path = os.path.join(output_dir, f"mcPFA_keymetabolites_{tissue}.csv")
    heatmap_df.to_csv(csv_path)
    logging.info(f"Filtered mcPFA data (key metabolites) saved: {csv_path}")


def create_flux_summary_table(tissue_flux_dfs: dict, output_dir: str):
    """
    Create summary table of flux values for all tissues.
    """
    os.makedirs(output_dir, exist_ok=True)

    summary_rows = []
    for tissue, df in tissue_flux_dfs.items():
        if df is not None and not df.empty:
            df_copy = df.copy()
            df_copy['Tissue'] = tissue
            summary_rows.append(df_copy)

    if not summary_rows:
        logging.warning("No flux data to create summary table")
        return

    summary_table = pd.concat(summary_rows, ignore_index=True)
    summary_table = summary_table[['Tissue', 'Pathways', '∑FG_flux', '∑SG_flux', 'Difference']]
    summary_table = summary_table.sort_values(['Tissue', 'Pathways'])

    for col in ['∑FG_flux', '∑SG_flux', 'Difference']:
        summary_table[col] = summary_table[col].round(1)

    out_path = os.path.join(output_dir, 'pathway_flux_summary_all_tissues.csv')
    summary_table.to_csv(out_path, index=False)
    logging.info(f"Summary flux table saved: {out_path}")

    return summary_table
