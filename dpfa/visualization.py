import os
import logging
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import math
import matplotlib.transforms as mtransforms

def _explode_pathways(df: pd.DataFrame) -> pd.DataFrame:
    """
    Explode pathways column and remove duplicate pathway entries per reaction.

    After pathway merging, a reaction may have duplicate pathway names.
    We deduplicate to avoid counting reactions multiple times in DRF bar plots.
    """
    if "Pathways" not in df.columns:
        return df

    out = df.copy()
    out["Pathways"] = out["Pathways"].fillna("").astype(str)
    out["Pathways"] = (out["Pathways"]
                       .str.replace("|", ";", regex=False)
                       .str.replace(" / ", ";", regex=False))
    out["Pathways"] = out["Pathways"].str.split(";")
    out = out.explode("Pathways")
    out["Pathways"] = out["Pathways"].str.strip()
    out = out[out["Pathways"] != ""]

    if "reaction_id" in out.columns:
        out = out.drop_duplicates(subset=["reaction_id", "Pathways"], keep="first")

    return out

def plot_regulation_counts(df: pd.DataFrame, output_dir: str, tissue: str,
                           pathways_filter=None, pathway_flux_diff=None,
                           flux_diff_threshold: float = 0.0,
                           xlim: tuple = None):
    """
    Plot DRF histogram by pathways with variable bar width.
    Bar width depends on |flux_difference| (normalized by maximum).
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

    if 'Pathways' not in df.columns or 'DRF_category' not in df.columns:
        logging.warning("No pathways or DRF_category for plotting.")
        return

    df = _explode_pathways(df)

    if pathways_filter is None:
        pathways_filter = UNIFIED_PATHWAYS

    df = df[df['Pathways'].isin(pathways_filter)]
    df = df[df['DRF_category'].isin(color_dict.keys())]

    if pathway_flux_diff is not None:
        flux_table = pathway_flux_diff[['Pathways', 'flux_fast', 'flux_slow', 'flux_difference']].copy()
        flux_table = flux_table.rename(columns={
            'flux_fast': '∑FG_flux',
            'flux_slow': '∑SG_flux',
            'flux_difference': 'Difference'
        })
        flux_table_path = os.path.join(output_dir, f'pathway_flux_values_{tissue}.csv')
        flux_table.to_csv(flux_table_path, index=False)
        logging.info(f"Flux values table saved: {flux_table_path}")

    df['Pathways label'] = df['Pathways']
    counts = df.groupby(['Pathways label', 'DRF_category']).size().unstack(fill_value=0)
    pathway_order = [p for p in pathways_filter if p in counts.index]
    counts = counts.reindex(pathway_order)

    if flux_diff_threshold > 0 and pathway_flux_diff is not None:
        flux_diff_dict = dict(zip(
            pathway_flux_diff['Pathways'],
            abs(pathway_flux_diff['flux_difference'])
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
        flux_diff_dict = dict(zip(
            pathway_flux_diff['Pathways'],
            abs(pathway_flux_diff['flux_difference'])
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

    if False:
        fig, ax = plt.subplots(figsize=(6, 2.2))
    else:
        fig, ax = plt.subplots(figsize=(6, 2.5))

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
        flux_diff_dict = dict(zip(
            pathway_flux_diff['Pathways'],
            pathway_flux_diff['flux_difference']
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

def analyze_pathway_flux_difference(merged_df: pd.DataFrame, output_dir: str,
                                   tissue: str, threshold: float = 1.0) -> pd.DataFrame:
    """Analyze flux difference by pathways (sum of absolute fluxes)."""
    if not {'Pathways', 'flux_slow', 'flux_fast'}.issubset(merged_df.columns):
        logging.warning("Pathway flux analysis skipped (missing flux data).")
        return None

    merged_df = _explode_pathways(merged_df)

    path_sums = merged_df.groupby('Pathways')[['flux_slow', 'flux_fast']].apply(
        lambda x: x.abs().sum()
    ).reset_index()

    path_sums['flux_slow'] = path_sums['flux_slow'].where(
        path_sums['flux_slow'] >= 1e-6, 0
    )
    path_sums['flux_fast'] = path_sums['flux_fast'].where(
        path_sums['flux_fast'] >= 1e-6, 0
    )

    path_sums['flux_difference'] = path_sums['flux_fast'] - path_sums['flux_slow']

    filtered = path_sums[abs(path_sums['flux_difference']) > threshold]
    filtered = filtered.sort_values(by='flux_difference', ascending=False)

    out_path = os.path.join(output_dir, f'pathway_flux_diff_{tissue}.csv')
    path_sums.to_csv(out_path, index=False)
    logging.info(f"Pathway flux difference saved: {out_path}")

    return filtered


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
                                merge_identical: bool = True):
    """
    Enhanced heatmap for mcPFA with unified pathways, empty metabolite filtering,
    identical metabolite merging, and Times New Roman font for compact display.
    """
    if df is None or df.empty:
        print(f"[plot_fluxsum_log2fc_heatmap] Empty DataFrame for {tissue} — skip.")
        return

    os.makedirs(output_dir, exist_ok=True)

    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times New Roman', 'Times', 'DejaVu Serif']

    label_fontsize = 7.5
    axis_title_fontsize = 9
    annot_fontsize = 8
    cbar_label_fontsize = 8
    cbar_tick_fontsize = 7

    pathway_column = "Pathway Group" if "Pathway Group" in df.columns else "Metabolic Pathways"

    heatmap_df = df.pivot_table(
        index=pathway_column,
        columns="Metabolites",
        values="log2FC_fluxsum",
        aggfunc="mean"
    )

    if heatmap_df.empty:
        print(f"[plot_fluxsum_log2fc_heatmap] Empty pivot for {tissue} — skip.")
        return

    if pathways_filter is None:
        pathways_filter = UNIFIED_PATHWAYS

    available_pathways = [p for p in pathways_filter if p in heatmap_df.index]
    if available_pathways:
        heatmap_df = heatmap_df.reindex(available_pathways)
    else:
        logging.warning(f"[{tissue}] No pathways matched filter; using all pathways.")

    heatmap_df = filter_metabolites_by_content(heatmap_df, min_non_empty=min_metabolite_values)

    if heatmap_df.empty:
        print(f"[plot_fluxsum_log2fc_heatmap] No data after filtering for {tissue} — skip.")
        return

    if merge_identical:
        logging.info(f"[{tissue}] Before merging: {len(heatmap_df.columns)} metabolites")
        heatmap_df = merge_identical_metabolites(heatmap_df)
        logging.info(f"[{tissue}] After merging: {len(heatmap_df.columns)} metabolite groups")

    if str.lower(tissue) == "liver":
        fig_width = 7
        fig_height = 2.5
    elif str.lower(tissue) == "breast":
        fig_width = 7
        fig_height = 3
    else:
        fig_width = 7
        fig_height = 2.7

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

    out_png = os.path.join(output_dir, f"mcPFA_heatmap_{tissue}.png")
    plt.savefig(out_png, dpi=300, bbox_inches='tight')
    plt.close()

    logging.info(f"mcPFA heatmap saved: {out_png}")

    csv_path = os.path.join(output_dir, f"mcPFA_filtered_{tissue}.csv")
    heatmap_df.to_csv(csv_path)
    logging.info(f"Filtered mcPFA data saved: {csv_path}")



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

    pathway_column = "Pathway Group" if "Pathway Group" in df.columns else "Metabolic Pathways"

    heatmap_df = df.pivot_table(
        index=pathway_column,
        columns="Metabolites",
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
