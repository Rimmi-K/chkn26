import os
import logging
from collections import defaultdict
from typing import Dict, Tuple, Optional, List
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.scale import FuncScale
from matplotlib.lines import Line2D
from matplotlib.patches import Wedge, Circle
from .utils.pathway_colors import canonicalize_pathway, palette_for_df, get_color_for_pathway
from .utils.gpr_utils import build_gpr_rule_map, make_gpr_rule_with_values, build_reaction_pathways_map
from .utils.flux_utils import FLUX_ON, FLUX_OFF, FLUX_REVERSED
import matplotlib.ticker as mticker


def parse_pathways(pw_str) -> List[str]:
    if pd.isna(pw_str) or pw_str == "Unknown":
        return []

    pw_str = str(pw_str)
    if pw_str.startswith("[") and pw_str.endswith("]"):
        inner = pw_str[1:-1].strip()
        parts = [p.strip(" '\"") for p in inner.split(",") if p.strip(" '\"")]
        return [canonicalize_pathway(p) for p in parts]

    for delim in [";", "|", " / "]:
        if delim in pw_str:
            paths = [p.strip() for p in pw_str.split(delim) if p.strip()]
            return [canonicalize_pathway(p) for p in paths]

    paths = [pw_str.strip()] if pw_str.strip() else []

    return [canonicalize_pathway(p) for p in paths]

def draw_pie_marker(ax, x, y, colors: List[str], size: float,
                   n_slices: int = None, edgecolor='black', linewidth=0.5, alpha=0.95):
    """
    Draws a pie-chart marker with multiple colors
    """
    if not colors or len(colors) < 2:
        ax.scatter([x], [y], s=size, c=colors[0] if colors else 'gray',
                  edgecolors=edgecolor, linewidth=linewidth, alpha=alpha)
        return

    fig = ax.get_figure()

    radius_points = np.sqrt(size / np.pi)
    xy_display = ax.transData.transform([(x, y)])[0]

    dpi = fig.dpi
    radius_pixels = radius_points * (dpi / 72.0)

    n = n_slices if n_slices else len(colors)
    angle_step = 360 / n

    for i in range(n):
        color = colors[i % len(colors)]
        theta1 = i * angle_step
        theta2 = (i + 1) * angle_step

        wedge = Wedge(
            xy_display, radius_pixels,
            theta1, theta2,
            facecolor=color, edgecolor=edgecolor,
            linewidth=linewidth, alpha=alpha,
            transform=None
        )
        wedge.set_transform(ax.transData.inverted() + ax.transData)
        ax.add_patch(wedge)


def draw_split_marker(ax, x, y, colors: List[str], size: float,
                     split_type: str = 'vertical', edgecolor='black',
                     linewidth=0.5, alpha=0.95):
    """
    Draws a split marker (for 2 pathways)
    """
    if len(colors) < 2:
        ax.scatter([x], [y], s=size, c=colors[0] if colors else 'gray',
                  edgecolors=edgecolor, linewidth=linewidth, alpha=alpha)
        return

    fig = ax.get_figure()
    radius_points = np.sqrt(size / np.pi)
    xy_display = ax.transData.transform([(x, y)])[0]

    dpi = fig.dpi
    radius_pixels = radius_points * (dpi / 72.0)

    if split_type == 'vertical':
        wedge1 = Wedge(xy_display, radius_pixels, 90, 270,
                      facecolor=colors[0], edgecolor=edgecolor,
                      linewidth=linewidth, alpha=alpha,
                      transform=None)
        wedge1.set_transform(ax.transData.inverted() + ax.transData)

        wedge2 = Wedge(xy_display, radius_pixels, 270, 450,
                      facecolor=colors[1], edgecolor=edgecolor,
                      linewidth=linewidth, alpha=alpha,
                      transform=None)
        wedge2.set_transform(ax.transData.inverted() + ax.transData)

    elif split_type == 'horizontal':
        wedge1 = Wedge(xy_display, radius_pixels, 0, 180,
                      facecolor=colors[0], edgecolor=edgecolor,
                      linewidth=linewidth, alpha=alpha,
                      transform=None)
        wedge1.set_transform(ax.transData.inverted() + ax.transData)

        wedge2 = Wedge(xy_display, radius_pixels, 180, 360,
                      facecolor=colors[1], edgecolor=edgecolor,
                      linewidth=linewidth, alpha=alpha,
                      transform=None)
        wedge2.set_transform(ax.transData.inverted() + ax.transData)

    else:  # diagonal
        wedge1 = Wedge(xy_display, radius_pixels, 45, 225,
                      facecolor=colors[0], edgecolor=edgecolor,
                      linewidth=linewidth, alpha=alpha,
                      transform=None)
        wedge1.set_transform(ax.transData.inverted() + ax.transData)

        wedge2 = Wedge(xy_display, radius_pixels, 225, 405,
                      facecolor=colors[1], edgecolor=edgecolor,
                      linewidth=linewidth, alpha=alpha,
                      transform=None)
        wedge2.set_transform(ax.transData.inverted() + ax.transData)

    ax.add_patch(wedge1)
    ax.add_patch(wedge2)




def _size_from_q_discrete_thresholds(q: float, thresholds: list[float],
                                     sizes: list[float], default_size: float) -> float:
    if q is None or pd.isna(q) or not thresholds or not sizes or len(sizes) != len(thresholds):
        return float(default_size)
    t, v = sorted(map(float, thresholds)), list(map(float, sizes))
    q = float(q)
    for i, thr in enumerate(t):
        if q <= thr:
            return v[i]
    return float(default_size)


def _collect_used_pathways_fallback(df: pd.DataFrame,
                                    rxn_lfc_thr: float,
                                    flux_log2_thr: float) -> pd.DataFrame:
    """Collects pathways for fallback case"""
    m = df["deg_log2fc"].abs() >= float(rxn_lfc_thr)
    m &= df["log2_ratio"].abs() >= float(flux_log2_thr)
    pw = (df.loc[m, "Pathways"]
            .dropna()
            .map(canonicalize_pathway))
    if pw.empty:
        return pd.DataFrame(columns=["Pathways","n_highlighted","color_hex"])
    used = (pw.value_counts()
              .rename_axis("Pathways")
              .reset_index(name="n_highlighted"))
    used["color_hex"] = used["Pathways"].map(get_color_for_pathway)
    return used


def make_band_squash_transform(B: float = 1.0, k: float = 8.0) -> Tuple:
    C = (1 - 1/k) * B
    def f(x):
        x = np.asarray(x)
        m = (x >= -B) & (x <= B)
        out = np.where(m, x / k, np.sign(x) * (np.abs(x) - C))
        return out
    def finv(y):
        y = np.asarray(y)
        m = (y >= -B/k) & (y <= B/k)
        out = np.where(m, y * k, np.sign(y) * (np.abs(y) + C))
        return out
    return f, finv


def _to_log2_ratio(v: float) -> Optional[float]:
    if pd.isna(v): return None
    if v == FLUX_ON: return 3.0
    if v == FLUX_OFF: return -3.0
    r = np.sign(v) * max(abs(round(float(v), 2)), 0.01)
    try: return float(np.log2(r))
    except Exception: return None


def _size_from_q(q: Optional[float], s_min: float = 40.0, s_max: float = 120.0,
                q_floor: float = 1e-6, q_ceil: float = 0.25) -> float:
    if q is None or pd.isna(q): return s_min
    q_clamped = min(max(float(q), q_floor), q_ceil)
    lo, hi = -np.log10(q_ceil), -np.log10(q_floor)
    t = np.clip((-np.log10(q_clamped) - lo) / (hi - lo), 0.0, 1.0)
    return s_min + t * (s_max - s_min)


def _prepare_scatter_dataframe(
    merged_df: pd.DataFrame,
    rxn_lfc_map: Dict[str, float],
    rxn_sig_map: Dict[str, bool],
    rxn_p_map: Optional[Dict[str, float]],
    padj_rxn_map: Optional[Dict[str, float]],
    model_for_paths,
    pathway_merging: Optional[Dict[str, List[str]]],
    pathway_filter: Optional[List[str]],
    rxn_lfc_thr: float,
    flux_log2_thr: float,
    model_for_gpr_rules,
    deg_df_for_labels,
    tissue: str,
) -> pd.DataFrame:
    """Prepares the scatter plot dataframe with all necessary columns"""
    total_all = len(merged_df)
    logging.info(f"[{tissue}] total reactions in merged_df: {total_all}")

    base_cols = ["reaction_id", "flux_ratio"]
    df = merged_df.loc[merged_df["flux_ratio"] != FLUX_REVERSED, base_cols].copy()
    removed_reversed = total_all - len(df)
    logging.info(f"[{tissue}] removed REVERSED: {removed_reversed}")

    df["log2_ratio"] = df["flux_ratio"].apply(_to_log2_ratio)
    df["deg_log2fc"] = df["reaction_id"].map(rxn_lfc_map)

    df["rxn_p"] = df["reaction_id"].map(rxn_p_map) if rxn_p_map is not None else np.nan
    df["padj_rxn"] = df["reaction_id"].map(padj_rxn_map) if padj_rxn_map is not None else np.nan
    df["rxn_sig_by_padj"] = df["reaction_id"].map(rxn_sig_map).fillna(False).astype(bool)

    if "Pathways" in merged_df.columns:
        path_map = dict(zip(merged_df["reaction_id"], merged_df["Pathways"]))
        df["Pathways"] = df["reaction_id"].map(path_map)
    elif model_for_paths is not None:
        rxn2path = build_reaction_pathways_map(model_for_paths)
        df["Pathways"] = df["reaction_id"].map(rxn2path)
    else:
        df["Pathways"] = "Unknown"

    if "Pathways_list" in merged_df.columns:
        list_map = dict(zip(merged_df["reaction_id"], merged_df["Pathways_list"]))
        df["pathway_list"] = df["reaction_id"].map(list_map)
        df["pathway_list"] = df["pathway_list"].apply(
            lambda v: v if isinstance(v, list) else parse_pathways(v)
        )
    else:
        df["pathway_list"] = df["Pathways"].apply(parse_pathways)

    if pathway_merging:
        reverse_mapping = {}
        for new_name, old_names in pathway_merging.items():
            for old_name in old_names:
                reverse_mapping[old_name] = new_name

        def merge_pathway_list(pw_list):
            if not isinstance(pw_list, list):
                return pw_list
            merged = [reverse_mapping.get(pw, pw) for pw in pw_list]
            seen = set()
            result = []
            for pw in merged:
                if pw not in seen:
                    seen.add(pw)
                    result.append(pw)
            return result

        df["pathway_list"] = df["pathway_list"].apply(merge_pathway_list)
        logging.info(f"[{tissue}] Applied pathway merging: {len(pathway_merging)} groups")

    df["n_pathways"] = df["pathway_list"].apply(len)

    logging.info(f"[{tissue}] Pathway count distribution:")
    for n, count in df["n_pathways"].value_counts().sort_index().items():
        logging.info(f"  {n} pathways: {count} reactions")

    if pathway_filter is not None:
        def filter_pathway_list(pw_list):
            if not isinstance(pw_list, list):
                return pw_list
            filtered = [pw for pw in pw_list if pw in pathway_filter]
            return filtered

        before_filter = len(df)
        df["pathway_list"] = df["pathway_list"].apply(filter_pathway_list)
        df = df[df["pathway_list"].apply(len) > 0]
        df["n_pathways"] = df["pathway_list"].apply(len)
        logging.info(f"[{tissue}] Applied pathway filter: {before_filter} → {len(df)} reactions")

    before_dropna = len(df)
    df = df.dropna(subset=["deg_log2fc", "log2_ratio"])
    logging.info(f"[{tissue}] dropped NaNs: {before_dropna - len(df)}; remain: {len(df)}")

    if df.empty:
        return df

    df["pass_lfc"] = df["deg_log2fc"].abs() >= float(rxn_lfc_thr)
    df["pass_flux"] = df["log2_ratio"].abs() >= float(flux_log2_thr)
    df["highlight"] = df["pass_lfc"] & df["pass_flux"]

    logging.info(f"[{tissue}] highlighted: {int(df['highlight'].sum())} of {len(df)}")

    if model_for_gpr_rules is not None:
        rxn2rule = build_gpr_rule_map(model_for_gpr_rules)
        df["GPR_rule"] = df["reaction_id"].map(rxn2rule).fillna("")
        if deg_df_for_labels is not None and {"gene","log2fc"}.issubset(deg_df_for_labels.columns):
            lfc_map = dict(zip(deg_df_for_labels["gene"], deg_df_for_labels["log2fc"]))
            rxn2genes = {rxn.id: [g.id for g in rxn.genes] for rxn in model_for_gpr_rules.reactions}
            def _rule_with_vals(rid: str) -> str:
                rule = rxn2rule.get(rid, "")
                genes = rxn2genes.get(rid, [])
                return make_gpr_rule_with_values(rule, genes, lfc_map,
                                                na_token="NA", fmt="{:.2f}")
            df["GPR_rule_log2FC"] = df["reaction_id"].map(_rule_with_vals)
        else:
            df["GPR_rule_log2FC"] = ""

    return df


def _calculate_marker_sizes(
    df: pd.DataFrame,
    size_mode: str,
    size_thresholds: Optional[list[float]],
    size_values: Optional[list[float]],
    size_default: float,
) -> pd.DataFrame:
    """Calculates marker sizes based on p-values"""
    df["size_s"] = np.nan
    hi_mask = df["highlight"]

    if size_mode == "discrete" and size_thresholds and size_values and len(size_thresholds) == len(size_values):
        df.loc[hi_mask, "size_s"] = df.loc[hi_mask, "padj_rxn"].apply(
            lambda q: _size_from_q_discrete_thresholds(q, size_thresholds, size_values, size_default)
        )
    else:
        df.loc[hi_mask, "size_s"] = df.loc[hi_mask, "padj_rxn"].apply(
            lambda q: _size_from_q(q, s_min=40.0, s_max=120.0, q_floor=1e-6, q_ceil=0.25)
        )

    return df


def _draw_scatter_with_markers(ax, df: pd.DataFrame, tissue: str) -> List[dict]:
    """Draws scatter plot with single and multi-pathway markers"""
    df_lo = df[~df["highlight"]]
    if not df_lo.empty:
        sns.scatterplot(
            data=df_lo, x="deg_log2fc", y="log2_ratio",
            color="lightgrey", alpha=0.35, s=22, edgecolor="black",
            linewidth=0.3, ax=ax
        )

    df_hi = df[df["highlight"]].copy()
    multi_markers_data = []

    if not df_hi.empty:
        df_hi["size_s"] = df_hi["size_s"].fillna(63.0)

        df_single = df_hi[df_hi["n_pathways"] == 1].copy()
        df_multi = df_hi[df_hi["n_pathways"] >= 2].copy()

        if not df_single.empty:
            for pathway in df_single["pathway_list"].apply(lambda x: x[0] if x else "Unknown").unique():
                mask = df_single["pathway_list"].apply(lambda x: x[0] if x else "Unknown") == pathway
                df_pw = df_single[mask]
                color = get_color_for_pathway(pathway)
                ax.scatter(
                    df_pw["deg_log2fc"],
                    df_pw["log2_ratio"],
                    s=df_pw["size_s"],
                    c=color,
                    alpha=0.95,
                    edgecolors='black',
                    linewidth=0.5
                )

        if not df_multi.empty:
            logging.info(f"[{tissue}] Sample multi-pathway reactions:")
            for idx, row in df_multi.head(5).iterrows():
                x = row["deg_log2fc"]
                y = row["log2_ratio"]
                size = row["size_s"]
                pathways = row["pathway_list"]
                rxn_id = row["reaction_id"]

                top_pathways = pathways[:3]
                colors = [get_color_for_pathway(pw) for pw in top_pathways]

                logging.info(f"  {rxn_id}: pathways={top_pathways}, colors={colors}")

                multi_markers_data.append({
                    'x': x, 'y': y, 'size': size, 'colors': colors
                })

            for idx, row in df_multi.iloc[5:].iterrows():
                x = row["deg_log2fc"]
                y = row["log2_ratio"]
                size = row["size_s"]
                pathways = row["pathway_list"]

                top_pathways = pathways[:3]
                colors = [get_color_for_pathway(pw) for pw in top_pathways]

                multi_markers_data.append({
                    'x': x, 'y': y, 'size': size, 'colors': colors
                })

    return multi_markers_data


def _configure_axes(ax, tissue: str, rxn_lfc_thr: float, flux_log2_thr: float):
    """Configures axis scales, ticks, and grid"""
    fx, fxi = make_band_squash_transform(B=1.0, k=8.0)
    ax.set_xscale(FuncScale(ax, (fx, fxi)))
    ax.set_yscale(FuncScale(ax, (fx, fxi)))

    xt = np.arange(-1, 4, 1)
    yt = np.array([-4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7])

    ax.set_xticks(xt)
    ax.set_yticks(yt)
    ax.xaxis.set_major_locator(mticker.FixedLocator(xt))
    ax.yaxis.set_major_locator(mticker.FixedLocator(yt))

    ax.xaxis.set_major_formatter(mticker.FuncFormatter(
        lambda v, pos: f"{int(v)}" if float(v).is_integer() else f"{v}"
    ))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(
        lambda v, pos: f"{int(v)}" if float(v).is_integer() else f"{v}"
    ))

    ax.xaxis.set_minor_locator(mticker.NullLocator())
    ax.yaxis.set_minor_locator(mticker.NullLocator())

    if tissue == "leg":
        ax.set_ylim(bottom=-3.5, top=7)
    else:
        ax.set_ylim(bottom=-2, top=7)

    ax.axvline(0, ls="--", lw=1, color="k", alpha=0.35)
    ax.axhline(0, ls="--", lw=1, color="k", alpha=0.35)
    if rxn_lfc_thr > 0:
        ax.axvline(+rxn_lfc_thr, ls="--", lw=1, color="k", alpha=0.5)
        ax.axvline(-rxn_lfc_thr, ls="--", lw=1, color="k", alpha=0.5)
    if flux_log2_thr > 0:
        ax.axhline(+flux_log2_thr, ls="--", lw=1, color="k", alpha=0.5)
        ax.axhline(-flux_log2_thr, ls="--", lw=1, color="k", alpha=0.5)

    ax.grid(which="major", linestyle="--", linewidth=0.6, alpha=0.35)

    ax.set_xlabel("DEG-derived reaction log2FC (via GPR)", fontsize=12)
    ax.set_ylabel("log2(FG/SG flux ratio)", fontsize=12)
    ax.set_title(f"DEG vs flux ratio — {tissue}", fontsize=13)
    ax.tick_params(axis='both', labelsize=10)


def _add_scatter_annotations(ax, df_hi: pd.DataFrame):
    """Adds annotations for highlighted reactions"""
    xs = df_hi["deg_log2fc"].to_numpy(float)
    ys = df_hi["log2_ratio"].to_numpy(float)
    xy_pix = ax.transData.transform(np.column_stack([xs, ys]))
    parent = np.arange(len(xy_pix))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[rj] = ri

    cell_x, cell_y = 10.0, 5.0
    cells = np.column_stack([
        np.floor(xy_pix[:,0] / cell_x).astype(int),
        np.floor(xy_pix[:,1] / cell_y).astype(int)
    ])
    buckets = defaultdict(list)
    for idx, c in enumerate(map(tuple, cells)):
        buckets[c].append(idx)

    neighbor_shifts = [(dx,dy) for dx in (-1,0,1) for dy in (-1,0,1)]
    for (cx,cy), idxs in buckets.items():
        cand = []
        for dx,dy in neighbor_shifts:
            cand += buckets.get((cx+dx, cy+dy), [])
        cand = np.array(cand, int)
        for i in idxs:
            dx = np.abs(xy_pix[i,0] - xy_pix[cand,0])
            dy = np.abs(xy_pix[i,1] - xy_pix[cand,1])
            close = cand[(dx <= 10.0) & (dy <= 5.0)]
            for j in close:
                union(i, j)

    clusters = defaultdict(list)
    ids = df_hi["reaction_id"].tolist()
    for i in range(len(xy_pix)):
        clusters[find(i)].append(i)

    for comp in clusters.values():
        comp = sorted(comp)
        i_left = comp[np.argmin(xy_pix[comp, 0])]
        anchor_pix = xy_pix[i_left].copy()
        anchor_data = ax.transData.inverted().transform(anchor_pix)
        cluster_ids = [ids[i] for i in comp]
        text = ", ".join(cluster_ids) if len(cluster_ids) <= 20 \
               else ", ".join(cluster_ids[:20]) + f", +{len(cluster_ids)-20}"
        ax.annotate(
            text,
            xy=anchor_data, xytext=(3, 2), textcoords="offset points",
            fontsize=11, alpha=0.95, clip_on=True,
            ha="left", va="bottom"
        )


def _draw_multicolor_markers(ax, multi_markers_data: List[dict]):
    """Draws multi-pathway pie markers"""
    for marker_info in multi_markers_data:
        x, y = marker_info['x'], marker_info['y']
        size, colors = marker_info['size'], marker_info['colors']
        draw_pie_marker(ax, x, y, colors, size)


def _save_legends(output_dir: str, tissue: str, size_mode: str,
                  size_thresholds: Optional[list[float]],
                  size_values: Optional[list[float]],
                  size_default: float):
    """Saves size and pathway legends"""
    out_dir = os.path.join(output_dir, "plots")

    # Size legend
    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    handles = []

    if size_mode == "discrete" and size_thresholds and size_values:
        for th, sz in zip(size_thresholds, size_values):
            handles.append(
                plt.scatter([], [], s=sz, color="gray", label=f"q ≤ {th:.2f}")
            )
        handles.append(
            plt.scatter([], [], s=size_default, color="gray",
                      label=f"q > {size_thresholds[-1]:.2f}")
        )
    else:
        q_values = [0.25, 0.1, 0.01, 0.001]
        for q in q_values:
            s = _size_from_q(q, s_min=40.0, s_max=120.0, q_floor=1e-6, q_ceil=0.25)
            handles.append(
                plt.scatter([], [], s=s, color="gray", label=f"q = {q}")
            )

    ax.legend(
        handles=handles,
        loc="center", frameon=False,
        title="Marker size (reaction padj)", fontsize=10, title_fontsize=11
    )
    ax.axis("off")
    out_size_legend = os.path.join(out_dir, f"size_legend_{tissue}.pdf")
    plt.savefig(out_size_legend, dpi=300, bbox_inches="tight")
    plt.close()
    logging.info(f"Size legend saved to {out_size_legend}")

    # Pathway legend
    fig, ax = plt.subplots(figsize=(4, 3))

    ax.text(0.1, 0.9, "Single pathway:", fontsize=11, weight='bold',
           transform=ax.transAxes)
    ax.scatter([0.05], [0.82], s=100, c='#008000', transform=ax.transAxes)
    ax.text(0.12, 0.82, "Example pathway", fontsize=10,
           transform=ax.transAxes, va='center')

    ax.text(0.1, 0.65, "Multiple pathways:", fontsize=11, weight='bold',
           transform=ax.transAxes)

    demo_colors = ['#008000', '#ff7f0e']
    demo_x, demo_y = 0.05, 0.57

    wedge1 = Wedge((demo_x, demo_y), 0.03, 0, 180,
                  facecolor=demo_colors[0], transform=ax.transAxes)
    wedge2 = Wedge((demo_x, demo_y), 0.03, 180, 360,
                  facecolor=demo_colors[1], transform=ax.transAxes)
    ax.add_patch(wedge1)
    ax.add_patch(wedge2)

    ax.text(0.12, 0.57, "2+ pathways (top 2 shown)", fontsize=10,
           transform=ax.transAxes, va='center')

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    out_multi_legend = os.path.join(out_dir, f"multicolor_legend_{tissue}.pdf")
    plt.savefig(out_multi_legend, dpi=300, bbox_inches="tight")
    plt.close()
    logging.info(f"Multi-color legend saved to {out_multi_legend}")


def _collect_used_pathways(df_hi: pd.DataFrame, df: pd.DataFrame,
                           tissue: str, rxn_lfc_thr: float,
                           flux_log2_thr: float) -> pd.DataFrame:
    """Collects and returns used pathways dataframe"""
    if not df_hi.empty:
        all_pathways = []
        for pw_list in df_hi["pathway_list"]:
            all_pathways.extend(pw_list)

        if all_pathways:
            pw_series = pd.Series(all_pathways)
            used_pw = (pw_series.value_counts()
                      .rename_axis("Pathways")
                      .reset_index(name="n_highlighted"))
            used_pw["color_hex"] = used_pw["Pathways"].map(get_color_for_pathway)
            used_pw["tissue"] = tissue
        else:
            used_pw = pd.DataFrame(columns=["Pathways","n_highlighted","color_hex","tissue"])
    else:
        used_pw = pd.DataFrame(columns=["Pathways","n_highlighted","color_hex","tissue"])

    if used_pw.empty:
        used_fb = _collect_used_pathways_fallback(df, rxn_lfc_thr, flux_log2_thr)
        if not used_fb.empty:
            used_fb["tissue"] = tissue
            used_pw = used_fb

    return used_pw


def make_scatter_deg_vs_flux(
    merged_df: pd.DataFrame,
    rxn_lfc_map: Dict[str, float],
    rxn_sig_map: Dict[str, bool],
    tissue: str,
    output_dir: str,
    model_for_paths=None,

    rxn_p_map: Optional[Dict[str, float]] = None,
    padj_rxn_map: Optional[Dict[str, float]] = None,

    size_mode: str = "continuous",
    size_thresholds: Optional[list[float]] = None,
    size_values: Optional[list[float]] = None,
    size_default: float = 40.0,

    pathway_merging: Optional[Dict[str, List[str]]] = None,
    pathway_filter: Optional[List[str]] = None,

    rxn_lfc_thr: float = 1.0,
    flux_log2_thr: float = 1.0,

    model_for_gpr_rules=None,
    deg_df_for_labels=None,
):
    """
    Draws scatter plot of DEG vs flux ratio with multi-colored markers
    for reactions belonging to multiple pathways
    """
    os.makedirs(os.path.join(output_dir, "plots"), exist_ok=True)

    # Prepare data
    df = _prepare_scatter_dataframe(
        merged_df, rxn_lfc_map, rxn_sig_map, rxn_p_map, padj_rxn_map,
        model_for_paths, pathway_merging, pathway_filter,
        rxn_lfc_thr, flux_log2_thr,
        model_for_gpr_rules, deg_df_for_labels, tissue
    )

    if df.empty:
        logging.warning("Scatter skipped: no data after mapping and filtering.")
        return pd.DataFrame(columns=["Pathways","n_highlighted","color_hex","tissue"])

    df = _calculate_marker_sizes(df, size_mode, size_thresholds, size_values, size_default)

    csv_out = os.path.join(output_dir, f"scatter_deg_vs_flux_{tissue}.csv")
    df.to_csv(csv_out, index=False)

    fig, ax = plt.subplots(figsize=(7, 6))

    multi_markers_data = _draw_scatter_with_markers(ax, df, tissue)

    _configure_axes(ax, tissue, rxn_lfc_thr, flux_log2_thr)

    df_hi = df[df["highlight"]]
    if not df_hi.empty and multi_markers_data:
        _draw_multicolor_markers(ax, multi_markers_data)

    if not df_hi.empty:
        _add_scatter_annotations(ax, df_hi)

    out_dir = os.path.join(output_dir, "plots")
    out_pdf = os.path.join(out_dir, f"deg_vs_flux_scatter_{tissue}.pdf")
    out_svg = os.path.join(out_dir, f"deg_vs_flux_scatter_{tissue}.svg")

    plt.tight_layout()
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.savefig(out_svg, dpi=300, bbox_inches="tight")
    plt.close()
    logging.info(f"Scatter saved to {out_pdf}")

    if not df_hi.empty:
        _save_legends(output_dir, tissue, size_mode, size_thresholds, size_values, size_default)

    used_pw = _collect_used_pathways(df_hi, df, tissue, rxn_lfc_thr, flux_log2_thr)
    out_used = os.path.join(output_dir, f"used_pathways_{tissue}.csv")
    used_pw.to_csv(out_used, index=False)

    return used_pw
