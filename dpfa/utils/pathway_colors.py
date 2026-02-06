"""Pathway color management for DPFA visualization."""

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib as mpl
from typing import List, Dict

mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

GLOBAL_PATHWAY_COLORS: Dict[str, str] = {}


def _generate_golden_angle_palette(n_colors: int = 60) -> List[str]:
    """Generate color palette using golden-angle HSV."""
    import colorsys

    golden_angle = 137.5
    colors = []
    brightness_levels = [0.92, 0.74, 0.58]
    saturation = 0.75

    hue_index = 0
    brightness_index = 0

    for i in range(n_colors):
        hue = (hue_index * golden_angle % 360) / 360.0
        brightness = brightness_levels[brightness_index % len(brightness_levels)]

        r, g, b = colorsys.hsv_to_rgb(hue, saturation, brightness)
        hex_color = f"#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}"
        colors.append(hex_color)

        hue_index += 1
        if (i + 1) % 3 == 0:
            brightness_index += 1

    return colors


def _initialize_color_pool() -> List[str]:
    """Create pool of contrasting colors."""
    colors = _generate_golden_angle_palette(n_colors=60)
    seen = set()
    dedup = []
    for h in colors:
        h_lower = h.lower()
        if h_lower not in seen:
            seen.add(h_lower)
            dedup.append(h)
    return dedup


_COLOR_POOL: List[str] = _initialize_color_pool()
_color_ptr = 0


def canonicalize_pathway(pw: str) -> str:
    """Normalize pathway name."""
    if not pw or not str(pw).strip():
        return "Unknown"
    return re.sub(r'\s+', ' ', str(pw).strip())


def get_color_for_pathway(pw: str) -> str:
    """Return color for pathway (auto-generated from golden angle palette)."""
    global _color_ptr

    pw = canonicalize_pathway(pw)

    if pw == "Unknown":
        return "#cccccc"

    if pw in GLOBAL_PATHWAY_COLORS:
        return GLOBAL_PATHWAY_COLORS[pw]

    new_color = _COLOR_POOL[_color_ptr % len(_COLOR_POOL)]
    GLOBAL_PATHWAY_COLORS[pw] = new_color
    _color_ptr += 1

    return new_color


def palette_for_df(df: pd.DataFrame) -> Dict[str, str]:
    """Create color palette for DataFrame with 'Pathways' column."""
    cats = sorted([canonicalize_pathway(x) for x in df["Pathways"].dropna().unique()])
    return {pw: get_color_for_pathway(pw) for pw in cats}


def save_legend_from_df(used_df: pd.DataFrame,
                        output_dir: str,
                        fname: str = "global_legend_pathways.pdf",
                        ncol: int = 2) -> None:
    """Draw legend from DataFrame with 'Pathways' and 'color_hex' columns."""
    os.makedirs(os.path.join(output_dir, "plots"), exist_ok=True)

    if used_df is None or used_df.empty:
        return

    df = used_df.copy()
    df["Pathways"] = df["Pathways"].map(canonicalize_pathway)
    df = df.dropna(subset=["Pathways"]).copy()

    if "color_hex" not in df.columns:
        df["color_hex"] = None

    df["color_hex"] = df.apply(
        lambda row: row["color_hex"] if pd.notna(row["color_hex"]) and row["color_hex"].strip()
                   else get_color_for_pathway(row["Pathways"]),
        axis=1
    )

    df = df.drop_duplicates(subset=["Pathways"], keep="first")

    items = sorted(
        df[["Pathways", "color_hex"]].values.tolist(),
        key=lambda x: x[0].lower()
    )

    if not items:
        return

    handles = [Line2D([0],[0], marker='o', linestyle='', color=col,
                     label=lab, markersize=12, markeredgewidth=0)
               for lab, col in items]

    fig, ax = plt.subplots(figsize=(8, max(1.5, 0.35*len(handles)/ncol)))
    ax.legend(handles=handles, loc="center left", frameon=True,
             title="Pathways", ncol=ncol, fontsize=11, title_fontsize=13)
    ax.axis("off")

    out_path = os.path.join(output_dir, "plots", fname)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    pd.DataFrame(items, columns=["Pathways", "color_hex"]).to_csv(
        os.path.join(output_dir, "pathway_colors_used.csv"), index=False
    )
