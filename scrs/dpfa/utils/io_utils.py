import pandas as pd
from typing import Optional, Set


def load_deg_table(deg_csv: str,
                   padj_col: str = "padj",
                   pval_col: str = "pvalue",
                   lfc_col: str = "log2FoldChange",
                   id_col: str = "gene",
                   fdr_thr: Optional[float] = None) -> pd.DataFrame:
    """
    Load DEG table and optionally filter by padj < fdr_thr.
    Returns columns: gene, log2fc, pvalue, padj
    """
    deg = pd.read_csv(deg_csv)
    need = {id_col, lfc_col, pval_col, padj_col}
    missing = need - set(deg.columns)
    if missing:
        raise ValueError(f"DEG CSV lacks required columns: {missing}")

    deg = deg[[id_col, lfc_col, pval_col, padj_col]].rename(
        columns={id_col: "gene", lfc_col: "log2fc", pval_col: "pvalue", padj_col: "padj"}
    )

    if fdr_thr is not None:
        deg = deg[deg["padj"] < float(fdr_thr)].reset_index(drop=True)

    deg["log2fc"] = pd.to_numeric(deg["log2fc"], errors="coerce")
    deg["pvalue"] = pd.to_numeric(deg["pvalue"], errors="coerce")
    deg["padj"]   = pd.to_numeric(deg["padj"],   errors="coerce")

    return deg


def set_from_deg(deg: pd.DataFrame, direction: str) -> Set[str]:
    """
    Return set of genes by direction (up/down).
    Assumes deg is already filtered by padj if needed.
    """
    if direction == "up":
        return set(deg.loc[deg.log2fc > 0, "gene"])
    elif direction == "down":
        return set(deg.loc[deg.log2fc < 0, "gene"])
    else:
        raise ValueError("direction must be 'up' or 'down'")
