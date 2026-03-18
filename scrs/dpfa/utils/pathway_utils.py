"""
Utility functions for pathway manipulation

Provides functions for merging and filtering pathways in DataFrames.
"""

import pandas as pd
import logging
from typing import List, Dict, Optional


def apply_pathway_merging_to_df(df: pd.DataFrame,
                                 merging_dict: Dict[str, List[str]],
                                 pathway_column: str = "Pathways") -> pd.DataFrame:
    """
    Apply pathway merging to a DataFrame column

    """
    if not merging_dict or pathway_column not in df.columns:
        return df

    reverse_mapping = {}
    for new_name, old_names in merging_dict.items():
        for old_name in old_names:
            reverse_mapping[old_name] = new_name

    df = df.copy()
    df[pathway_column] = df[pathway_column].map(
        lambda x: reverse_mapping.get(x, x) if pd.notna(x) else x
    )

    return df


def apply_pathway_merging_to_list_column(df: pd.DataFrame,
                                         merging_dict: Dict[str, List[str]],
                                         list_column: str = "pathway_list") -> pd.DataFrame:
    """
    Apply pathway merging to a DataFrame column containing lists of pathways
    """
    if not merging_dict or list_column not in df.columns:
        return df

    # Create reverse mapping
    reverse_mapping = {}
    for new_name, old_names in merging_dict.items():
        for old_name in old_names:
            reverse_mapping[old_name] = new_name

    # Apply merging to lists
    df = df.copy()

    def merge_pathway_list(pw_list):
        if not isinstance(pw_list, list):
            return pw_list
        merged = [reverse_mapping.get(pw, pw) for pw in pw_list]
        # Remove duplicates while preserving order
        seen = set()
        result = []
        for pw in merged:
            if pw not in seen:
                seen.add(pw)
                result.append(pw)
        return result

    df[list_column] = df[list_column].apply(merge_pathway_list)

    return df


def filter_pathways_in_df(df: pd.DataFrame,
                          mode: str = "none",
                          whitelist: Optional[List[str]] = None,
                          pathway_column: str = "Pathways") -> pd.DataFrame:
    """
    Filter DataFrame rows based on pathway names
    """
    if mode == "none" or pathway_column not in df.columns:
        return df

    df = df.copy()

    if mode == "whitelist" and whitelist is not None:
        df = df[df[pathway_column].isin(whitelist)]
    elif mode == "blacklist" and blacklist is not None:
        df = df[~df[pathway_column].isin(blacklist)]
    else:
        logging.warning(f"Unknown filter mode: {mode}. No filtering applied.")

    return df


def aggregate_pathways_in_df(df: pd.DataFrame,
                             merging_dict: Dict[str, List[str]],
                             pathway_column: str = "Pathways",
                             value_columns: Optional[List[str]] = None,
                             aggregation: str = "sum") -> pd.DataFrame:
    """
    Merge pathways and aggregate numeric columns

    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with pathway and value columns
    merging_dict : dict
        Dictionary {new_name: [old_name1, old_name2, ...]}
    pathway_column : str
        Name of column containing pathway names
    value_columns : list, optional
        Columns to aggregate (default: all numeric columns)
    aggregation : str
        Aggregation method: "sum", "mean", "median", "max", "min"

    Returns:
    --------
    pd.DataFrame
        DataFrame with merged and aggregated pathways
    """
    if not merging_dict:
        return df

    df = apply_pathway_merging_to_df(df, merging_dict, pathway_column)

    if value_columns is None:
        value_columns = df.select_dtypes(include=['number']).columns.tolist()
        if pathway_column in value_columns:
            value_columns.remove(pathway_column)

    if not value_columns:
        return df.drop_duplicates(subset=[pathway_column], keep="first")

    agg_dict = {col: aggregation for col in value_columns}

    # Preserve non-numeric columns (take first value)
    non_numeric_cols = [col for col in df.columns
                       if col not in value_columns and col != pathway_column]
    for col in non_numeric_cols:
        agg_dict[col] = 'first'

    grouped = df.groupby(pathway_column, as_index=False).agg(agg_dict)

    return grouped


def get_merged_pathway_name(pathway: str,
                            merging_dict: Dict[str, List[str]]) -> str:
    """
    Get merged name for a single pathway
    """
    if not merging_dict:
        return pathway

    for new_name, old_names in merging_dict.items():
        if pathway in old_names:
            return new_name

    return pathway


