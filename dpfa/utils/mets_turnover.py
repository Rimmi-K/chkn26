import os
import re
from collections import defaultdict
from typing import Dict, Iterable, Optional, Tuple, Union, List

import numpy as np
import pandas as pd

try:
    from cobra import Model  
except Exception:
    Model = object

try:
    from .pathway_database import PathwayDatabase
except Exception:
    PathwayDatabase = object


def smart_capitalize(name: str) -> str:
    """Capitalize first letter in string"""
    m = re.search(r'([a-zA-Z])', str(name))
    if m:
        i = m.start()
        return name[:i] + name[i].upper() + name[i+1:]
    return str(name)


def shorten_names(name: str, custom_abbrev: Optional[Dict[str, str]] = None) -> str:
    """
    Shorten metabolite names using abbreviation dictionary

    Parameters:
    -----------
    name : str
        Metabolite name to shorten
    custom_abbrev : dict, optional
        Custom abbreviation dictionary from input_parameters.yaml

    Returns:
    --------
    str
        Shortened metabolite name
    """
    if custom_abbrev is None:
        custom_abbrev = {}

    s = str(name).strip().lower()
    for full, abbr in custom_abbrev.items():
        if full.lower() in s:
            s = s.replace(full.lower(), abbr)
    return smart_capitalize(s)


def compute_metabolite_turnover_by_subsystem(
    model: "Model",
    pathway_db: "PathwayDatabase",
    excluded_subsystems: Optional[Iterable[str]] = None,
    return_rxns: bool = False,
    store_contribs: bool = True,
) -> Union[
    Dict[Tuple[str, str], float],
    Tuple[Dict[Tuple[str, str], float], Dict[Tuple[str, str], List[Tuple[str, float]]]]
]:
    solution = model.optimize()
    fluxes = solution.fluxes

    flux_sum = defaultdict(float)
    rxn_details = defaultdict(list)  # (met_id, pathway) -> list[(rxn_id, contrib)]

    excluded_set = set(excluded_subsystems) if excluded_subsystems else None

    for met in model.metabolites:
        for rxn in met.reactions:
            pathways = pathway_db.get_pathways_for_reaction(rxn.id)

            if not pathways:
                continue

            v = float(fluxes.get(rxn.id, 0.0))
            contrib = abs(v)

            if contrib == 0.0:
                continue

            for pathway in pathways:
                if excluded_set and pathway in excluded_set:
                    continue

                key = (met.id, pathway)
                flux_sum[key] += contrib

                if return_rxns:
                    if store_contribs:
                        rxn_details[key].append((rxn.id, contrib))
                    else:
                        rxn_details[key].append((rxn.id, 1.0))  

    if return_rxns:
        return dict(flux_sum), dict(rxn_details)
    return dict(flux_sum)


def analyze_mets_turnover(
    model: "Model",
    model_fast: "Model",
    model_slow: "Model",
    pathway_db: "PathwayDatabase",
    output_dir: str,
    tissue: str = "unknown",
    exclude_subsystems: Optional[Iterable[str]] = ("Unknown",),
    log2fc_threshold: float = 0.5,
    min_flux: float = 1e-5,
    pathways_filter: Optional[list] = None,
    merge_compartments: bool = True,
    metabolite_shortcuts: Optional[Dict[str, str]] = None,
    pathway_merging: Optional[Dict[str, List[str]]] = None,
) -> pd.DataFrame:

    import os
    import logging

    if tissue == "liver":
        log2fc_threshold = 1.5

    os.makedirs(output_dir, exist_ok=True)

    flux_fast, rxns_fast = compute_metabolite_turnover_by_subsystem(
        model_fast,
        pathway_db,
        exclude_subsystems,
        return_rxns=True,
        store_contribs=True,   
    )
    flux_slow, rxns_slow = compute_metabolite_turnover_by_subsystem(
        model_slow,
        pathway_db,
        exclude_subsystems,
        return_rxns=True,
        store_contribs=True,
    )

    logging.info(
        f"[{tissue}] mcPFA raw entries: fast={len(flux_fast)}, slow={len(flux_slow)}"
    )

    all_keys = set(flux_fast.keys()) | set(flux_slow.keys())
    logging.info(f"[{tissue}] mcPFA shared entries: {len(all_keys)}")

    rows = []
    kept_min_flux = 0
    kept_log2fc = 0
    log2fc_vals = []

    def _extract_compartment(met_obj, met_id_str: str) -> str:
        comp = getattr(met_obj, "compartment", None) if met_obj is not None else None
        if comp:
            return str(comp)
        parts = met_id_str.rsplit("_", 1)
        if len(parts) == 2 and parts[1]:
            return parts[1]
        return "?"

    def _format_reactions_for_key(key: Tuple[str, str], top_n: Optional[int] = None) -> str:

        lst_f = rxns_fast.get(key, [])
        lst_s = rxns_slow.get(key, [])

        score = defaultdict(float)
        for rid, c in lst_f:
            if c > score[rid]:
                score[rid] = c
        for rid, c in lst_s:
            if c > score[rid]:
                score[rid] = c

        ordered = sorted(score.items(), key=lambda x: x[1], reverse=True)

        if top_n is not None:
            ordered = ordered[:top_n]

        return ";".join([rid for rid, _ in ordered])


    for met_id, subsystem in all_keys:
        fs = float(flux_fast.get((met_id, subsystem), 0.0))
        sl = float(flux_slow.get((met_id, subsystem), 0.0))


        if fs < min_flux and sl < min_flux:
            continue
        kept_min_flux += 1
        

        if fs < min_flux or sl < min_flux:
            log2fc = 6.0 if sl < fs else -6.0
        else:
            log2fc = float(np.log2((fs / (sl))))

        if abs(log2fc) <= log2fc_threshold:
            continue

        kept_log2fc += 1
        log2fc_vals.append(log2fc)

        met = None
        met_original_name = met_id
        try:
            met = model.metabolites.get_by_id(met_id)
            met_original_name = met.name
        except Exception:
            pass

        # Get compartment
        comp = _extract_compartment(met, met_id)

        # Original metabolite name with compartment
        original_met_label = f"{met_original_name} ({comp})" if comp != "?" else met_original_name

        # Shortened metabolite name (for grouping/visualization)
        shortened_name = shorten_names(met_original_name, metabolite_shortcuts)
        if merge_compartments:
            met_group = shortened_name
        else:
            met_group = f"{shortened_name} ({comp})"

        key = (met_id, subsystem)
        rxn_str = _format_reactions_for_key(key, top_n=None)

        rows.append(
            {
                "metabolite_id": met_id,
                "Metabolite": original_met_label,
                "Metabolite Group": met_group,
                "Metabolic Pathways": subsystem,
                "log2FC_fluxsum": log2fc,
                "Reactions": rxn_str,
            }
        )

    if kept_min_flux == 0:
        logging.warning(f"[{tissue}] mcPFA no entries passed min_flux={min_flux}")
    if kept_log2fc == 0:
        logging.warning(
            f"[{tissue}] mcPFA no entries passed log2fc_threshold={log2fc_threshold}"
        )
    if log2fc_vals:
        logging.info(
            f"[{tissue}] mcPFA log2fc stats: min={min(log2fc_vals):.3f}, "
            f"max={max(log2fc_vals):.3f}, mean={sum(log2fc_vals)/len(log2fc_vals):.3f}"
        )

    df = pd.DataFrame(rows)

    if pathway_merging and "Metabolic Pathways" in df.columns:
        reverse_mapping = {}
        for new_name, old_names in pathway_merging.items():
            for old_name in old_names:
                reverse_mapping[old_name] = new_name

        original_pathways = df["Metabolic Pathways"].nunique()
        df["Pathway Group"] = df["Metabolic Pathways"].map(
            lambda x: reverse_mapping.get(x, x) if pd.notna(x) else x
        )
        merged_pathways = df["Pathway Group"].nunique()

        logging.info(
            f"[{tissue}] Applied pathway merging to mcPFA: "
            f"{original_pathways} original pathways → {merged_pathways} pathway groups"
        )
    else:
        df["Pathway Group"] = df["Metabolic Pathways"]

    # Reorder columns
    cols = ["metabolite_id", "Metabolite", "Metabolite Group",
            "Metabolic Pathways", "Pathway Group",
            "log2FC_fluxsum", "Reactions"]
    df = df[cols]

    out_csv = os.path.join(output_dir, "fluxsum_log2fc_by_subsystem.csv")
    df.to_csv(out_csv, index=False)

    if df.empty:
        logging.warning(
            f"[{tissue}] mcPFA empty after filtering "
            f"(min_flux={min_flux}, log2fc_threshold={log2fc_threshold})."
        )
        return df

    logging.info(f"[{tissue}] mcPFA rows after filtering: {len(df)}")
    logging.info(f"[{tissue}] mcPFA table saved: {out_csv}")

    # Import here to avoid circular dependency
    from ..visualization import plot_fluxsum_log2fc_heatmap

    plot_fluxsum_log2fc_heatmap(
        df,
        tissue=tissue,
        output_dir=output_dir,
        pathways_filter=pathways_filter,
    )

    return df
