import os
import logging
import pandas as pd
import cobra

from dpfa.utils.pathway_colors import save_legend_from_df
from dpfa.utils.pathway_database import PathwayDatabase
from .visualization import (
    plot_regulation_counts,
    analyze_pathway_flux_difference,
    create_rdpfa_summary
)
from dpfa.utils.flux_utils import (
    process_models,
    ratio_fluxes,
    create_DRF,
    get_reaction_formula
)
from dpfa.utils.mets_turnover import analyze_mets_turnover


def run_analysis(model, slow_model, fast_model, tissue,
                 pathways_csv=None,
                 output_dir="results_analysis",
                 pathways_filter=None,
                 flux_diff_threshold=0.0,
                 global_pathway_colors: dict = None,
                 secretion_threshold: float = 1e-6,
                 pathway_merging: dict = None,
                 pathway_filter_mode: str = "none",
                 metabolite_filter: dict = None,
                 metabolite_shortcuts: dict = None,
                 flux_diff_plot_threshold: float = 0.0,
                 drf_xlim_max: int = None,
                 drf_threshold: float = 0.3,
                 mcpfa_log2fc_threshold: float = 1.5,
                 mcpfa_min_pathways: int = 2):
    """
    Performs differential pathway flux analysis (DPFA) of metabolic models.

    Includes:
    - rDPFA (reaction-centric differential pathway flux analysis)
    - mcDPFA (metabolite-centric differential pathway flux analysis)
    - DRF histogram plots
    - Metabolite turnover heatmaps
   """
    os.makedirs(output_dir, exist_ok=True)

    slow_flux_path, fast_flux_path = process_models(slow_model, fast_model, tissue, output_dir)
    ratio_path = os.path.join(output_dir, f'ratio_{tissue}.csv')
    ratio_df = ratio_fluxes(slow_flux_path, fast_flux_path, ratio_path)
    drf_df = create_DRF(ratio_df, drf_threshold)
    if "reaction_id" in drf_df.columns:
        drf_df["reaction_id"] = drf_df["reaction_id"].astype(str).str.strip()

    model_slow = None
    model_fast = None
    try:
        model_slow = cobra.io.load_json_model(slow_model)
        model_fast = cobra.io.load_json_model(fast_model)
    except Exception as e:
        logging.exception(f"Model loading failed: {e}")


    drf_df["Formula"] = drf_df["reaction_id"].apply(lambda rxn_id: get_reaction_formula(model, rxn_id))


    pathway_db = PathwayDatabase('data/models/subsystem_matrix.csv')
    merged_df = drf_df.copy()

    def _get_pathways_for_rxn(rxn_id):
        """Get list of pathways for reaction from KEGG database"""
        pathways = pathway_db.get_pathways_for_reaction(rxn_id)
        return pathways if pathways else ["Unknown"]

    # Store original pathways
    merged_df["Pathways_list"] = merged_df["reaction_id"].apply(_get_pathways_for_rxn)
    merged_df["Pathways"] = merged_df["Pathways_list"].apply(lambda lst: "; ".join(lst))

    multi_pathway_rxns = merged_df[merged_df["Pathways_list"].apply(len) > 1]
    if len(multi_pathway_rxns) > 0:
        sample = multi_pathway_rxns.iloc[0]
        logging.info(
            f"[{tissue}] multi-pathway sample: "
            f"{sample['reaction_id']} -> {sample['Pathways_list']}"
        )

    for rid in ("GHMT2rm", "GHMT2r"):
        pathways = pathway_db.get_pathways_for_reaction(rid)
        if pathways:
            logging.info(f"[{tissue}] {rid} pathways: {pathways}")

    # Create pathway groups for visualization (keep original pathways intact)
    if pathway_merging:
        from dpfa.utils.pathway_utils import apply_pathway_merging_to_list_column

        if "Pathways_list" in merged_df.columns:
            # Create a copy for grouping
            merged_df_grouped = apply_pathway_merging_to_list_column(
                merged_df.copy(), pathway_merging, "Pathways_list"
            )
            merged_df["Subsystems_list"] = merged_df_grouped["Pathways_list"]
            merged_df["Subsystems"] = merged_df["Subsystems_list"].apply(lambda lst: "; ".join(lst))
        else:
            merged_df["Subsystems_list"] = merged_df["Pathways_list"]
            merged_df["Subsystems"] = merged_df["Pathways"]

        logging.info(f"[{tissue}] Applied pathway merging to analysis data")
        logging.info(f"[{tissue}] Original pathways: {merged_df['Pathways'].nunique()}, Subsystems: {merged_df['Subsystems'].nunique()}")
    else:
        # No merging - subsystems are same as original pathways
        merged_df["Subsystems_list"] = merged_df["Pathways_list"]
        merged_df["Subsystems"] = merged_df["Pathways"]

    logging.info(f"\n[{tissue}] Pathways in merged_df: {merged_df['Pathways'].nunique()}")
    logging.info(f"[{tissue}] DRF categories: {merged_df['DRF_category'].value_counts().to_dict()}")

    # Create unified rDPFA summary file (replaces separate drf, pathway_flux_values, pathwaygroup_flux_values files)
    rdpfa_summary = create_rdpfa_summary(merged_df, output_dir, tissue)

    filter_pathways_list = None

    if pathway_filter_mode == "whitelist" and pathways_filter:
        filter_pathways_list = pathways_filter
        logging.info(f"[{tissue}] Pathway filter: whitelist ({len(filter_pathways_list)} pathways)")
    elif pathway_filter_mode == "blacklist" and pathway_blacklist:
        all_pathways = set()
        if "Pathways" in merged_df.columns:
            all_pathways.update(merged_df["Pathways"].dropna().unique())
        filter_pathways_list = [p for p in all_pathways if p not in pathway_blacklist]
        logging.info(f"[{tissue}] Pathway filter: blacklist ({len(pathway_blacklist)} excluded, {len(filter_pathways_list)} remaining)")
    else:
        logging.info(f"[{tissue}] Pathway filter: none (all pathways)")

    try:
        if model_slow is None or model_fast is None:
            raise RuntimeError("Models are not loaded")
        analyze_mets_turnover(
            model=model,
            model_fast=model_fast,
            model_slow=model_slow,
            pathway_db=pathway_db,
            output_dir=output_dir,
            tissue=tissue,
            exclude_subsystems=('Unknown',),
            min_flux=1e-6,
            pathways_filter=filter_pathways_list,
            merge_compartments=True,
            metabolite_shortcuts=metabolite_shortcuts,
            pathway_merging=pathway_merging,
            metabolite_filter=metabolite_filter,
            log2fc_threshold=mcpfa_log2fc_threshold,
            min_pathways_with_change=mcpfa_min_pathways,
        )
    except Exception as e:
        logging.exception(f"Metabolite turnover analysis failed: {e}")

    pathway_flux_diff = analyze_pathway_flux_difference(
        merged_df, output_dir, tissue,
        threshold=flux_diff_threshold if flux_diff_threshold is not None else 0.0
    )

    plot_regulation_counts(
        merged_df,
        output_dir,
        tissue,
        pathways_filter=filter_pathways_list,
        pathway_flux_diff=pathway_flux_diff,
        flux_diff_threshold=flux_diff_plot_threshold,
        xlim=(0, drf_xlim_max) if drf_xlim_max else None
    )

    logging.info(f"[{tissue}] DPFA analysis completed")

    # Return pathway colors used in analysis for global legend
    return {}

