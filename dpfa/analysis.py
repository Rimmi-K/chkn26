import os
import logging
import pandas as pd
import cobra

from .utils.pathway_colors import save_legend_from_df
from .utils.pathway_database import PathwayDatabase
from .utils.gpr_utils import (
    build_gpr_clauses,
    reaction_log2fc_via_gpr,
    reaction_pvalues_via_gpr,
)
from .utils.io_utils import load_deg_table
from .visualization import (
    plot_regulation_counts,
    analyze_pathway_flux_difference
)
from .scatter_plot import make_scatter_deg_vs_flux
from .utils.flux_utils import (
    process_models,
    ratio_fluxes,
    create_DRF
)
from .utils.mets_turnover import analyze_mets_turnover


def run_analysis(model, slow_model, fast_model, tissue,
                 pathways_csv=None,
                 output_dir="results_analysis",
                 pathways_filter=None,
                 flux_diff_threshold=0.0,
                 deg_csv=None,
                 gene_id_map=None,
                 fdr_thr=0.05,
                 scatter_rxn_lfc_thr: float = 0.5,
                 scatter_flux_thr: float = 0.5,
                 scatter_require_sig: bool = False,
                 gpr_missing_strategy: str = "impute",
                 gpr_impute_value: float = 0.0,
                 global_pathway_colors: dict = None,
                 secretion_threshold: float = 1e-6,
                 pathway_merging: dict = None,
                 pathway_filter_mode: str = "none",
                 metabolite_shortcuts: dict = None,
                 flux_diff_plot_threshold: float = 0.0,
                 drf_xlim_max: int = None,
                 drf_threshold: float = 0.3):
    """
    Performs comprehensive differential pathway flux analysis (DPFA) of metabolic models and DEG data.
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

    pathway_db = PathwayDatabase('data/models/subsystem_matrix.csv')
    merged_df = drf_df.copy()

    def _get_pathways_for_rxn(rxn_id):
        """Get list of pathways for reaction from KEGG database"""
        pathways = pathway_db.get_pathways_for_reaction(rxn_id)
        return pathways if pathways else ["Unknown"]

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

    if pathway_merging:
        from .utils.pathway_utils import apply_pathway_merging_to_list_column

        if "Pathways_list" in merged_df.columns:
            merged_df = apply_pathway_merging_to_list_column(
                merged_df, pathway_merging, "Pathways_list"
            )
            merged_df["Pathways"] = merged_df["Pathways_list"].apply(lambda lst: "; ".join(lst))

        logging.info(f"[{tissue}] Applied pathway merging to analysis data")
        logging.info(f"[{tissue}] Merged pathways count: {merged_df['Pathways'].nunique()}")

    merged_df.to_csv(os.path.join(output_dir, f'drf_{tissue}.csv'), index=False)

    logging.info(f"\n[{tissue}] Pathways in merged_df: {merged_df['Pathways'].nunique()}")
    logging.info(f"[{tissue}] DRF categories: {merged_df['DRF_category'].value_counts().to_dict()}")

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
            min_flux=1e-3,
            pathways_filter=filter_pathways_list,
            merge_compartments=True,
            metabolite_shortcuts=metabolite_shortcuts,
            pathway_merging=pathway_merging
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

    used_df = pd.DataFrame(columns=["Pathways", "n_highlighted", "color_hex", "tissue"])

    if deg_csv:
        model_for_gpr = cobra.io.load_json_model(fast_model)
        clauses = build_gpr_clauses(model_for_gpr)
        deg = load_deg_table(deg_csv)

        if 'padj' not in deg.columns:
            if 'pvalue' in deg.columns:
                logging.warning(f"[{tissue}] DEG table doesn't have 'padj' column. "
                              "Will use 'pvalue' but recommend running proper BH correction first!")
                from statsmodels.stats.multitest import multipletests
                _, padj_values, _, _ = multipletests(deg['pvalue'], method='fdr_bh')
                deg['padj'] = padj_values
            else:
                raise ValueError("DEG table must have either 'padj' or 'pvalue' column")

        if gene_id_map is not None:
            deg["gene"] = deg["gene"].map(gene_id_map).fillna(deg["gene"])

        try:
            rxn_lfc_map = reaction_log2fc_via_gpr(
                clauses=clauses,
                deg_df=deg,
                missing_strategy=gpr_missing_strategy,
                impute_value=gpr_impute_value
            )

            rxn_p_map = reaction_pvalues_via_gpr(
                clauses=clauses,
                deg_df=deg,
                p_col="padj",
                require_all_genes=True
            )

            rxn_q_map = rxn_p_map
            rxn_sig_map = {rid: (q < float(fdr_thr)) for rid, q in rxn_q_map.items()}

            md_for_paths = model_for_gpr if "Pathways" not in merged_df.columns else None

            used_df = make_scatter_deg_vs_flux(
                merged_df=merged_df,
                rxn_lfc_map=rxn_lfc_map,
                rxn_sig_map=rxn_sig_map if scatter_require_sig else {rid: True for rid in rxn_lfc_map},
                tissue=tissue,
                output_dir=output_dir,
                model_for_paths=md_for_paths,
                rxn_lfc_thr=scatter_rxn_lfc_thr,
                flux_log2_thr=scatter_flux_thr,
                pathway_merging=pathway_merging,
                pathway_filter=filter_pathways_list,
                model_for_gpr_rules=model_for_gpr,
                deg_df_for_labels=deg,
                rxn_p_map=rxn_p_map,
                padj_rxn_map=rxn_q_map,
                size_mode="discrete",
                size_thresholds=[0.001, 0.01, 0.05, 0.25],
                size_values=[175, 150, 100, 50],
                size_default=25
            )

            logging.info(f"[{tissue}] Reaction statistics:")
            logging.info(f"  - Total reactions with log2FC: {len(rxn_lfc_map)}")
            logging.info(f"  - Total reactions with p-value: {len(rxn_p_map)}")
            logging.info(f"  - Significant reactions (padj<{fdr_thr}): {sum(rxn_sig_map.values())}")

        except Exception as e:
            logging.exception(f"Scatter DEG vs flux failed: {e}")

    if isinstance(used_df, pd.DataFrame) and len(used_df) > 0:
        return dict(zip(used_df["Pathways"], used_df["color_hex"]))
    else:
        return {}

if __name__ == "__main__":
    from .config_loader import load_config
    from .utils.pathway_utils import apply_pathway_merging_to_df

    logging.info("="*70)
    logging.info("STARTING METABOLIC FLUX ANALYSIS")
    logging.info("="*70)

    config_path = "input_parameters.yaml"
    try:
        config = load_config(config_path)
        config.validate()
        logging.info(f"Loaded configuration from {config_path}")
    except FileNotFoundError:
        logging.error(f"ERROR: Configuration file not found {config_path}")
        logging.error("Please create input_parameters.yaml or use default settings")
        raise
    except Exception as e:
        logging.error(f"ERROR: Configuration file - {e}")
        raise

    base_model_path = config.get_base_model_path()
    model = cobra.io.read_sbml_model(base_model_path)

    pathway_db_path = config.get_pathway_database_path()
    pathway_db = PathwayDatabase(pathway_db_path)
    all_pathways = pathway_db.get_all_pathways(exclude_empty=True)

    global_pathway_colors = pathway_db.assign_global_colors(pathways=all_pathways)

    tissue_configs = config.get_tissue_configs()
    logging.info(f"Processing {len(tissue_configs)} tissues: {[t['tissue'] for t in tissue_configs]}")

    fdr_thr = config.get_fdr_threshold()
    flux_diff_thr = config.analysis.get("flux_diff_threshold")
    gpr_strategy = config.analysis.get("gpr_missing_strategy")
    gpr_impute = config.analysis.get("gpr_impute_value")
    secretion_thr = config.analysis.get("secretion_threshold")
    drf_threshold = config.analysis.get("drf_threshold")

    flux_diff_plot_thr = config.get_flux_diff_threshold()
    metabolite_shortcuts = config.get_metabolite_name_shortcuts()
    drf_xlim_max = config.visualization.get("drf_xlim_max")
    scatter_params = config.scatter

    filter_mode = config.pathway_filter.get("mode", "none")
    if filter_mode == "whitelist":
        whitelist = config.pathway_filter.get("whitelist", [])
        logging.info(f"  - Whitelist: {len(whitelist)} pathways")

    used_all = pd.DataFrame(columns=["Pathways", "color_hex"])

    for cfg in tissue_configs:
        tissue = cfg["tissue"]
        logging.info(f"\n{'='*70}")
        logging.info(f"PROCESSING TISSUE: {tissue.upper()}")
        logging.info(f"{'='*70}")

        used_colors = run_analysis(
            model=model,
            slow_model=cfg["slow_model"],
            fast_model=cfg["fast_model"],
            tissue=tissue,
            fdr_thr=fdr_thr,
            flux_diff_threshold=flux_diff_thr,
            output_dir=cfg.get("output_dir"),
            deg_csv=cfg.get("deg_csv"),
            gpr_missing_strategy=gpr_strategy,
            gpr_impute_value=gpr_impute,
            secretion_threshold=secretion_thr,
            global_pathway_colors=global_pathway_colors,
            pathway_merging=config.pathway_merging,
            pathway_filter_mode=filter_mode,
            pathways_filter=config.pathway_filter.get("whitelist") if filter_mode == "whitelist" else None,
            scatter_rxn_lfc_thr=scatter_params['rxn_lfc_threshold'],
            scatter_flux_thr=scatter_params['flux_threshold'],
            scatter_require_sig=scatter_params.get('require_significance', False),
            metabolite_shortcuts=metabolite_shortcuts,
            flux_diff_plot_threshold=flux_diff_plot_thr,
            drf_xlim_max=drf_xlim_max,
            drf_threshold=drf_threshold
        )

        if isinstance(used_colors, dict) and used_colors:
            part = (pd.DataFrame.from_dict(used_colors, orient="index", columns=["color_hex"])
                    .reset_index().rename(columns={"index": "Pathways"}))
        elif isinstance(used_colors, pd.DataFrame) and not used_colors.empty:
            cols = [c for c in used_colors.columns if c in ("Pathways", "color_hex")]
            part = used_colors[cols].copy()
        else:
            continue

        if config.pathway_merging:
            part = apply_pathway_merging_to_df(part, config.pathway_merging, "Pathways")

        part = part.dropna(subset=["Pathways"]).drop_duplicates(subset=["Pathways", "color_hex"])
        used_all = pd.concat([used_all, part], ignore_index=True)

    base_output_dir = config.get_base_output_dir()

    if not used_all.empty:
        used_all = (used_all
                    .dropna(subset=["Pathways"])
                    .drop_duplicates(subset=["Pathways"], keep="first"))

        legend_fname = config.visualization.get("legend_filename", "global_legend_pathways.pdf")
        legend_ncol = config.visualization.get("legend_ncol", 2)

        save_legend_from_df(used_all, output_dir=base_output_dir,
                          fname=legend_fname, ncol=legend_ncol)
        logging.info(f"Saved global legend: {base_output_dir}/{legend_fname}")
    else:
        logging.info("[main] No pathways to draw in legend")


