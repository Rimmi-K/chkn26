"""
Entry point for TFCA (Transcript-Flux Concordance Analysis) module.

Usage:
    python -m scrs.tfca
"""
import logging
import pandas as pd
import cobra
from statsmodels.stats.multitest import multipletests

from scrs.config_loader import load_config
from scrs.dpfa.utils.pathway_database import PathwayDatabase
from scrs.dpfa.utils.gpr_utils import (
    build_gpr_clauses,
    reaction_log2fc_via_gpr,
    reaction_pvalues_via_gpr,
)
from scrs.dpfa.utils.io_utils import load_deg_table
from scrs.dpfa.utils.flux_utils import (
    process_models,
    ratio_fluxes,
    create_DRF,
)
from .analysis import make_scatter_deg_vs_flux


def main():
    logging.info("="*70)
    logging.info("STARTING TRANSCRIPT-FLUX CONCORDANCE ANALYSIS (TFCA)")
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

    pathway_db_path = config.get_pathway_database_path()
    pathway_db = PathwayDatabase(pathway_db_path)

    tissue_configs = config.get_tissue_configs()
    logging.info(f"Processing {len(tissue_configs)} tissues: {[t['tissue'] for t in tissue_configs]}")

    fdr_thr = config.get_fdr_threshold()
    gpr_strategy = config.analysis.get("gpr_missing_strategy")
    gpr_impute = config.analysis.get("gpr_impute_value")
    drf_threshold = config.analysis.get("drf_threshold")

    scatter_params = config.scatter
    filter_mode = config.pathway_filter.get("mode", "none")

    for cfg in tissue_configs:
        tissue = cfg["tissue"]
        logging.info(f"\n{'='*70}")
        logging.info(f"PROCESSING TISSUE: {tissue.upper()}")
        logging.info(f"{'='*70}")

        slow_model = cfg["slow_model"]
        fast_model = cfg["fast_model"]
        deg_csv = cfg.get("deg_csv")
        output_dir = cfg.get("output_dir")
        figure_sizes = cfg.get("figure_sizes", {})
        scatter_axis_limits = cfg.get("scatter_axis_limits", {})

        if not deg_csv:
            logging.warning(f"[{tissue}] No DEG file specified, skipping TFCA")
            continue

        # Process flux models
        slow_flux_path, fast_flux_path = process_models(slow_model, fast_model, tissue, output_dir)
        ratio_df = ratio_fluxes(slow_flux_path, fast_flux_path,
                               f"{output_dir}/ratio_{tissue}.csv")
        drf_df = create_DRF(ratio_df, drf_threshold)

        if "reaction_id" in drf_df.columns:
            drf_df["reaction_id"] = drf_df["reaction_id"].astype(str).str.strip()

        # Add pathways to drf_df
        def _get_pathways_for_rxn(rxn_id):
            pathways = pathway_db.get_pathways_for_reaction(rxn_id)
            return pathways if pathways else ["Unknown"]

        drf_df["Pathways_list"] = drf_df["reaction_id"].apply(_get_pathways_for_rxn)
        drf_df["Pathways"] = drf_df["Pathways_list"].apply(lambda lst: "; ".join(lst))

        # Apply pathway merging if configured
        if config.pathway_merging:
            from scrs.dpfa.utils.pathway_utils import apply_pathway_merging_to_list_column

            if "Pathways_list" in drf_df.columns:
                drf_df_grouped = apply_pathway_merging_to_list_column(
                    drf_df.copy(), config.pathway_merging, "Pathways_list"
                )
                drf_df["Subsystems_list"] = drf_df_grouped["Pathways_list"]
                drf_df["Subsystems"] = drf_df["Subsystems_list"].apply(lambda lst: "; ".join(lst))
        else:
            drf_df["Subsystems_list"] = drf_df["Pathways_list"]
            drf_df["Subsystems"] = drf_df["Pathways"]

        # Prepare pathway filter
        filter_pathways_list = None
        if filter_mode == "whitelist":
            filter_pathways_list = config.pathway_filter.get("whitelist", [])

        # Load and process DEG data
        model_for_gpr = cobra.io.load_json_model(fast_model)
        clauses = build_gpr_clauses(model_for_gpr)
        deg = load_deg_table(deg_csv)

        if 'padj' not in deg.columns:
            if 'pvalue' in deg.columns:
                logging.warning(f"[{tissue}] DEG table doesn't have 'padj' column. "
                              "Applying BH correction to pvalues...")
                _, padj_values, _, _ = multipletests(deg['pvalue'], method='fdr_bh')
                deg['padj'] = padj_values
            else:
                raise ValueError("DEG table must have either 'padj' or 'pvalue' column")

        # Calculate reaction-level log2FC and p-values from GPR
        rxn_lfc_map = reaction_log2fc_via_gpr(
            clauses=clauses,
            deg_df=deg,
            missing_strategy=gpr_strategy,
            impute_value=gpr_impute
        )

        rxn_p_map = reaction_pvalues_via_gpr(
            clauses=clauses,
            deg_df=deg,
            p_col="padj",
            require_all_genes=True
        )

        rxn_q_map = rxn_p_map
        rxn_sig_map = {rid: (q < float(fdr_thr)) for rid, q in rxn_q_map.items()}

        # Run TFCA scatter plot analysis
        logging.info(f"[{tissue}] Running transcript-flux concordance analysis...")

        used_df = make_scatter_deg_vs_flux(
            merged_df=drf_df,
            rxn_lfc_map=rxn_lfc_map,
            rxn_sig_map=rxn_sig_map if scatter_params.get('require_significance', False) else {rid: True for rid in rxn_lfc_map},
            tissue=tissue,
            output_dir=output_dir,
            model_for_paths=None,  # Pathways already in merged_df
            rxn_lfc_thr=scatter_params['rxn_lfc_threshold'],
            flux_log2_thr=scatter_params['flux_threshold'],
            pathway_merging=config.pathway_merging,
            pathway_filter=filter_pathways_list,
            model_for_gpr_rules=model_for_gpr,
            deg_df_for_labels=deg,
            rxn_p_map=rxn_p_map,
            padj_rxn_map=rxn_q_map,
            size_mode=scatter_params.get('size_mode', 'discrete'),
            size_thresholds=scatter_params.get('size_thresholds', [0.001, 0.01, 0.05, 0.25]),
            size_values=scatter_params.get('size_values', [175, 150, 100, 50]),
            size_default=scatter_params.get('size_default', 25),
            figsize=figure_sizes.get('scatter'),
            xlim=scatter_axis_limits.get('xlim'),
            ylim=scatter_axis_limits.get('ylim')
        )

        logging.info(f"[{tissue}] TFCA analysis completed")
        logging.info(f"[{tissue}] Reaction statistics:")
        logging.info(f"  - Total reactions with log2FC: {len(rxn_lfc_map)}")
        logging.info(f"  - Total reactions with p-value: {len(rxn_p_map)}")
        logging.info(f"  - Significant reactions (padj<{fdr_thr}): {sum(rxn_sig_map.values())}")

    logging.info("="*70)
    logging.info("TFCA ANALYSIS COMPLETED")
    logging.info("="*70)


if __name__ == "__main__":
    main()
