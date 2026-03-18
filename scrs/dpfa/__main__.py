"""
Entry point for DPFA (Differential Pathway Flux Analysis) module.

Usage:
    python -m scrs.dpfa
"""
import logging
import pandas as pd
import cobra

from scrs.config_loader import load_config
from scrs.dpfa.utils.pathway_database import PathwayDatabase
from .analysis import run_analysis


def main():
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

    flux_diff_thr = config.analysis.get("flux_diff_threshold")
    secretion_thr = config.analysis.get("secretion_threshold")
    drf_threshold = config.analysis.get("drf_threshold")

    flux_diff_plot_thr = config.get_flux_diff_threshold()
    metabolite_shortcuts = config.get_metabolite_name_shortcuts()
    metabolite_filter = config.get_metabolite_filter()
    drf_xlim_max = config.visualization.get("drf_xlim_max")

    filter_mode = config.pathway_filter.get("mode", "none")
    if filter_mode == "whitelist":
        whitelist = config.pathway_filter.get("whitelist", [])
        logging.info(f"  - Whitelist: {len(whitelist)} pathways")

    for cfg in tissue_configs:
        tissue = cfg["tissue"]
        logging.info(f"\n{'='*70}")
        logging.info(f"PROCESSING TISSUE: {tissue.upper()}")
        logging.info(f"{'='*70}")

        figure_sizes = cfg.get("figure_sizes", {})

        run_analysis(
            model=model,
            slow_model=cfg["slow_model"],
            fast_model=cfg["fast_model"],
            tissue=tissue,
            flux_diff_threshold=flux_diff_thr,
            output_dir=cfg.get("output_dir"),
            secretion_threshold=secretion_thr,
            global_pathway_colors=global_pathway_colors,
            pathway_merging=config.pathway_merging,
            pathway_filter_mode=filter_mode,
            pathways_filter=config.pathway_filter.get("whitelist") if filter_mode == "whitelist" else None,
            metabolite_shortcuts=metabolite_shortcuts,
            metabolite_filter=metabolite_filter,
            flux_diff_plot_threshold=flux_diff_plot_thr,
            drf_xlim_max=drf_xlim_max,
            drf_threshold=drf_threshold,
            mcpfa_log2fc_threshold=config.visualization.get('mcpfa_log2fc_threshold', 1.0),
            mcpfa_min_pathways=config.visualization.get('mcpfa_min_pathways', 2),
            figure_sizes=figure_sizes
        )

    logging.info("\n" + "="*70)
    logging.info("DPFA ANALYSIS COMPLETED")
    logging.info("="*70)


if __name__ == "__main__":
    main()
