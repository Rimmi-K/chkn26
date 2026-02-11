"""
Configuration loader for metabolic flux analysis

Loads parameters from input_parameters.yaml and provides utilities
for pathway merging and filtering.
"""

import os
import logging
from typing import List, Dict, Any, Optional

try:
    import yaml
    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False
    logging.warning("PyYAML not available. Install with: pip install pyyaml")


class AnalysisConfig:
    """Configuration manager for metabolic flux analysis"""

    def __init__(self, config_path: str = "input_parameters.yaml"):
        """
        Load configuration from YAML file

        Parameters:
        -----------
        config_path : str
            Path to configuration YAML file
        """
        self.config_path = config_path
        self.config = self._load_config()

        # Extract main sections
        self.paths = self.config.get("paths", {})
        self.tissues = self.config.get("tissues", [])
        self.analysis = self.config.get("analysis", {})
        self.scatter = self.config.get("scatter", {})
        self.pathway_merging = self.config.get("pathway_merging", {})
        self.pathway_filter = self.config.get("pathway_filter", {})
        self.visualization = self.config.get("visualization", {})
        self.metabolite_name_shortcuts = self.config.get("metabolite_name_shortcuts", {})

    def _load_config(self) -> Dict[str, Any]:
        """Load YAML configuration file"""
        if not YAML_AVAILABLE:
            raise ImportError("PyYAML is required. Install with: pip install pyyaml")

        if not os.path.exists(self.config_path):
            raise FileNotFoundError(f"Config file not found: {self.config_path}")

        with open(self.config_path, 'r') as f:
            config = yaml.safe_load(f)

        logging.info(f"Loaded configuration from {self.config_path}")
        return config

    def apply_pathway_merging(self, pathways_list: List[str]) -> List[str]:
        """
        Apply pathway merging to a list of pathways

        Parameters:
        -----------
        pathways_list : list
            List of pathway names

        Returns:
        --------
        list
            List of pathways after merging
        """
        if not self.pathway_merging:
            return pathways_list

        # Create reverse mapping: old_name -> new_name
        reverse_mapping = {}
        for new_name, old_names in self.pathway_merging.items():
            for old_name in old_names:
                reverse_mapping[old_name] = new_name

        # Apply merging
        merged_pathways = []
        for pathway in pathways_list:
            merged_name = reverse_mapping.get(pathway, pathway)
            merged_pathways.append(merged_name)

        return merged_pathways

    def get_merged_pathway_name(self, pathway: str) -> str:
        """
        Get merged name for a single pathway

        Parameters:
        -----------
        pathway : str
            Pathway name

        Returns:
        --------
        str
            Merged name (or original if not found in merging dict)
        """
        if not self.pathway_merging:
            return pathway

        for new_name, old_names in self.pathway_merging.items():
            if pathway in old_names:
                return new_name

        return pathway

    def filter_pathways(self, pathways_list: List[str]) -> List[str]:
        """
        Filter pathways according to configuration

        Parameters:
        -----------
        pathways_list : list
            List of pathway names

        Returns:
        --------
        list
            Filtered list of pathways
        """
        mode = self.pathway_filter.get("mode", "none")

        if mode == "none":
            return pathways_list
        elif mode == "whitelist":
            whitelist = self.pathway_filter.get("whitelist", [])
            return [p for p in pathways_list if p in whitelist]
        elif mode == "blacklist":
            blacklist = self.pathway_filter.get("blacklist", [])
            return [p for p in pathways_list if p not in blacklist]
        else:
            logging.warning(f"Unknown filter mode: {mode}. No filtering applied.")
            return pathways_list

    def should_include_pathway(self, pathway: str) -> bool:
        """
        Check if pathway should be included in analysis

        Parameters:
        -----------
        pathway : str
            Pathway name

        Returns:
        --------
        bool
            True if pathway should be included
        """
        mode = self.pathway_filter.get("mode", "none")

        if mode == "none":
            return True
        elif mode == "whitelist":
            whitelist = self.pathway_filter.get("whitelist", [])
            return pathway in whitelist
        elif mode == "blacklist":
            blacklist = self.pathway_filter.get("blacklist", [])
            return pathway not in blacklist
        else:
            return True

    def validate(self) -> bool:
        """
        Validate configuration

        Returns:
        --------
        bool
            True if configuration is valid

        Raises:
        -------
        ValueError
            If configuration is invalid
        """
        errors = []

        # Check paths exist
        base_model = self.paths.get("base_model")
        if base_model and not os.path.exists(base_model):
            errors.append(f"Base model not found: {base_model}")

        pathway_db = self.paths.get("pathway_database")
        if pathway_db and not os.path.exists(pathway_db):
            errors.append(f"Pathway database not found: {pathway_db}")

        # Check tissue configs
        for cfg in self.tissues:
            tissue = cfg.get("tissue", "unknown")

            slow_model = cfg.get("slow_model")
            if slow_model and not os.path.exists(slow_model):
                errors.append(f"Slow model not found: {slow_model} ({tissue})")

            fast_model = cfg.get("fast_model")
            if fast_model and not os.path.exists(fast_model):
                errors.append(f"Fast model not found: {fast_model} ({tissue})")

            deg_csv = cfg.get("deg_csv")
            if deg_csv and not os.path.exists(deg_csv):
                errors.append(f"DEG file not found: {deg_csv} ({tissue})")

        # Check pathway filter mode
        filter_mode = self.pathway_filter.get("mode", "none")
        if filter_mode not in ["whitelist", "blacklist", "none"]:
            errors.append(f"Invalid pathway_filter mode: {filter_mode}")

        if errors:
            raise ValueError("Configuration validation failed:\n" +
                           "\n".join(f"  - {e}" for e in errors))

        logging.info("Configuration validation passed")
        return True

    def get_base_model_path(self) -> str:
        """Get path to base model"""
        return self.paths.get("base_model", "")

    def get_pathway_database_path(self) -> str:
        """Get path to pathway database"""
        return self.paths.get("pathway_database", "")

    def get_base_output_dir(self) -> str:
        """Get base output directory"""
        return self.paths.get("base_output_dir", "results_riptide")

    def get_tissue_configs(self) -> List[Dict[str, str]]:
        """
        Get list of tissue configurations with processed output directories.

        If output_dir is not specified for a tissue, it will be constructed
        from base_output_dir and tissue name.
        """
        base_output_dir = self.paths.get("base_output_dir", "results/fluxomics")

        processed_tissues = []
        for tissue_cfg in self.tissues:
            cfg = tissue_cfg.copy()

            # If output_dir not specified, construct it from base_output_dir + tissue name
            if "output_dir" not in cfg or not cfg["output_dir"]:
                tissue_name = cfg.get("tissue", "unknown")
                cfg["output_dir"] = os.path.join(base_output_dir, tissue_name)

            processed_tissues.append(cfg)

        return processed_tissues

    def get_fdr_threshold(self) -> float:
        """Get FDR threshold"""
        return self.analysis.get("fdr_threshold", 0.05)

    def get_scatter_params(self) -> Dict[str, Any]:
        """Get scatter plot parameters"""
        return self.scatter

    def get_metabolite_name_shortcuts(self) -> Dict[str, str]:
        """Get metabolite name shortcuts dictionary"""
        return self.metabolite_name_shortcuts

    def get_flux_diff_threshold(self) -> float:
        """Get flux difference threshold for DRF plots"""
        return self.visualization.get("flux_diff_threshold", 0.0)

    def print_summary(self):
        """Print configuration summary"""
        print("=" * 70)
        print("CONFIGURATION SUMMARY")
        print("=" * 70)

        print(f"\nConfig file: {self.config_path}")

        print("\n[PATHS]")
        for key, val in self.paths.items():
            print(f"  {key}: {val}")

        print("\n[TISSUES]")
        for cfg in self.tissues:
            print(f"  - {cfg['tissue']}")
            print(f"    slow: {cfg['slow_model']}")
            print(f"    fast: {cfg['fast_model']}")

        print("\n[ANALYSIS]")
        for key, val in self.analysis.items():
            print(f"  {key}: {val}")

        print("\n[PATHWAY MERGING]")
        if self.pathway_merging:
            for new_name, old_names in self.pathway_merging.items():
                print(f"  {new_name}:")
                for old in old_names:
                    print(f"    - {old}")
        else:
            print("  (disabled)")

        print("\n[PATHWAY FILTERING]")
        mode = self.pathway_filter.get("mode", "none")
        print(f"  mode: {mode}")
        if mode == "whitelist":
            print(f"  whitelist ({len(self.pathway_filter.get('whitelist', []))} pathways):")
            for pw in self.pathway_filter.get("whitelist", [])[:5]:
                print(f"    - {pw}")
            if len(self.pathway_filter.get("whitelist", [])) > 5:
                print(f"    ... and {len(self.pathway_filter.get('whitelist', [])) - 5} more")
        elif mode == "blacklist":
            print(f"  blacklist: {self.pathway_filter.get('blacklist', [])}")

        print("\n[SCATTER]")
        for key, val in self.scatter.items():
            print(f"  {key}: {val}")

        print("\n" + "=" * 70)


def load_config(config_path: str = "input_parameters.yaml") -> AnalysisConfig:
    """
    Load configuration from YAML file

    Parameters:
    -----------
    config_path : str
        Path to configuration file

    Returns:
    --------
    AnalysisConfig
        Configuration object
    """
    return AnalysisConfig(config_path)


if __name__ == "__main__":
    # Test configuration loading
    logging.basicConfig(level=logging.INFO)

    try:
        config = load_config("input_parameters.yaml")
        config.print_summary()

        print("\n[VALIDATION]")
        config.validate()
        print("✓ Configuration is valid")

        print("\n[TESTING PATHWAY MERGING]")
        test_pathways = [
            "Transport, extracellular",
            "Transport, mitochondrial",
            "Citric acid cycle",
            "Valine, leucine, and isoleucine metabolism"
        ]
        print(f"Original: {test_pathways}")
        merged = config.apply_pathway_merging(test_pathways)
        print(f"Merged:   {merged}")

        print("\n[TESTING PATHWAY FILTERING]")
        filtered = config.filter_pathways(merged)
        print(f"Filtered: {filtered}")

    except Exception as e:
        print(f"✗ Error: {e}")
        import traceback
        traceback.print_exc()
