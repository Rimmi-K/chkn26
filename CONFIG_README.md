# Configuration System for Metabolic Flux Analysis

## Overview

The analysis pipeline now uses a YAML-based configuration system that allows you to:
- **Centralize all parameters** in one file (`input_parameters.yaml`)
- **Merge metabolic pathways** under common names
- **Filter pathways** for analysis (whitelist/blacklist)
- **Customize visualization** and analysis thresholds
- **Manage multiple tissues** with different models

---

## Quick Start

### 1. Edit Configuration

All parameters are in `input_parameters.yaml`. Edit this file to configure your analysis:

```bash
nano input_parameters.yaml
```

### 2. Run Analysis

```bash
python -m scrs.dpfa.analysis
```

The script will automatically load `input_parameters.yaml` from the project root.

---

## Configuration File Structure

### Paths
```yaml
paths:
  base_model: "data/models/curated_model55.xml"
  pathway_database: "data/models/subsystem_matrix.csv"
  base_output_dir: "results_riptide"
```

### Tissue Configurations
```yaml
tissues:
  - tissue: "breast"
    slow_model: "data/models/context_specific/mcs_slow_breast_fraction_0.65.json"
    fast_model: "data/models/context_specific/mcs_fast_breast_fraction_0.90.json"
    deg_csv: "data/processed/transcriptomics/degs/degs_results_breast.csv"
    output_dir: "results_riptide/breast01"
```

### Analysis Parameters
```yaml
analysis:
  fdr_threshold: 0.99              # FDR for DEG significance
  flux_diff_threshold: 0.0         # Minimum flux difference
  gpr_missing_strategy: "impute"   # "impute" or "drop"
  gpr_impute_value: 0.0            # Value for missing genes
```

### Scatter Plot Parameters
```yaml
scatter:
  rxn_lfc_threshold: 0.5           # Reaction log2FC threshold
  flux_threshold: 0.5              # Flux log2 ratio threshold
  size_mode: "discrete"            # "discrete" or "continuous"
  size_thresholds: [0.001, 0.01, 0.05, 0.25]
  size_values: [175, 150, 100, 50]
```

---

## Pathway Merging

### Purpose
Merge multiple related pathways under a single name to reduce complexity and improve visualization.

### Configuration
```yaml
pathway_merging:
  "Transport, all":
    - "Transport, extracellular"
    - "Transport, mitochondrial"
    - "Transport, peroxisomal"
    - "Transport, endoplasmic reticular"

  "Nucleotide metabolism":
    - "Nucleotide interconversion"
    - "Purine metabolism"
    - "Pyrimidine metabolism"
```

### How It Works
1. All occurrences of old pathway names are replaced with the new merged name
2. Applies to **all analysis outputs**: DRF tables, scatter plots, heatmaps, legends
3. Flux values are aggregated (summed) for merged pathways

### Disabling
Set to empty dictionary:
```yaml
pathway_merging: {}
```

---

## Pathway Filtering

### Purpose
Control which pathways are analyzed in DRF bar plots, mcPFA heatmaps, and pathway-level statistics.

**Note**: Scatter plots always show **all pathways** (no filtering applied).

### Modes

#### Whitelist (analyze ONLY listed pathways)
```yaml
pathway_filter:
  mode: "whitelist"
  whitelist:
    - "Oxidative phosphorylation"
    - "Citric acid cycle"
    - "Glycolysis/gluconeogenesis"
    - "Transport, all"  # Use merged names here!
```

#### Blacklist (analyze ALL EXCEPT listed pathways)
```yaml
pathway_filter:
  mode: "blacklist"
  blacklist:
    - "Unknown"
    - "Other"
```

#### None (no filtering)
```yaml
pathway_filter:
  mode: "none"
```

### Important Notes
- Filtering is applied **AFTER** pathway merging
- Use **merged names** in whitelist/blacklist (e.g., "Transport, all" not "Transport, extracellular")
- Scatter plots ignore filtering (all pathways shown for comprehensive view)

---

## Using Configuration in Python

### Loading Configuration
```python
from scrs.dpfa.config_loader import load_config

config = load_config("input_parameters.yaml")
config.validate()  # Check that all paths exist
```

### Accessing Parameters
```python
# Paths
base_model = config.get_base_model_path()
pathway_db = config.get_pathway_database_path()

# Analysis params
fdr_thr = config.get_fdr_threshold()
gpr_strategy = config.analysis["gpr_missing_strategy"]

# Scatter params
scatter_params = config.get_scatter_params()
rxn_lfc_thr = scatter_params["rxn_lfc_threshold"]

# Tissue configs
for tissue_cfg in config.get_tissue_configs():
    print(f"Processing {tissue_cfg['tissue']}")
```

### Applying Pathway Merging
```python
from scrs.dpfa.utils.pathway_utils import (
    apply_pathway_merging_to_df,
    aggregate_pathways_in_df
)

# Merge pathway names in DataFrame
merged_df = apply_pathway_merging_to_df(
    df,
    config.pathway_merging,
    pathway_column="Pathways"
)

# Merge and aggregate flux values
aggregated_df = aggregate_pathways_in_df(
    df,
    config.pathway_merging,
    pathway_column="Pathways",
    value_columns=["flux", "log2fc"],
    aggregation="sum"  # or "mean", "median", etc.
)
```

### Applying Pathway Filtering
```python
from scrs.dpfa.utils.pathway_utils import filter_pathways_in_df

# Filter DataFrame based on config
filtered_df = filter_pathways_in_df(
    merged_df,  # Apply AFTER merging!
    mode=config.pathway_filter["mode"],
    whitelist=config.pathway_filter.get("whitelist"),
    blacklist=config.pathway_filter.get("blacklist"),
    pathway_column="Pathways"
)
```

### Using with scatter_plot.py
```python
from scrs.dpfa.scatter_plot import make_scatter_deg_vs_flux

make_scatter_deg_vs_flux(
    merged_df,
    rxn_lfc_map,
    rxn_sig_map,
    tissue="breast",
    output_dir="results/breast",

    # Pass pathway merging from config
    pathway_merging=config.pathway_merging,

    # Use scatter parameters from config
    **config.get_scatter_params()  # Unpacks all scatter params
)
```

---

## Validation

Validate your configuration before running analysis:

```bash
python -m scrs.dpfa.config_loader
```

This will:
- Check that all model files exist
- Check that all DEG files exist
- Verify pathway filter mode is valid
- Print configuration summary

---

## Example Workflow

See `example_config_usage.py` for a complete example:

```bash
python example_config_usage.py
```

This demonstrates:
- Loading configuration
- Applying pathway merging
- Filtering pathways
- Checking pathway inclusion
- Using with scatter plots

---

## Migrating Existing Code

### Old Way (Hard-coded Parameters)
```python
tissue_configs = [
    {
        "tissue": "breast",
        "slow_model": "data/models/.../mcs_slow_breast_fraction_0.65.json",
        "fast_model": "data/models/.../mcs_fast_breast_fraction_0.90.json",
        # ... more hard-coded paths
    }
]

run_analysis(
    model=model,
    fdr_thr=0.99,
    scatter_rxn_lfc_thr=0.5,
    # ... 20 more parameters
)
```

### New Way (Config-based)
```python
config = load_config("input_parameters.yaml")

for tissue_cfg in config.get_tissue_configs():
    run_analysis(
        model=model,
        **tissue_cfg,  # Unpack tissue config
        fdr_thr=config.get_fdr_threshold(),
        pathway_merging=config.pathway_merging,
        pathway_filter=config.pathway_filter,
        **config.analysis  # Unpack all analysis params
    )
```

---

## Troubleshooting

### "PyYAML not available"
Install PyYAML:
```bash
pip install pyyaml
```

### "Config file not found"
Ensure `input_parameters.yaml` is in the project root (same directory as where you run `python -m scrs.dpfa.analysis`).

### "Pathway not found after filtering"
- Check if pathway name matches exactly (case-sensitive)
- If using pathway merging, use the **merged name** in whitelist/blacklist
- Set `mode: "none"` to disable filtering temporarily

### "Duplicate colors in scatter plot"
This can happen if many pathways merge into one. The golden-angle HSV palette handles 60+ unique colors, so this should be rare.

---

## Advanced: Custom Aggregation

When pathways are merged, you can control how numeric values are combined:

```python
# Sum fluxes (default)
aggregated_df = aggregate_pathways_in_df(df, merging_dict, aggregation="sum")

# Average log2FC values
aggregated_df = aggregate_pathways_in_df(df, merging_dict, aggregation="mean")

# Maximum flux
aggregated_df = aggregate_pathways_in_df(df, merging_dict, aggregation="max")
```

---

## File Locations

```
chkn26/
├── input_parameters.yaml           # Main config file (EDIT THIS)
├── scrs/dpfa/
│   ├── config_loader.py            # Config loading module
│   ├── utils/
│   │   └── pathway_utils.py        # Pathway manipulation utilities
│   └── scatter_plot.py             # Supports pathway_merging parameter
├── example_config_usage.py         # Usage examples
└── CONFIG_README.md                # This file
```

---

## Support

For issues or questions:
1. Check `example_config_usage.py` for working examples
2. Validate config: `python -m scrs.dpfa.config_loader`
3. Review error messages - they usually indicate which parameter is wrong

---

**Happy analyzing!** 🧬🔬
