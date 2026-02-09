# chkn26

# Metabolic and transcriptomic analysis of fast- and slow-growing chickens

[![Python](https://img.shields.io/badge/Python-3.12-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-%E2%89%A54.3-blue.svg)](https://www.r-project.org/)
[![COBRApy](https://img.shields.io/badge/COBRApy-0.29.1-green.svg)](https://opencobra.github.io/cobrapy/)
[![RIPTiDe](https://img.shields.io/badge/RIPTiDe-3.4.81-green.svg)](https://github.com/mjenior/riptide)

This repository contains a multi-omics analysis of fast- and slow-growing chickens, integrating **transcriptomics**, **metabolomics**, **phenotypes**, and **context-specific genome-scale metabolic models (csGEMs)** to investigate tissue-specific metabolic differences.
The main focus is on **liver**, **leg muscle**, and **breast muscle** tissues.

## Table of Contents

- [Background](#background)
- [Research questions](#research-questions)
- [Data description](#data-description)
- [Repository structure](#repository-structure)
- [Methods overview](#methods-overview)
- [How to reproduce the analysis](#how-to-reproduce-the-analysis)
  - [Requirements](#requirements)
  - [Installation](#installation)
  - [Running the analysis](#running-the-analysis)
- [The dpfa package](#the-dpfa-package)
- [Configuration](#configuration)
- [Expected outputs](#expected-outputs)
- [Troubleshooting](#troubleshooting)
- [License](#license)
- [Citation](#citation)
- [Contact](#contact)
- [Acknowledgments](#acknowledgments)

---

## Background

Growth rate in poultry is a complex trait shaped by coordinated transcriptional and metabolic regulation across tissues. While transcriptomic differences between fast- and slow-growing chickens can be detected, understanding their **mechanistic metabolic consequences** requires a systems-level approach.

Genome-scale metabolic models (GEMs), combined with expression integration methods (e.g., RIPTiDe), provide a framework to translate expression patterns into **flux-level metabolic changes**, which can then be compared with metabolomic measurements and phenotypic traits.

---

## Research questions

This project addresses the following questions:

- Which metabolic pathways differ between fast- and slow-growing chickens?
- Are differences tissue-specific (liver vs muscle)?
- How do transcriptomic changes translate into metabolic flux alterations?
- Are metabolomic and phenotypic measurements consistent with model predictions?

---

## Data description

### Organism and experimental design
- **Species:** *Gallus gallus*
- **Age:** 9 weeks  
- **Groups:** fast-growing (n = 6) vs slow-growing (n = 6)
- **Tissues (main focus):** liver, leg muscle, breast muscle

### Data types
- **Transcriptomics:** CAGE-seq–derived gene expression; includes TMM-normalized matrices, DEGs, GSEA results
- **Metabolomics:** metabolite tables and concentrations per tissue; summary statistics and visualizations
- **Phenotypes:** biomass-related / growth-associated traits and statistical summaries
- **Models:** curated chicken GEM and a collection of context-specific models produced across multiple expression thresholds

All data are stored under `data/`. Generated outputs are stored under `results/`.

---

## Repository structure

High-level map of the repository:

- `data/` — raw and processed datasets used in analysis  
  - `data/transcriptomics/` — expression matrices, counts, DEGs, GSEA outputs  
  - `data/metabolomics/` — metabolite concentrations and metabolomics statistics inputs  
  - `data/phenotypes/` — phenotype tables and sample metadata  
  - `data/models/` — curated GEM + context-specific models + subsystem matrix

- `scrs/` — scripts for building/processing data and models  
  - `scrs/transcriptomics/` — TMM normalization, DEG analysis, GSEA, metadata creation  
  - `scrs/riptide_integration/` — RIPTiDe integration pipeline  
  - `scrs/models/` — model preparation, annotation enrichment, subsystem processing  
  - `scrs/metabolomics/` — metabolomics statistics scripts  
  - `scrs/phenotypes/` — phenotype statistics scripts (R)

- `dpfa/` — python package for downstream flux/pathway analyses and plotting utilities

- `results/` — generated figures/tables organized by data modality  
  - `results/fluxomics/` — flux-based outputs per tissue  
  - `results/metabolomics/` — metabolomics plots (e.g., forest plots)  
  - `results/phenotypes/` — phenotype analysis results and plots  
  - `results/transcriptomics/` — transcriptomic summaries (e.g., GSEA)

Other important files:
- `input_parameters.yaml` — main configuration for the pipeline
- `CONFIG_README.md` — explanation of configuration parameters
- `envs/requirements.txt` — python dependencies
- `riptide_log.txt`, `riptide_summary.csv` — RIPTiDe run logs/summary

---

## Methods overview

The analysis workflow includes:

1. **Transcriptomic preprocessing and normalization** (TMM)
2. **Differential expression analysis** (fast vs slow) per tissue
3. **Pathway enrichment analysis** (GSEA)
4. **Curation and annotation** of a chicken genome-scale metabolic model
5. **Integration of transcriptomic data into GEMs** using RIPTiDe
6. **Construction of context-specific models** across multiple expression thresholds (`fraction_0.10` … `fraction_0.95`)
7. **Flux analysis and pathway-level interpretation**
8. **Cross-validation** using metabolomics and phenotypes

---

## How to reproduce the analysis

### Requirements

- **Python:** ≥ 3.10 (recommended: 3.12)
- **R:** ≥ 4.3
- **RAM:** ≥ 8 GB recommended (genome-scale models can be memory-intensive)
- **Disk space:** ≥ 5 GB for data, models, and results
- Python dependencies: see `envs/requirements.txt`

> **Note:** R packages are installed separately via your R environment; key packages include `edgeR` for differential expression and various plotting/statistics libraries.

### Installation

#### 1. Clone the repository

```bash
git clone https://github.com/YOUR_USERNAME/chkn26.git
cd chkn26
```

#### 2. Set up Python environment

It is recommended to use a conda/mamba environment:

```bash
conda create -n chkn26 python=3.12
conda activate chkn26
```

Install Python dependencies:

```bash
pip install -r envs/requirements.txt
```

Key Python packages include:
- `cobra` (0.29.1) — genome-scale metabolic modeling
- `riptide` (3.4.81) — transcriptomic integration
- `pandas`, `numpy`, `scipy` — data manipulation and analysis
- `matplotlib`, `seaborn` — visualization
- `pydeseq2` — differential expression analysis

#### 3. Set up R environment

Install R (≥ 4.3) and required packages:

```r
install.packages(c("edgeR", "dplyr", "readr", "ggplot2", "tidyr", "effsize"))
```

Key R packages:
- `edgeR` — TMM normalization and differential expression
- `dplyr`, `tidyr` — data wrangling
- `ggplot2` — plotting
- `effsize` — effect size calculations

---

## Running the analysis

### Step 1: Transcriptomic preprocessing

Run TMM normalization and differential expression analysis using R scripts in `scrs/transcriptomics/`:

```bash
# TMM normalization
Rscript scrs/transcriptomics/TMM.R

# Differential expression analysis
Rscript scrs/transcriptomics/DEG_analysis.R

# Gene Set Enrichment Analysis
Rscript scrs/transcriptomics/GSEA.R
```

### Step 2: Model preparation

Prepare and annotate the genome-scale metabolic model:

```bash
python scrs/models/prepare_model.py
python scrs/models/enrich_annotations.py
python scrs/models/update_subsystem_matrix.py
```

### Step 3: RIPTiDe integration

Integrate transcriptomic data into metabolic models to generate context-specific models:

```bash
python scrs/riptide_integration/integrate_riptide.py
```

This will:
- Load the curated chicken GEM
- Read transcriptomic data for fast and slow-growing groups
- Generate context-specific models across multiple expression thresholds (fraction 0.10 to 0.95)
- Save models to `data/models/` with naming pattern: `mcs_{growth}_{tissue}_fraction_{threshold}.json`

Output logs are saved to `riptide_log.txt` and `riptide_summary.csv`.

### Step 4: Flux analysis

Use the `dpfa` package to analyze flux distributions and compare metabolic pathways:

```python
from dpfa import analysis, visualization

# Load configuration
from dpfa.config_loader import load_config
config = load_config('input_parameters.yaml')

# Perform flux analysis
# See dpfa/analysis.py for available functions
```

### Step 5: Statistical analysis

Run metabolomic and phenotypic statistical analyses:

```bash
# Metabolomics statistics
python scrs/metabolomics/metabolomics_stats.py

# Phenotype statistics
Rscript scrs/phenotypes/PhenStatAn.R
```

---

## The dpfa package

The `dpfa` (Downstream Pathway and Flux Analysis) package provides tools for:

- **Flux analysis** (`dpfa/analysis.py`) — differential flux analysis between conditions
- **Visualization** (`dpfa/visualization.py`) — plotting flux distributions, pathway activities
- **Scatter plots** (`dpfa/scatter_plot.py`) — correlation and comparison plots
- **Configuration** (`dpfa/config_loader.py`) — YAML parameter loading
- **Utilities** (`dpfa/utils/`) — pathway databases, GPR parsing, flux calculations, metabolite turnover

Example usage:

```python
from dpfa.analysis import compare_fluxes
from dpfa.visualization import plot_pathway_activity

# Load models
# ... load context-specific models ...

# Compare flux between fast and slow groups
flux_diff = compare_fluxes(model_fast, model_slow)

# Visualize pathway activity
plot_pathway_activity(flux_diff, save_path='results/fluxomics/')
```

---

## Configuration

Pipeline settings are centralized in `input_parameters.yaml`, which controls:

- **Paths:** locations of models, databases, and output directories
- **Tissue configurations:** model files, DEG data, and fraction thresholds for each tissue
- **Pathway settings:** pathway merging, filtering (whitelist/blacklist), and naming
- **Visualization parameters:** plot settings, statistical thresholds, color schemes
- **Analysis options:** objectives, optimization parameters, flux calculations

**Configuration file structure:**

```yaml
paths:
  base_model: "data/models/curated_model55.xml"
  pathway_database: "data/models/subsystem_matrix.csv"
  base_output_dir: "results_riptide"

tissues:
  - tissue: "breast"
    slow_model: "data/models/context_specific/mcs_slow_breast_fraction_0.65.json"
    fast_model: "data/models/context_specific/mcs_fast_breast_fraction_0.90.json"
    # ... additional tissue-specific settings

pathway_settings:
  merge_pathways:
    "Unified Pathway Name":
      - "Original pathway 1"
      - "Original pathway 2"
  # ... filtering and exclusion rules
```

See `CONFIG_README.md` for detailed documentation of all configuration options.

> **Important:** If you change parameters, keep a copy of the configuration file with your results for reproducibility.

---

## Expected outputs

After running the full pipeline, you should have:

- `data/models/` — context-specific metabolic models (JSON format)
- `results/transcriptomics/` — DEG lists, GSEA results, normalized expression matrices
- `results/fluxomics/` — flux analysis results and pathway comparisons
- `results/metabolomics/` — metabolomics plots (forest plots, volcano plots)
- `results/phenotypes/` — phenotype statistics and visualizations
- `riptide_log.txt` — detailed RIPTiDe integration log
- `riptide_summary.csv` — summary of RIPTiDe runs

**Computation time:** The full pipeline (especially RIPTiDe integration across all thresholds and tissues) may take several hours depending on your hardware.

---

## Troubleshooting

### Common issues

**1. Import errors with COBRApy or RIPTiDe**
- Ensure you're using Python 3.12 and have installed all dependencies from `envs/requirements.txt`
- Try creating a fresh conda environment

**2. R script errors (e.g., `edgeR` not found)**
- Install Bioconductor packages:
  ```r
  if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("edgeR")
  ```

**3. Memory errors during RIPTiDe integration**
- Reduce the number of fraction thresholds being tested
- Process tissues sequentially rather than in parallel
- Increase system swap space

**4. File path errors in R scripts**
- R scripts may contain hardcoded paths (e.g., `C:/Users/...`). Update these to match your directory structure
- Check paths in `TMM.R`, `DEG_analysis.R`, `GSEA.R`, and `PhenStatAn.R`

**5. Model optimization failures**
- Check that model files are not corrupted
- Verify that the objective function is set correctly in `input_parameters.yaml`
- Ensure metabolite/reaction IDs match between model and transcriptome data

For additional help, please open an issue on GitHub with:
- Error message and full traceback
- Python/R version
- Operating system
- Steps to reproduce

---

## License

This project is currently unlicensed. Please contact the authors for usage permissions.

---

## Citation

If you use this work in your research, please cite:

```
[Add citation information here when published]
```

---

## Contact

For questions or issues, please open an issue on GitHub or contact the authors.

---

## Acknowledgments

This work uses:
- [COBRApy](https://opencobra.github.io/cobrapy/) for constraint-based metabolic modeling
- [RIPTiDe](https://github.com/mjenior/riptide) for transcriptomic integration
- [edgeR](https://bioconductor.org/packages/edgeR/) for differential expression analysis
