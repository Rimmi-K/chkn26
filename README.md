# chkn26

# Metabolic and transcriptomic analysis of fast- and slow-growing chickens

[![Python](https://img.shields.io/badge/Python-3.12-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-%E2%89%A54.3-blue.svg)](https://www.r-project.org/)
[![COBRApy](https://img.shields.io/badge/COBRApy-0.29.1-green.svg)](https://opencobra.github.io/cobrapy/)
[![RIPTiDe](https://img.shields.io/badge/RIPTiDe-3.4.81-green.svg)](https://github.com/mjenior/riptide)

This repository contains a multi-omics analysis of fast- and slow-growing chickens, integrating **transcriptomics**, **metabolomics**, **phenotypes** to investigate tissue-specific metabolic differences.
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
- [License](#license)
- [Citation](#citation)
- [Contact](#contact)
- [Acknowledgments](#acknowledgments)

---

## Introduction 

Growth rate in poultry is a complex trait shaped by coordinated transcriptional and metabolic regulation across tissues. While transcriptomic differences between fast- and slow-growing chickens can be detected, understanding their **mechanistic metabolic consequences** requires a systems-level approach.

Genome-scale metabolic models (GEMs), combined with expression integration methods (e.g., RIPTiDe), provide a framework to translate expression patterns into **flux-level metabolic changes**, which can then be compared with metabolomic measurements and phenotypic traits.

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

All data are stored under `data/`. Results of analysis are stored under `results/`.

---

## Repository structure

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


#### 3. Set up R environment

Install R (≥ 4.3) and required packages:

```r
install.packages(c("edgeR", "dplyr", "readr", "ggplot2", "tidyr", "effsize"))
```

---

## Running the analysis

### Step 0: Statistical analysis

Run metabolomic and phenotypic statistical analyses:

```bash
# Metabolomics statistics
python scrs/metabolomics/metabolomics_stats.py

# Phenotype statistics
Rscript scrs/phenotypes/PhenStatAn.Rscript

```

### Step 1: Transcriptomic preprocessing

Run TMM normalization and differential expression analysis using R scripts in `scrs/transcriptomics/`:

```bash
Rscript scrs/transcriptomics/TMM.R
Rscript scrs/transcriptomics/DEG_analysis.R
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

First of all, check parameters in config file - input_parameters.yaml

Then, you can just exec dpfa:
```bash
python -m dpfa
```


## The dpfa package

The `dpfa` (Downstream Pathway and Flux Analysis) package provides tools for:

- **Flux analysis** (`dpfa/analysis.py`) — differential flux analysis between conditions
- **Visualization** (`dpfa/visualization.py`) — plotting flux distributions, pathway activities
- **Scatter plots** (`dpfa/scatter_plot.py`) — correlation and comparison plots
- **Configuration** (`dpfa/config_loader.py`) — YAML parameter loading
- **Utilities** (`dpfa/utils/`) — pathway databases, GPR parsing, flux calculations, metabolite turnover

## License


---

## Citation

If you use this work in your research, please cite:

---

## Contact

kaplanrimi@gmail.com

---

## Acknowledgments

This work uses:
- [COBRApy](https://opencobra.github.io/cobrapy/) for constraint-based metabolic modeling
- [RIPTiDe](https://github.com/mjenior/riptide) for transcriptomic integration
- [edgeR](https://bioconductor.org/packages/edgeR/) for differential expression analysis
