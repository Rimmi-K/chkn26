# chkn26

[![Python](https://img.shields.io/badge/Python-3.12-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-%E2%89%A54.3-blue.svg)](https://www.r-project.org/)
[![COBRApy](https://img.shields.io/badge/COBRApy-0.29.1-green.svg)](https://opencobra.github.io/cobrapy/)

Multi-omics analysis integrating transcriptomics, metabolomics, and constraint-based modeling to investigate tissue-specific metabolic adaptations in fast- vs slow-growing chickens.

---

## Repository structure

```
data/
├── transcriptomics/    # TMM-normalised matrices, DESeq2 DEGs, GSEA results
├── metabolomics/       # NMR concentrations (nmol/g wet tissue)
├── phenotypes/         # Body & organ weights, serum biochemistry
└── models/             # iES1470 GEM (SBML) + RIPTiDe context-specific models

scrs/
├── transcriptomics/    # CAGE-seq QC, DESeq2, GSEA (R scripts)
├── metabolomics/       # Hodges–Lehmann analysis, bootstrap CIs
├── dpfa/               # Differential Pathway Flux Analysis package
├── tfca/               # Transcript-Flux Concordance Analysis
├── riptide_integration/# RIPTiDe integration
└── models/             # scripts for refine model

results/                # Generated figures, tables, supplementary files
```

---

## Methods

- TMM normalization & DESeq2 differential expression
- GSEA pathway enrichment (GO:BP)
- Metabolic model curation (iES1470: 1,470 genes, 2,698 reactions)
- RIPTiDe integration (fraction 0.10–0.95)
- Differential pathway flux analysis (DPFA: r-DPFA + mc-DPFA)

---

## The DPFA Package

Differential Pathway Flux Analysis toolkit for comparing context-specific flux distributions.

**Modules:**
- `dpfa.analysis` — DRF classification, r-DPFA, mc-DPFA
- `dpfa.visualization` — pathway bar charts, metabolite heatmaps
- `dpfa.scatter_plot` — transcript–flux concordance analysis: compares transcript-derived reaction effects (aggregated via GPR rules) with predicted flux ratios to identify concordant, flux-dominant, and transcript-dominant regulation
- `dpfa.utils` — GPR rule parsing, metabolite turnover calculations

**Usage:** Configure `input_parameters.yaml`, run `python -m scrs.dpfa`

---

## Reproducibility

### Computational environment

Analysis performed using three conda environments:
- **`base`**: Python 3.12 + COBRApy 0.29.1 (flux analysis, DPFA, TFCA)
- **`r_gsea_bio`**: R 4.3+ with Bioconductor (DESeq2, clusterProfiler, GSEA)
- **`riptide`**: RIPTiDe context-specific model integration

Environment specifications provided in `environment_*.yml` files. To recreate environments:

```bash
conda env create -f environment_base.yml
conda env create -f environment_r_gsea.yml
conda env create -f environment_riptide.yml
```

### Analysis workflow

**1. Transcriptomics** (R, `r_gsea_bio` env)
```bash
conda activate r_gsea_bio
cd scrs/transcriptomics/
Rscript TMM.R              # TMM normalization
Rscript DEG_analysis.R     # DESeq2 differential expression
Rscript GSEA.R             # Gene set enrichment analysis
```

**2. Metabolomics** (Python, `base` env)
```bash
conda activate base
python scrs/metabolomics/metabolomics_stats.py
```

**3. Context-specific models** (Python, `riptide` env)
```bash
conda activate riptide
# Generate tissue-specific models using RIPTiDe
# (scripts in scrs/riptide_integration/)
```

**4. Flux analysis** (Python, `base` env)
```bash
conda activate base
# Configure tissue-specific parameters in input_parameters.yaml
python -m scrs.dpfa        # Differential pathway flux analysis
python -m scrs.tfca        # Transcript-flux concordance analysis
```

---

## Citation


```bibtex
@article{kaplan2026chkn,
  title   = {Genome-scale and Omics-Driven Modelling Reveals Metabolic
             Differences Between Fast- and Slow-Growing Chicken Groups},
  author  = {Kaplan, Vladimir and Evshin, Ivan and Osik, Nataliya and
             Yanshole, Lyudmila and Tsentalovich, Yuri and
             Shagimardanova, Elena and Gusev, Oleg and Kolpakov, Fedor and
             Kulyashov, Mikhail and Akberdin, Ilya},
  journal = {submitted},
  year    = {2026}
}
```

---

## Contact

kaplanrimi@gmail.com
---
