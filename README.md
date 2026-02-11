# chkn26

[![Python](https://img.shields.io/badge/Python-3.12-blue.svg)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-%E2%89%A54.3-blue.svg)](https://www.r-project.org/)
[![COBRApy](https://img.shields.io/badge/COBRApy-0.29.1-green.svg)](https://opencobra.github.io/cobrapy/)

Multi-omics analysis integrating transcriptomics, metabolomics, and constraint-based modeling to investigate tissue-specific metabolic adaptations in fast- vs slow-growing chickens.

---

**Requirements:** Python ≥3.10, R ≥4.3, 8GB RAM  
**Output:** Results → `results/` (organized by tissue/analysis type)

---

## Data

**Design:** 9-week *Gallus gallus*, fast-growing (n=6) vs slow-growing (n=6)  
**Tissues:** Liver, leg muscle, breast muscle  
**Modalities:** CAGE-seq, NMR metabolomics, phenotypes, context-specific GEMs
```
data/
├── transcriptomics/  # TMM matrices, DEGs, GSEA
├── metabolomics/     # NMR concentrations
├── phenotypes/       # Growth traits
└── models/           # iES1470 GEM + context-specific models
```

---

## Methods

- TMM normalization & DESeq2 differential expression
- GSEA pathway enrichment (GO:BP)
- Metabolic model curation (iES1470: 1,470 genes, 2,698 reactions)
- RIPTiDe integration (fraction 0.10–0.95)
- Differential pathway flux analysis (DPFA: r-DPFA + mc-DPFA)

**Full methods:** See [publication](link) or `scrs/` scripts

---

## Repository Structure
```
data/          # Omics data & models
scrs/          # Analysis scripts
  ├── transcriptomics/
  ├── metabolomics/
  ├── riptide_integration/
  └── models/
dpfa/          # Flux analysis package
results/       # Generated outputs
```

---

## The dpfa Package

Differential pathway flux analysis toolkit:
- `dpfa.analysis` — flux comparisons (FG vs SG)
- `dpfa.visualization` — pathway/metabolite heatmaps
- `dpfa.utils` — GPR parsing, metabolite turnover

**Usage:** Configure `input_parameters.yaml`, run `python -m dpfa`

---

## Citation

```bibtex
@article{chkn2026,
  title={Genome-scale and Omics-Driven Modelling Reveals Metabolic Differences Between Fast- and Slow-Growing Chicken Groups},
  author={...},
  journal={...},
  year={2026}
}
```

---

## Contact

📧 kaplanrimi@gmail.com  
---
