# Environment Setup

## Python Environment

### Using conda (recommended):
```bash
conda create -n chicken-gsm python=3.12
conda activate chicken-gsm
conda install -c conda-forge cobra=0.29.1
pip install -r requirements.txt
```

### Using pip:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt
```

## R Environment

Required R packages are installed via R scripts in `scrs/data_preparation/transcriptomics/`:

```R
# Install Bioconductor packages
install.packages('BiocManager')
BiocManager::install(c('DESeq2', 'edgeR', 'clusterProfiler', 'org.Gg.eg.db'))

# Install CRAN packages
install.packages(c('ggplot2', 'dplyr', 'tidyr'))
```

## Key Dependencies

- **Python**: 3.12+
- **R**: 4.3.2+
- **COBRApy**: 0.29.1
- **RIPTiDe**: 3.4.81
- **Solver**: GLPK (installed via conda/pip)
