# SpatialBridge

> HPC-native, multi-platform spatial transcriptomics pipeline for tumor microenvironment analysis

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A522.10.0-brightgreen.svg)](https://www.nextflow.io/)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

SpatialBridge is a reproducible, HPC-deployable computational pipeline for end-to-end spatial transcriptomics analysis in oncology. It processes 10x Genomics Visium data through quality control, normalization, clustering, cell type deconvolution, and tumor microenvironment (TME) characterization — all within a single Nextflow-orchestrated workflow.

The pipeline addresses a key infrastructure gap: no existing open-source tool combines multi-platform spatial transcriptomics analysis with HPC-native execution, automated GPU-based deconvolution via cloud API, and clinically interpretable TME output.

---

## Pipeline Architecture

```
Raw Visium Data (.h5 + spatial/)
          ↓
    QC Module (Scanpy)
    - Adaptive MAD-based thresholding
    - Mitochondrial % filtering
    - Tissue spot filtering
          ↓
  Normalization Module
    - Library size normalization (10k counts)
    - Log1p transformation
    - Highly Variable Gene selection (3000 HVGs)
          ↓
  Clustering Module (Scanpy + Squidpy)
    - PCA (50 components)
    - Neighborhood graph (k=15)
    - Leiden clustering (resolution=0.5)
    - UMAP visualization
          ↓
  Cell Type Deconvolution (Cell2location)
    - GPU execution via Kaggle API
    - Reference: Wu et al. 2021 Breast Cancer scRNA-seq Atlas
    - Bayesian spatial mapping (200 + 300 epochs)
    - Auto-fetches results back to HPC
          ↓
  TME Analysis (Squidpy)
    - Dominant cell type per spot
    - Immune infiltration scoring
    - Neighborhood enrichment analysis
    - Moran's I spatial autocorrelation
          ↓
  Results (CSV + PNG outputs)
```

---

## Key Features

- **Nextflow DSL2 orchestration** — modular, resumable, reproducible
- **CDAC HPC native** — tested on CDAC ICE Cloud (8 vCPU, 32 GB RAM)
- **Kaggle GPU integration** — automated Cell2location training via Kaggle API; results auto-fetched to HPC
- **Adaptive QC** — MAD-based thresholding (no manual cutoffs)
- **Spatially-aware TME analysis** — Moran's I autocorrelation for spatial pattern detection
- **Nextflow `-resume`** — checkpoint-based resumption; no recomputation on restart

---

## Results on Human Breast Cancer Dataset

**Dataset:** 10x Genomics Visium V1 Human Breast Cancer Block A Section 1

| Metric | Value |
|---|---|
| Total spots (raw) | 3,798 |
| Spots after QC | 2,851 |
| Spots removed | 947 (24.9%) |
| Genes analyzed | 36,601 |
| Clusters identified | 10 |

**Cell Type Deconvolution Results (Cell2location):**

| Cell Type | Mean Abundance |
|---|---|
| Cancer Epithelial | 0.973 |
| Normal Epithelial | 0.894 |
| CAFs | 0.820 |
| Endothelial | 0.810 |
| Myeloid | 0.801 |
| T-cells | 0.797 |
| B-cells | 0.760 |
| Plasmablasts | 0.565 |

**Spatial Autocorrelation (Moran's I):**

| Cell Type | Moran's I | p-value |
|---|---|---|
| Normal Epithelial | 0.795 | < 0.001 |
| Cancer Epithelial | 0.780 | < 0.001 |
| PVL | 0.715 | < 0.001 |

High Moran's I values confirm non-random spatial organization — tumor and immune cells form distinct spatial domains within the tissue.

---

## Requirements

- Nextflow >= 22.10.0
- Python 3.10+
- CDAC ICE Cloud or any Linux HPC (8 cores, 32 GB RAM minimum)
- Kaggle account with GPU access (free T4 x2)
- Internet access on HPC

**Python packages:**
```
scanpy>=1.9
squidpy>=1.3
cell2location>=0.1.5
leidenalg
igraph
pandas
numpy
scipy
matplotlib
```

---

## Installation

```bash
git clone https://github.com/TusharParihar128/SpatialBridge.git
cd SpatialBridge
pip install scanpy squidpy cell2location leidenalg igraph --break-system-packages
```

**Kaggle CLI setup:**
```bash
pip install kaggle --break-system-packages
mkdir -p ~/.kaggle
cp kaggle.json ~/.kaggle/
chmod 600 ~/.kaggle/kaggle.json
```

---

## Usage

**Full pipeline:**
```bash
nextflow run main.nf -profile cdac -resume
```

**Individual modules:**
```bash
# QC only
python3 modules/qc/qc.py data/visium_breast/V1_*.h5 data/visium_breast/spatial results/qc

# Normalization
python3 modules/normalize/normalize.py results/qc/filtered_qc.h5ad results/normalize

# Clustering
python3 modules/cluster/cluster.py results/normalize/normalized.h5ad results/cluster

# TME Analysis
python3 modules/tme/tme.py results/cluster/clustered.h5ad results/deconv/cell_proportions.csv results/tme
```

---

## Project Structure

```
SpatialBridge/
├── main.nf                         # Nextflow DSL2 pipeline
├── nextflow.config                 # HPC resource configuration
├── modules/
│   ├── qc/
│   │   ├── qc.py                  # Adaptive MAD-based QC
│   │   └── qc.nf                  # Nextflow QC module
│   ├── normalize/
│   │   └── normalize.py           # Scran-style normalization
│   ├── cluster/
│   │   └── cluster.py             # Leiden clustering + UMAP
│   ├── deconv/
│   │   └── deconv.py              # Cell2location deconvolution
│   └── tme/
│       └── tme.py                 # TME + Moran's I analysis
└── kaggle_deconv/
    ├── deconv_kaggle.ipynb        # GPU deconvolution notebook
    └── kernel-metadata.json       # Kaggle kernel config
```

---

## HPC + Kaggle GPU Architecture

A key design challenge in this pipeline is that Cell2location requires GPU for stable Bayesian inference (CPU execution causes numerical precision failures in GammaPoisson distributions). Since CDAC HPC does not provide GPU allocation by default, SpatialBridge implements an automated Kaggle API workflow:

```
CDAC HPC
    ↓
kaggle kernels push     ← Submit GPU job
    ↓
Kaggle T4 x2 GPU        ← Cell2location trains (200+300 epochs)
    ↓
kaggle kernels output   ← Auto-fetch results back to HPC
    ↓
CDAC HPC continues      ← TME analysis with deconvolved data
```

This pattern enables GPU-dependent steps within an otherwise CPU-only HPC environment — a novel infrastructure approach not implemented in existing spatial transcriptomics pipelines.

---

## Biological Interpretation

The pipeline was validated on human breast cancer FFPE tissue. Key findings:

- **Cancer Epithelial cells** dominate 74.9% of tissue spots, consistent with high tumor burden
- **T-cell infiltration** (88 spots, abundance 0.797) indicates active immune response — "inflamed tumor" phenotype associated with immunotherapy response
- **High Moran's I** (>0.77) for tumor and stromal cells confirms spatially structured tissue organization rather than random cell distribution

---

## Limitations and Future Work

- Currently supports **Visium** platform only; Xenium and CosMx support planned
- Deconvolution requires GPU (Kaggle API used as workaround)
- Docker/Singularity containerization in progress for full reproducibility

---

## References

1. Cable et al. (2022). Robust decomposition of cell type mixtures in spatial transcriptomics. *Nature Biotechnology.*
2. Wu et al. (2021). A single-cell and spatially resolved atlas of human breast cancers. *Nature Genetics.*
3. Palla et al. (2022). Squidpy: a scalable framework for spatial omics analysis. *Nature Methods.*
4. Lun et al. (2016). Pooling across cells to normalize single-cell RNA sequencing data. *Genome Biology.*

---

## Author

**Tushar Parihar**
MSc Bioinformatics, Savitribai Phule Pune University (SPPU)
CDAC ICE Cloud | Spatial Transcriptomics | Pipeline Engineering
