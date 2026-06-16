# Data

This folder contains input data files required to run the pipeline. Data is not tracked in this repository due to file size.

## Expected files

| File | Description |
|------|-------------|
| `dds_fitted.rds` | DESeqDataSet object, already fitted via `DESeq()`. Load with `readRDS("data/dds_fitted.rds")`. |

## Data source

Raw RNA-seq data is publicly available from:

> Dantas Machado, A.C. et al. (2022). *Diet and gut microbiota composition shape the ileal transcriptome in mice.*

Preprocessing steps applied before this pipeline:
1. Quality control — FastQC
2. Pseudoalignment — Kallisto
3. Import and summarization — tximport
4. Normalization and fitting — DESeq2 (`DESeq()`)
