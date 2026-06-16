# ============================================================
# RNA-seq Analysis — Main Script
# Author: Santi Isgrò
# Dataset: Dantas Machado et al. 2022
# Description: Applies run_rnaseq_pipeline() across all three
#              dietary group comparisons (FA, FT, NA)
# ============================================================

# ---- DEPENDENCIES ----
library(DESeq2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)

# ---- LOAD PIPELINE FUNCTION ----
source("R/pipeline.R")

# ---- LOAD DATA ----
# dds: DESeqDataSet object fitted via DESeq()
# Expected to be pre-loaded or imported from a saved .rds file.
# Example: dds <- readRDS("data/dds_fitted.rds")

# ---- RUN PIPELINE: ALL THREE COMPARISONS ----
# Dietary groups:
#   FA = Fat diet with Akkermansia muciniphila
#   FT = Fat diet without Akkermansia (control)
#   NA = Normal diet with Akkermansia muciniphila

FA_NA <- run_rnaseq_pipeline(dds, group1 = "FA", group2 = "NA")
FT_NA <- run_rnaseq_pipeline(dds, group1 = "FT", group2 = "NA")
FA_FT <- run_rnaseq_pipeline(dds, group1 = "FA", group2 = "FT")

# ---- ACCESS RESULTS ----

# Significant DEGs
FA_NA$results_significant
FT_NA$results_significant
FA_FT$results_significant

# Enrichment plots
FA_NA$go_plot
FA_NA$kegg_plot

# Save plots to file
ggsave("results/plots/GO_FA_vs_NA.png",  plot = FA_NA$go_plot,   width = 10, height = 8)
ggsave("results/plots/GO_FT_vs_NA.png",  plot = FT_NA$go_plot,   width = 10, height = 8)
ggsave("results/plots/GO_FA_vs_FT.png",  plot = FA_FT$go_plot,   width = 10, height = 8)
ggsave("results/plots/KEGG_FA_vs_NA.png", plot = FA_NA$kegg_plot, width = 10, height = 8)
ggsave("results/plots/KEGG_FT_vs_NA.png", plot = FT_NA$kegg_plot, width = 10, height = 8)
ggsave("results/plots/KEGG_FA_vs_FT.png", plot = FA_FT$kegg_plot, width = 10, height = 8)
