# ============================================================
# RNA-seq Analysis Pipeline — Core Function
# Author: Santi Isgrò
# Dataset: Dantas Machado et al. 2022
# Description: Reusable pipeline for differential expression
#              analysis and functional enrichment (GO/KEGG)
# ============================================================

#' Run RNA-seq Differential Expression and Enrichment Pipeline
#'
#' Performs DESeq2-based differential expression analysis followed by
#' GO (Biological Process) and KEGG pathway enrichment analysis.
#' Returns a named list of results and plots for downstream use.
#'
#' @param dds A DESeqDataSet object (already fitted via DESeq()).
#' @param group1 Character. Name of the first condition (numerator in contrast).
#' @param group2 Character. Name of the second condition (denominator in contrast).
#' @param padj_cutoff Numeric. Adjusted p-value threshold (default: 0.05).
#' @param lfc_cutoff Numeric. Absolute log2 fold-change threshold (default: 1).
#'
#' @return A named list containing:
#'   \item{results_full}{Full DESeq2 results object}
#'   \item{results_significant}{Filtered significant DEGs}
#'   \item{upregulated}{Upregulated DEGs (LFC > lfc_cutoff)}
#'   \item{downregulated}{Downregulated DEGs (LFC < -lfc_cutoff)}
#'   \item{go_enrichment}{GO enrichment results (clusterProfiler)}
#'   \item{kegg_enrichment}{KEGG enrichment results (clusterProfiler)}
#'   \item{go_plot}{GO dotplot (ggplot2 object)}
#'   \item{kegg_plot}{KEGG dotplot (ggplot2 object)}

run_rnaseq_pipeline <- function(dds,
                                 group1,
                                 group2,
                                 padj_cutoff = 0.05,
                                 lfc_cutoff = 1) {

  # ---- STEP 1: DIFFERENTIAL EXPRESSION ANALYSIS ----
  cat("Running DESeq2 comparison:", group1, "vs", group2, "\n")

  res  <- results(dds, contrast = c("condition", group1, group2))
  sig  <- subset(res, padj < padj_cutoff)
  up   <- subset(sig, log2FoldChange >  lfc_cutoff)
  down <- subset(sig, log2FoldChange < -lfc_cutoff)

  cat("Total DEGs:     ", nrow(sig),  "\n")
  cat("Upregulated:    ", nrow(up),   "\n")
  cat("Downregulated:  ", nrow(down), "\n\n")

  # ---- STEP 2: GO ENRICHMENT (Biological Process) ----
  cat("Running GO enrichment...\n")

  # Strip version suffixes from Ensembl IDs (e.g. ENSMUSG00000001.5 -> ENSMUSG00000001)
  gene_list  <- sub("\\..*", "", rownames(sig))
  gene_entrez <- bitr(gene_list,
                      fromType = "ENSEMBL",
                      toType   = "ENTREZID",
                      OrgDb    = org.Mm.eg.db)

  go_res <- enrichGO(gene          = gene_entrez$ENTREZID,
                     OrgDb         = org.Mm.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = padj_cutoff,
                     readable      = TRUE)

  cat("GO terms enriched:",
      nrow(go_res@result[go_res@result$p.adjust < padj_cutoff, ]), "\n\n")

  # ---- STEP 3: KEGG PATHWAY ENRICHMENT ----
  cat("Running KEGG enrichment...\n")

  kegg_res <- enrichKEGG(gene         = gene_entrez$ENTREZID,
                          organism     = "mmu",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = padj_cutoff)

  cat("KEGG pathways enriched:",
      nrow(kegg_res@result[kegg_res@result$p.adjust < padj_cutoff, ]), "\n\n")

  # ---- STEP 4: VISUALIZATION ----
  cat("Generating plots...\n")

  go_plot <- dotplot(go_res,
                     showCategory = 15,
                     title = paste("GO Enrichment —", group1, "vs", group2))

  kegg_plot <- dotplot(kegg_res,
                       showCategory = 15,
                       title = paste("KEGG Enrichment —", group1, "vs", group2))

  print(go_plot)
  print(kegg_plot)

  # ---- RETURN ----
  return(list(
    results_full        = res,
    results_significant = sig,
    upregulated         = up,
    downregulated       = down,
    go_enrichment       = go_res,
    kegg_enrichment     = kegg_res,
    go_plot             = go_plot,
    kegg_plot           = kegg_plot
  ))
}
