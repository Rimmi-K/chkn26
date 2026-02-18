library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(GOSemSim)
library(ggplot2)
library(readr)
library(org.Gg.eg.db)
library(BiocParallel)
library(ggtree)
library(tidyverse)
library(patchwork)
library(openxlsx)

register(SerialParam())

.extract_gene_col <- function(df) {
  cn <- colnames(df)

  cand <- cn[grepl("^gene$|^Gene$|^SYMBOL$", cn)]
  if (length(cand) >= 1) return(df[[cand[1]]])

  cand <- cn[grepl("^X1$|^[.]{3}1$", cn, perl = TRUE)]
  if (length(cand) >= 1) return(df[[cand[1]]])

  rn <- rownames(df)
  if (!is.null(rn) && length(rn) == nrow(df)) return(rn)

  stop("Cannot identify gene column in data frame.")
}

run_gsea_and_save <- function(tissue, pval_cutoff = 0.05, show_categories = 15, ontology = "BP") {
  message("== ", tissue, " ==")

  deg_file <- paste0("../deg/archive/ggrsw1_v01/degs_results_", tissue, ".csv")
  deg_df   <- suppressMessages(read_csv(deg_file, show_col_types = FALSE)) |> as.data.frame()

  if (!"log2FoldChange" %in% colnames(deg_df)) stop("Column 'log2FoldChange' not found.")

  gene_list <- tibble(
    gene   = .extract_gene_col(deg_df),
    log2FC = deg_df$log2FoldChange
  ) |>
    filter(!is.na(gene), is.finite(log2FC)) |>
    group_by(gene) |>
    summarise(log2FC = log2FC[which.max(abs(log2FC))], .groups = "drop") |>
    arrange(desc(log2FC)) |>
    deframe()

  message("Genes in geneList: ", length(gene_list))
  if (length(gene_list) < 100) warning("Few genes in geneList (", length(gene_list), ").")

  gse_go <- gseGO(
    geneList      = gene_list,
    OrgDb         = org.Gg.eg.db,
    ont           = ontology,
    keyType       = "SYMBOL",
    minGSSize     = 5,
    maxGSSize     = 1000,
    pvalueCutoff  = pval_cutoff,
    pAdjustMethod = "BH",
    verbose       = FALSE
  )

  if (is.null(gse_go) || nrow(as.data.frame(gse_go)) == 0) {
    warning("GSEA returned empty result for ", tissue, ".")
    return(NULL)
  }

  message("Significant pathways: ", nrow(as.data.frame(gse_go)))

  gse_sim <- tryCatch(
    pairwise_termsim(gse_go),
    error = function(e) { message("pairwise_termsim error: ", e$message); gse_go }
  )

  write.xlsx(
    as.data.frame(gse_go) |> arrange(pvalue),
    paste0(tissue, "_gsea_results.xlsx"),
    overwrite = TRUE
  )

  p_tree <- treeplot(gse_sim, showCategory = show_categories) +
    ggtitle(paste("GSEA GO", ontology, "-", tissue)) +
    theme(legend.position = "bottom")

  p_dot <- dotplot(gse_sim, showCategory = show_categories, split = ".sign") +
    facet_grid(. ~ .sign) +
    ggtitle(paste("Dotplot -", tissue))

  ggsave(paste0(tissue, "_gsea_treeplot1.pdf"), p_tree, width = 12, height = 8)
  ggsave(paste0(tissue, "_gsea_dotplot1.pdf"),  p_dot,  width = 12, height = 8)

  message("Done: ", tissue)
  list(gse_result = gse_go, sim_result = gse_sim, tree_plot = p_tree, dot_plot = p_dot)
}

setwd("/home/rimmi/WorkStation/drvoms-chkn/work_with_genes/gsea/")

breast <- run_gsea_and_save("breast", pval_cutoff = 0.10)
leg    <- run_gsea_and_save("leg",    pval_cutoff = 0.30)
liver  <- run_gsea_and_save("liver",  pval_cutoff = 0.30)
