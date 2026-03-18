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

run_gsea_and_save <- function(tissue, pval_cutoff = 0.05, show_categories = 15,
                              ontology = "BP", minGSSize = 5, nPermSimple = 1000,
                              rank_metric = "stat") {

  message("== ", tissue, " (all genes, rank: ", rank_metric, ") ==")

  # Load DEG results
  deg_file <- paste0("data/transcriptomics/degs/degs_results_", tissue, ".csv")
  deg_df   <- suppressMessages(read_csv(deg_file, show_col_types = FALSE)) |> as.data.frame()

  if (!"log2FoldChange" %in% colnames(deg_df)) stop("Column 'log2FoldChange' not found.")

  # Calculate ranking metric
  if (rank_metric == "log2FC") {
    gene_list <- tibble(
      gene   = .extract_gene_col(deg_df),
      rank_value = deg_df$log2FoldChange
    )
  } else if (rank_metric == "stat") {
    if (!"stat" %in% colnames(deg_df)) stop("Column 'stat' not found for Wald statistic ranking.")
    gene_list <- tibble(
      gene   = .extract_gene_col(deg_df),
      rank_value = deg_df$stat
    )
  } else if (rank_metric == "signed_pvalue") {
    if (!"pvalue" %in% colnames(deg_df)) stop("Column 'pvalue' not found.")
    gene_list <- tibble(
      gene   = .extract_gene_col(deg_df),
      rank_value = sign(deg_df$log2FoldChange) * -log10(deg_df$pvalue + 1e-300)
    )
  } else {
    stop("Unknown rank_metric: ", rank_metric)
  }

  # Filter and deduplicate
  gene_list <- gene_list |>
    filter(!is.na(gene), is.finite(rank_value)) |>
    group_by(gene) |>
    summarise(rank_value = rank_value[which.max(abs(rank_value))], .groups = "drop") |>
    arrange(desc(rank_value)) |>
    deframe()

  message("Genes in geneList: ", length(gene_list))
  if (length(gene_list) < 100) warning("Few genes in geneList (", length(gene_list), ").")

  gse_go <- gseGO(
    geneList      = gene_list,
    OrgDb         = org.Gg.eg.db,
    ont           = ontology,
    keyType       = "SYMBOL",
    minGSSize     = minGSSize,
    maxGSSize     = 1000,
    pvalueCutoff  = pval_cutoff,
    pAdjustMethod = "BH",
    nPermSimple   = nPermSimple,
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

  # Save results to results/transcriptomics/gsea/
  output_dir <- "results/transcriptomics/gsea"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  output_suffix <- ""
  if (rank_metric != "log2FC") {
    output_suffix <- paste0("_", rank_metric)
  }

  # Save as CSV with ontology metadata
  results_df <- as.data.frame(gse_go) |> arrange(pvalue)
  results_df$ontology <- ontology  # Add ontology column
  write_csv(
    results_df,
    file.path(output_dir, paste0(tissue, "_gsea_results", output_suffix, ".csv"))
  )

  n_pathways <- nrow(as.data.frame(gse_go))
  actual_show_categories <- min(show_categories, n_pathways)

  # Only create treeplot if we have enough pathways (need at least 3 for clustering)
  plot_suffix <- if (rank_metric != "log2FC") paste0("_", rank_metric) else ""

  if (n_pathways >= 3) {
    p_tree <- treeplot(gse_sim, showCategory = actual_show_categories) +
      ggtitle(paste("GSEA GO", ontology, "-", tissue, "(rank:", rank_metric, ")")) +
      theme(legend.position = "bottom")
    ggsave(file.path(output_dir, paste0(tissue, "_gsea_treeplot", plot_suffix, ".svg")),
           p_tree, width = 12, height = 8)
  } else {
    message("Skipping treeplot for ", tissue, " (only ", n_pathways, " pathways found, need at least 3)")
    p_tree <- NULL
  }

  p_dot <- dotplot(gse_sim, showCategory = actual_show_categories, split = ".sign") +
    facet_grid(. ~ .sign) +
    ggtitle(paste("Dotplot -", tissue, "(rank:", rank_metric, ")"))

  ggsave(file.path(output_dir, paste0(tissue, "_gsea_dotplot", plot_suffix, ".svg")),
         p_dot, width = 12, height = 8)

  message("Done: ", tissue)
  list(gse_result = gse_go, sim_result = gse_sim, tree_plot = p_tree, dot_plot = p_dot)
}




breast <- run_gsea_and_save("breast", pval_cutoff = 0.25, minGSSize = 15, nPermSimple = 100000,
                            rank_metric = "stat")

leg    <- run_gsea_and_save("leg",    pval_cutoff = 0.25, minGSSize = 15, nPermSimple = 100000,
                            rank_metric = "stat")

liver  <- run_gsea_and_save("liver",  pval_cutoff = 0.25, minGSSize = 15, nPermSimple = 100000,
                            rank_metric = "stat")


