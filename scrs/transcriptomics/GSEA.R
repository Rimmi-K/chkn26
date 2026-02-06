setwd("/home/rimmi/WorkStation/drvoms-chkn/work_with_genes/gsea/")

library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(GOSemSim)
library(ggplot2)
library(readr)
library(org.Gg.eg.db)
library(BiocParallel)
register(SerialParam())
library(ggtree)
library(tidyverse)
library(patchwork)
library(openxlsx)

# helper: достать колонку с именами генов из CSV, где они могли быть в первом столбце
.extract_gene_col <- function(df) {
  cn <- colnames(df)
  # 1) если уже есть колонка gene / Gene / SYMBOL
  cand <- cn[grepl("^gene$|^Gene$|^SYMBOL$", cn)]
  if (length(cand) >= 1) return(df[[cand[1]]])
  
  # 2) частый случай после write.csv(row.names=TRUE): X1, ...1
  cand <- cn[grepl("^X1$|^...1$", cn, perl = TRUE)]
  if (length(cand) >= 1) return(df[[cand[1]]])
  
  # 3) если пришло как data.frame с rownames
  rn <- rownames(df)
  if (!is.null(rn) && length(rn) == nrow(df)) return(rn)
  
  stop("Не нашёл колонку с именами генов: ожидаю 'gene' или первый столбец (X1/…1), либо rownames.")
}

run_gsea_and_save <- function(tissue, PVAL_CUTOFF = 0.05, SHOW_CATEGORIES = 15, ONTOLOGY = "BP") {
  message("== ", tissue, " ==")
  
  # 1. Загрузка DESeq2 результатов
  deg_file <- paste0("../deg/archive/ggrsw1_v01/degs_results_", tissue, ".csv")
  deg_unfiltered <- suppressMessages(read_csv(deg_file, show_col_types = FALSE)) |> as.data.frame()
  
  # 2. Подготовка geneList
  genes <- .extract_gene_col(deg_unfiltered)
  if (!"log2FoldChange" %in% colnames(deg_unfiltered)) {
    stop("В файле нет колонки 'log2FoldChange'. Проверь путь и формат CSV.")
  }
  
  df <- tibble(
    gene = genes,
    log2FC = deg_unfiltered$log2FoldChange
  ) |>
    filter(!is.na(gene), !is.na(log2FC), is.finite(log2FC))
  
  # Агрегация дублей (взять максимум по модулю)
  df <- df |>
    group_by(gene) |>
    summarise(log2FC = log2FC[which.max(abs(log2FC))], .groups = "drop")
  
  # Именованный вектор, отсортированный по убыванию
  gene_list <- df$log2FC
  names(gene_list) <- df$gene
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  message("Генов в geneList: ", length(gene_list))
  
  if (length(gene_list) < 100) {
    warning("Мало генов для GSEA (", length(gene_list), "). Результаты могут быть ненадёжными.")
  }
  
  # 3. GSEA (GO) с расширенными границами
  gse_go <- gseGO(
    geneList      = gene_list,
    OrgDb         = org.Gg.eg.db,
    ont           = ONTOLOGY,
    keyType       = "SYMBOL",
    minGSSize     = 5,                # включаем маленькие пути
    maxGSSize     = 1000,             # включаем большие пути
    pvalueCutoff  = PVAL_CUTOFF,
    pAdjustMethod = "BH",
    verbose       = FALSE
  )
  
  if (is.null(gse_go) || nrow(as.data.frame(gse_go)) == 0) {
    warning("GSEA вернула пустой результат для ", tissue, ".")
    return(NULL)
  }
  
  message("Найдено значимых путей: ", nrow(as.data.frame(gse_go)))
  
  # 4. Семантическое сходство
  gse_sim <- tryCatch(
    pairwise_termsim(gse_go), 
    error = function(e) { 
      message("pairwise_termsim предупредил: ", e$message)
      gse_go 
    }
  )
  
  # 5. Сохранение результатов
  res_df <- as.data.frame(gse_go) |>
    arrange(pvalue)  # лучше сортировать по p-value, а не NES
  write.xlsx(res_df, paste0(tissue, "_gsea_results.xlsx"), overwrite = TRUE)
  
  # 6. Визуализация
  p_tree <- treeplot(gse_sim, showCategory = SHOW_CATEGORIES) +
    ggtitle(paste("GSEA GO", ONTOLOGY, "-", tissue)) +
    theme(legend.position = "bottom")
  
  p_dot <- dotplot(gse_sim, showCategory = SHOW_CATEGORIES, split = ".sign") +
    facet_grid(. ~ .sign) +
    ggtitle(paste("Dotplot -", tissue))
  
  ggsave(paste0(tissue, "_gsea_treeplot1.pdf"), p_tree, width = 12, height = 8)
  ggsave(paste0(tissue, "_gsea_dotplot1.pdf"),  p_dot,  width = 12, height = 8)
  
  message("Готово: ", tissue)
  list(gse_result = gse_go, sim_result = gse_sim, tree_plot = p_tree, dot_plot = p_dot)
}

# Теперь можно вернуться к нормальным p-value cutoff
breast <- run_gsea_and_save("breast", 0.10, 15, "BP")
leg    <- run_gsea_and_save("leg",    0.30, 15, "BP") 
liver  <- run_gsea_and_save("liver",  0.30, 15, "BP") 
