library(DESeq2)
library(ggplot2)

setwd('C:/Users/pleas/YandexDisk/work_with_genes')
data <- read.csv('./counts/counts_from_peaks_softCPM.csv', row.names = 1, check.names = FALSE)

tissues <- c("leg", "liver", "breast")
dir.create("./deg", showWarnings = FALSE)

for (tissue in tissues) {
  # новые паттерны: <tissue>_f<digits> и <tissue>_s<digits>
  fast_cols <- grep(paste0("^", tissue, "_f\\d+$"), colnames(data), ignore.case = TRUE)
  slow_cols <- grep(paste0("^", tissue, "_s\\d+$"), colnames(data), ignore.case = TRUE)
  
  if (length(fast_cols) == 0 || length(slow_cols) == 0) {
    message(sprintf("Пропускаю '%s': нет fast или slow колонок (fast=%d, slow=%d).",
                    tissue, length(fast_cols), length(slow_cols)))
    next
  }
  
  counts <- data[, c(fast_cols, slow_cols), drop = FALSE]
  
  # страховки для DESeq2
  counts[is.na(counts)] <- 0
  counts <- as.matrix(round(counts))
  storage.mode(counts) <- "integer"
  
  condition <- factor(c(rep("fast", length(fast_cols)), rep("slow", length(slow_cols))),
                      levels = c("slow", "fast"))
  coldata <- data.frame(row.names = colnames(counts), condition = condition)
  
  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
  
  # лёгкая префильтрация низких сумм счётов
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  
  dds <- DESeq(dds, test = "Wald")
  
  res <- results(dds, contrast = c("condition", "fast", "slow"), alpha = 0.05)
  res_full <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]
  
  lfc_thr <- 1.0
  res_subset <- res_full[abs(res_full$log2FoldChange) > lfc_thr & res_full$padj < 0.05, ]
  
  write.csv(as.data.frame(res_full),
            sprintf("./deg/full_ggrsw1_degs_%s_v9.csv", tissue),
            row.names = TRUE)
  write.csv(as.data.frame(res_subset),
            sprintf("./deg/filt_ggrsw1_degs_%s_v9.csv", tissue),
            row.names = TRUE)
  
  sf <- data.frame(sample = colnames(dds), sizeFactor = sizeFactors(dds))
  write.csv(sf, sprintf("./deg/sizefactors_%s_v8.csv", tissue), row.names = FALSE)
  
  cat("\n==", tissue, "==\n",
      "samples: fast=", length(fast_cols), ", slow=", length(slow_cols), "\n",
      "DEG (FDR<0.05): ", sum(res_full$padj < 0.05), "\n",
      "DEG (|LFC|>", lfc_thr, " & FDR<0.05): ", nrow(res_subset), "\n", sep = "")
  
  # --- Volcano (надёжная версия) ---
  lfc_thr <- 1.0  # если выше не объявлен
  
  res_df <- as.data.frame(res)
  # Защита от 0 и NA в padj
  res_df$padj_safe <- pmax(res_df$padj, .Machine$double.eps)
  
  # Оставляем только строки, где координаты для графика корректны
  finite_idx <- is.finite(res_df$log2FoldChange) & is.finite(-log10(res_df$padj_safe))
  vdf <- res_df[finite_idx, ]
  
  # Метка значимости рассчитывается по тем же строкам
  vdf$significant <- ifelse(vdf$padj < 0.05 & abs(vdf$log2FoldChange) > lfc_thr,
                            "Significant", "Not Significant")
  
  # Подготовка финального датафрейма для ggplot
  volcano_df <- data.frame(
    log2FC = vdf$log2FoldChange,
    neglog10FDR = -log10(vdf$padj_safe),
    gene = rownames(vdf),
    significant = vdf$significant,
    stringsAsFactors = FALSE
  )
  
  p <- ggplot(volcano_df, aes(x = log2FC, y = neglog10FDR)) +
    geom_point(aes(color = significant), alpha = 0.6, size = 2) +
    scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
    theme_minimal() +
    labs(title = paste("Volcano Plot -", tissue),
         x = "log2 Fold Change",
         y = "-log10(FDR)") +
    theme(legend.title = element_blank())
  
  ggsave(sprintf("./deg/volcano_%s_v9.png", tissue), plot = p, width = 6, height = 4, dpi = 300)
  
}
