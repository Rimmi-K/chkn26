library(DESeq2)
library(ggplot2)

data <- read.csv('./counts/counts_from_peaks_softCPM.csv', row.names = 1, check.names = FALSE)

tissues <- c("leg", "liver", "breast")
dir.create("./deg", showWarnings = FALSE)

LFC_THR <- 1.0

for (tissue in tissues) {
  fast_cols <- grep(paste0("^", tissue, "_f\\d+$"), colnames(data), ignore.case = TRUE)
  slow_cols <- grep(paste0("^", tissue, "_s\\d+$"), colnames(data), ignore.case = TRUE)

  if (length(fast_cols) == 0 || length(slow_cols) == 0) {
    message(sprintf("Skipping '%s': no fast or slow columns (fast=%d, slow=%d).",
                    tissue, length(fast_cols), length(slow_cols)))
    next
  }

  counts <- data[, c(fast_cols, slow_cols), drop = FALSE]
  counts[is.na(counts)] <- 0
  storage.mode(counts) <- "integer"
  counts <- as.matrix(round(counts))

  condition <- factor(
    c(rep("fast", length(fast_cols)), rep("slow", length(slow_cols))),
    levels = c("slow", "fast")
  )
  coldata <- data.frame(row.names = colnames(counts), condition = condition)

  dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~condition)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds, test = "Wald")

  res      <- results(dds, contrast = c("condition", "fast", "slow"), alpha = 0.05)
  res_full <- as.data.frame(res[!is.na(res$log2FoldChange) & !is.na(res$padj), ])
  res_filt <- res_full[abs(res_full$log2FoldChange) > LFC_THR & res_full$padj < 0.05, ]

  write.csv(res_full, sprintf("./deg/full_ggrsw1_degs_%s_v9.csv", tissue), row.names = TRUE)
  write.csv(res_filt, sprintf("./deg/filt_ggrsw1_degs_%s_v9.csv", tissue), row.names = TRUE)
  write.csv(
    data.frame(sample = colnames(dds), sizeFactor = sizeFactors(dds)),
    sprintf("./deg/sizefactors_%s_v8.csv", tissue),
    row.names = FALSE
  )

  cat(sprintf("\n== %s ==\n  fast=%d  slow=%d\n  DEG (FDR<0.05): %d\n  DEG (|LFC|>%g & FDR<0.05): %d\n",
              tissue, length(fast_cols), length(slow_cols),
              sum(res_full$padj < 0.05), LFC_THR, nrow(res_filt)))

  # Volcano plot
  vdf <- as.data.frame(res)
  vdf$padj_safe <- pmax(vdf$padj, .Machine$double.eps)
  vdf <- vdf[is.finite(vdf$log2FoldChange) & is.finite(-log10(vdf$padj_safe)), ]
  vdf$significant <- ifelse(vdf$padj < 0.05 & abs(vdf$log2FoldChange) > LFC_THR, "Significant", "Not Significant")

  p <- ggplot(vdf, aes(x = log2FoldChange, y = -log10(padj_safe))) +
    geom_point(aes(color = significant), alpha = 0.6, size = 2) +
    scale_color_manual(values = c("Not Significant" = "grey", "Significant" = "red")) +
    theme_minimal() +
    theme(legend.title = element_blank()) +
    labs(title = paste("Volcano -", tissue), x = "log2 Fold Change", y = "-log10(FDR)")

  ggsave(sprintf("./deg/volcano_%s_v9.png", tissue), plot = p, width = 6, height = 4, dpi = 300)
}
