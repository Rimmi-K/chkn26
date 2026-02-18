library(edgeR)

COUNTS_FILE <- "counts/counts_ggrsw1_v01_v2.csv"
OUTPUT_DIR  <- "TMM/ggrsw1_v01/"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

counts_data <- read.csv(COUNTS_FILE, row.names = 1, check.names = FALSE)
colnames(counts_data) <- make.unique(colnames(counts_data))

sample_info <- data.frame(
  row.names = colnames(counts_data),
  growth    = sub("-.*", "", colnames(counts_data)),
  tissue    = gsub("\\.$", "", sub(".*-", "", gsub("\\d+$", "", colnames(counts_data))))
)
sample_info$growth <- factor(sample_info$growth, levels = c("slow", "fast"))
sample_info$tissue <- factor(sample_info$tissue)

write_group_matrix <- function(norm_counts, group_samples, tissue, group_label, output_dir) {
  if (length(group_samples) == 0) return(invisible(NULL))

  mat <- norm_counts[, group_samples, drop = FALSE]

  write.table(
    data.frame(genes = rownames(mat), mat, check.names = FALSE),
    file  = paste0(output_dir, "TMM_matrix_", tissue, "_", group_label, ".csv"),
    quote = FALSE, sep = "\t", row.names = FALSE
  )

  cv <- apply(mat, 1, function(x) sd(x) / mean(x) * 100)
  if (any(cv > 15, na.rm = TRUE)) {
    warning(sprintf("High CV (>15%%) in %s %s: %s",
                    tissue, group_label,
                    paste(rownames(mat)[cv > 15], collapse = ", ")))
  }

  median_df <- data.frame(genes = rownames(mat), median = apply(mat, 1, median))
  write.table(median_df,
              file  = paste0(output_dir, "median_", tissue, "_", group_label, ".csv"),
              quote = FALSE, sep = "\t", row.names = FALSE)
}

for (tissue in levels(sample_info$tissue)) {
  tissue_samples <- rownames(sample_info[sample_info$tissue == tissue, ])
  tissue_info    <- sample_info[tissue_samples, ]

  dge <- DGEList(counts = counts_data[, tissue_samples], group = tissue_info$growth)
  dge <- dge[rowSums(cpm(dge) > 0.5) >= 3, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")

  norm_counts <- cpm(dge, normalized.lib.sizes = TRUE)

  write.table(
    data.frame(genes = rownames(norm_counts), norm_counts, check.names = FALSE),
    file  = paste0(OUTPUT_DIR, "TMM_full_matrix_", tissue, ".csv"),
    quote = FALSE, sep = "\t", row.names = FALSE
  )

  fast_samples <- colnames(norm_counts)[tissue_info$growth == "fast"]
  slow_samples <- colnames(norm_counts)[tissue_info$growth == "slow"]

  write_group_matrix(norm_counts, fast_samples, tissue, "fast", OUTPUT_DIR)
  write_group_matrix(norm_counts, slow_samples, tissue, "slow", OUTPUT_DIR)

  write.csv(
    data.frame(genes = rownames(dge$samples), dge$samples, check.names = FALSE),
    file      = paste0(OUTPUT_DIR, "norm_factors_", tissue, ".csv"),
    row.names = FALSE
  )
}
