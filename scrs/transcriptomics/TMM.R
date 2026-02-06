library(edgeR)
setwd('C:/Users/pleas/YandexDisk/work_with_genes')

# Указываем путь к файлу с коунтами
counts_file <- "counts/counts_ggrsw1_v01_v2.csv"
output_dir <- "TMM/ggrsw1_v01/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Загружаем данные
counts_data <- read.csv(counts_file, row.names = 1, check.names = FALSE)

# Уникализируем имена столбцов
colnames(counts_data) <- make.unique(colnames(counts_data))

# Создаем таблицу с метаданными
sample_info <- data.frame(
  row.names = colnames(counts_data),
  growth = sub("-.*", "", colnames(counts_data)),
  tissue = sub(".*-", "", gsub("\\d+$", "", colnames(counts_data)))
)

sample_info$growth <- factor(sample_info$growth, levels = c("slow", "fast"))
sample_info$tissue <- factor(sample_info$tissue)
sample_info$tissue <- gsub("\\.$", "", sample_info$tissue)
tissues <- unique(sample_info$tissue)

# Функция для проверки согласованности и агрегации
check_and_aggregate <- function(data, threshold_cv = 15) {
  cv <- apply(data, 1, function(x) sd(x) / mean(x) * 100)
  aggregated_data <- data.frame(
    gene = rownames(data),
    mean = apply(data, 1, mean),
    median = apply(data, 1, median),
    cv = cv,
    stringsAsFactors = FALSE
  )
  if (any(cv > threshold_cv)) {
    warning(paste("Высокий CV (>", threshold_cv, "%) обнаружен для генов:", 
                  paste(rownames(data)[cv > threshold_cv], collapse = ", ")))
  }
  return(aggregated_data)
}

# Проходим по каждой ткани
for (tissue in tissues) {
  # Отбираем образцы для текущей ткани
  tissue_samples <- rownames(sample_info[sample_info$tissue == tissue, ])
  tissue_counts <- counts_data[, tissue_samples]
  tissue_sample_info <- sample_info[tissue_samples, ]
  
  # Создание DGEList-объекта
  dge <- DGEList(counts = tissue_counts, group = tissue_sample_info$growth)
  
  # Фильтрация низкоэкспрессированных генов
  keep <- rowSums(cpm(dge) > 0.5) >= 3
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  # TMM нормализация
  dge <- calcNormFactors(dge, method = "TMM")
  
  # Нормализованные counts
  norm_counts <- cpm(dge, normalized.lib.sizes = TRUE)
  
  # Сохранение полной TMM матрицы для всех образцов ткани
  norm_counts_df <- data.frame(genes = rownames(norm_counts), norm_counts, check.names = FALSE)
  write.table(norm_counts_df, 
              file = paste0(output_dir, "TMM_full_matrix_", tissue, ".csv"),
              quote = FALSE, sep = "\t", row.names = FALSE)
  
  # Разделение на fast и slow
  fast_samples <- colnames(norm_counts)[tissue_sample_info$growth == "fast"]
  slow_samples <- colnames(norm_counts)[tissue_sample_info$growth == "slow"]
  
  # Обработка fast
  if (length(fast_samples) > 0) {
    fast_counts <- norm_counts[, fast_samples, drop = FALSE]
    # Сохранение TMM матрицы для fast
    fast_counts_df <- data.frame(genes = rownames(fast_counts), fast_counts, check.names = FALSE)
    write.table(fast_counts_df,
                file = paste0(output_dir, "TMM_matrix_", tissue, "_fast.csv"),
                quote = FALSE, sep = "\t", row.names = FALSE)
    # Медианы
    fast_aggregated <- check_and_aggregate(fast_counts)
    median_data <- fast_aggregated[, c("gene", "median")]
    colnames(median_data)[1] <- "genes"
    write.table(median_data, 
                file = paste0(output_dir, "median_", tissue, "_fast.csv"),
                quote = FALSE, row.names = FALSE, sep = "\t")
  }
  
  # Обработка slow
  if (length(slow_samples) > 0) {
    slow_counts <- norm_counts[, slow_samples, drop = FALSE]
    # Сохранение TMM матрицы для slow
    slow_counts_df <- data.frame(genes = rownames(slow_counts), slow_counts, check.names = FALSE)
    write.table(slow_counts_df,
                file = paste0(output_dir, "TMM_matrix_", tissue, "_slow.csv"),
                quote = FALSE, sep = "\t", row.names = FALSE)
    # Медианы
    slow_aggregated <- check_and_aggregate(slow_counts)
    median_data <- slow_aggregated[, c("gene", "median")]
    colnames(median_data)[1] <- "genes"
    write.table(median_data, 
                file = paste0(output_dir, "median_", tissue, "_slow.csv"),
                quote = FALSE, row.names = FALSE, sep = "\t")
  }
  
  # Сохранение факторов нормализации
  samples_df <- data.frame(genes = rownames(dge$samples), dge$samples, check.names = FALSE)
  write.csv(samples_df, 
            file = paste0(output_dir, "norm_factors_", tissue, ".csv"),
            row.names = FALSE)
}