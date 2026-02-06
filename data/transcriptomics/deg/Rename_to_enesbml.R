if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")

library(biomaRt)
setwd('C:/Users/pleas/YandexDisk/work_with_genes/deg')
# Подключение к Ensembl
ensembl <- useMart("ensembl")

# Выбираем dataset для курицы
ensembl <- useDataset("ggallus_gene_ensembl", mart = ensembl)

# Загрузка данных из CSV файла
data <- read.csv("filt_ggrsw1_degs_v7.csv")
head(data)

# Выведем колонку с символами генов
gene_symbols <- data$Genes
head(gene_symbols)

# Конвертация символов генов в Ensembl IDs
gene_info <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                   filters = 'external_gene_name',
                   values = gene_symbols,
                   mart = ensembl)

# Проверим результат
head(gene_info)

# Объединение таблиц по символам генов
merged_data <- merge(data, gene_info, by.x = "Genes", by.y = "external_gene_name", all.x = TRUE)

# Сохраним обновленный файл с добавленными Ensembl IDs
write.csv(merged_data, "filt_ggrsw1_degs_v7_with_ensembl.csv", row.names = FALSE)

# Выведем первые строки для проверки
head(merged_data)
