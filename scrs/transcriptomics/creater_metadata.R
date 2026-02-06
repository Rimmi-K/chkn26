# Загружаем файл с коунтами
counts_file <- "counts/counts_ggrsw1_v01_v2.csv" # Укажите ваш файл
counts_data <- read.csv(counts_file, row.names = 1)

# Создаем метаданные
sample_names <- colnames(counts_data)

# Разбиваем названия образцов
sample_info <- data.frame(
  sample = sample_names,
  growth = sub("-.*", "", sample_names),          # Извлекаем "fast" или "slow"
  tissue = sub(".*-", "", gsub("\\d+$", "", sample_names)) # Извлекаем "breast", "leg", и т.д.
)

# Устанавливаем sample как rownames
rownames(sample_info) <- sample_info$sample

# Выводим таблицу метаданных
print(sample_info)