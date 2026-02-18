counts_data <- read.csv("counts/counts_ggrsw1_v01_v2.csv", row.names = 1)

sample_names <- colnames(counts_data)

sample_info <- data.frame(
  row.names = sample_names,
  sample    = sample_names,
  growth    = sub("-.*",  "", sample_names),
  tissue    = gsub("\\d+$", "", sub(".*-", "", sample_names))
)

print(sample_info)
