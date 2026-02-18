library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(effsize)

growth_df <- read_csv("samples.csv", col_names = FALSE)
colnames(growth_df) <- c("id", "animal_id", "tissue", "growth_rate")

fast_ids <- unique(growth_df |> filter(growth_rate == "High growth rate") |> pull(animal_id))
slow_ids <- unique(growth_df |> filter(growth_rate == "Low growth rate")  |> pull(animal_id))

pheno_df <- read_csv("phenotype_base.csv", locale = locale(encoding = "CP1251")) |>
  mutate(
    `Ind. №` = as.numeric(`Ind. №`),
    growth_group = case_when(
      `Ind. №` %in% fast_ids ~ "fast",
      `Ind. №` %in% slow_ids ~ "slow",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(growth_group)) |>
  select(-`Ind. №`, -`№`)

pheno_long <- pheno_df |>
  select(where(is.numeric), growth_group) |>
  pivot_longer(-growth_group, names_to = "parameter", values_to = "value")

hist_plot <- ggplot(pheno_long, aes(x = value, fill = growth_group)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, position = "identity", alpha = 0.5) +
  geom_density(alpha = 0.7) +
  facet_wrap(~parameter, scales = "free", ncol = 3) +
  theme_minimal(base_size = 18) +
  theme(
    strip.text    = element_text(size = 20),
    axis.title    = element_text(size = 20),
    axis.text     = element_text(size = 16),
    legend.title  = element_text(size = 18),
    legend.text   = element_text(size = 16),
    plot.title    = element_text(size = 22)
  ) +
  labs(title = "Trait distributions by growth group", x = "Value", y = "Density", fill = "Growth group")

ggsave("hist_density_plot.png", hist_plot, width = 12, height = 8, dpi = 300)

wilcox_results <- pheno_long |>
  group_by(parameter) |>
  summarise(p_value = wilcox.test(value ~ growth_group)$p.value, .groups = "drop") |>
  mutate(p_adj = p.adjust(p_value, method = "BH"))

cliff_results <- pheno_long |>
  group_by(parameter) |>
  summarise(
    cliff_delta = cliff.delta(pull(pick(value)), pull(pick(growth_group)))$estimate,
    .groups = "drop"
  )

results <- left_join(wilcox_results, cliff_results, by = "parameter") |>
  arrange(p_adj)

signif_params <- results |> filter(p_adj < 0.05) |> pull(parameter)

box_plot <- pheno_long |>
  filter(parameter %in% signif_params) |>
  ggplot(aes(x = growth_group, y = value, fill = growth_group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  facet_wrap(~parameter, scales = "free_y", ncol = 4) +
  theme_minimal(base_size = 18) +
  theme(
    strip.text    = element_text(size = 20),
    axis.title    = element_text(size = 20),
    axis.text     = element_text(size = 16),
    axis.text.x   = element_text(angle = 45, hjust = 1),
    legend.title  = element_text(size = 18),
    legend.text   = element_text(size = 16),
    plot.title    = element_text(size = 22)
  ) +
  labs(title = "Significant parameters (padj < 0.05)", x = "Growth group", y = "Value", fill = "Growth group")

ggsave("boxplot_significant_params.pdf", box_plot, width = 15, height = 15, dpi = 300)

write_csv(results,  "results_statas.csv")
write_csv(pheno_df, "pheno_df.csv")
