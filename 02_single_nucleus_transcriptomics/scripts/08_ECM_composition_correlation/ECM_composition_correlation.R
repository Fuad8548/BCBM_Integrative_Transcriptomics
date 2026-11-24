# ============================================================
# MILESTONE 8
# Stromal Core ECM Score vs Cellular Composition
# ============================================================

library(tidyverse)

score_file <- "results/ECM_signature_scoring/ECM_signature_scores_long.csv"

composition_file <- "results/cellular_composition/sample_cell_proportions.csv"

output_dir <- "results/ECM_composition_association"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------

scores <- read.csv(score_file)
composition <- read.csv(composition_file)

# ------------------------------------------------------------
# 2. Extract stromal Core ECM scores
# ------------------------------------------------------------

ecm_scores <- scores %>%
  filter(
    cell_type == "Stromal cells",
    signature == "Core_stromal_ECM"
  ) %>%
  select(sample, signature_score)

# Verify
cat("ECM score samples:", nrow(ecm_scores), "\n")
print(head(ecm_scores))

# ------------------------------------------------------------
# 3. Prepare composition table
# ------------------------------------------------------------

composition_clean <- composition %>%
  rename(sample = X) %>%
  select(
    sample,
    Stromal.cells,
    Pericytes,
    Endothelial.cells,
    Tumor.cells,
    Myeloid.cells,
    T.NK.cells
  ) %>%
  pivot_longer(
    cols = -sample,
    names_to = "cell_type",
    values_to = "proportion"
  )

# ------------------------------------------------------------
# 4. Join by sample ID
# ------------------------------------------------------------

merged_data <- ecm_scores %>%
  inner_join(composition_clean, by = "sample")

cat("Merged rows:", nrow(merged_data), "\n")
cat("Matched samples:", n_distinct(merged_data$sample), "\n")

stopifnot(nrow(merged_data) > 0)

write.csv(
  merged_data,
  file.path(output_dir, "Stromal_CoreECM_composition_merged.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. Spearman correlation
# ------------------------------------------------------------

cor_results <- merged_data %>%
  group_by(cell_type) %>%
  summarise(
    n = sum(is.finite(signature_score) & is.finite(proportion)),
    rho = cor(
      signature_score,
      proportion,
      method = "spearman",
      use = "complete.obs"
    ),
    p_value = cor.test(
      signature_score,
      proportion,
      method = "spearman",
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(FDR)

write.csv(
  cor_results,
  file.path(output_dir, "Stromal_CoreECM_composition_correlations.csv"),
  row.names = FALSE
)

print(cor_results)

# ------------------------------------------------------------
# 6. Correlation plot
# ------------------------------------------------------------

p_cor <- ggplot(
  cor_results,
  aes(x = reorder(cell_type, rho), y = rho)
) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  labs(
    title = "Stromal Core ECM Score vs Cellular Composition",
    x = "Cell population",
    y = "Spearman rho"
  )

ggsave(
  file.path(output_dir, "ECM_composition_correlation_plot.png"),
  p_cor,
  width = 8,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 7. Scatterplots
# ------------------------------------------------------------

for (cell in unique(merged_data$cell_type)) {
  plot_data <- merged_data %>%
    filter(cell_type == cell)

  rho <- cor(
    plot_data$signature_score,
    plot_data$proportion,
    method = "spearman"
  )

  p <- ggplot(
    plot_data,
    aes(x = proportion, y = signature_score)
  ) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", se = TRUE) +
    theme_bw() +
    labs(
      title = paste("Stromal Core ECM vs", cell),
      subtitle = paste("Spearman rho =", round(rho, 3)),
      x = "Cell proportion",
      y = "Stromal Core ECM score"
    )

  ggsave(
    file.path(
      output_dir,
      paste0("scatter_ECM_vs_", cell, ".png")
    ),
    p,
    width = 6,
    height = 5,
    dpi = 300
  )
}

cat("\nMilestone 8 analysis completed successfully.\n")
