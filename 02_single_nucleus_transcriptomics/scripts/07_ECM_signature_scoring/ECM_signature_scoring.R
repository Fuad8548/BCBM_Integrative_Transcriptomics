# ============================================================
# 06_ECM_signature_scoring_v2.R
# Corrected sample-level ECM signature scoring
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

input_file <- "results/ECM_pseudobulk/ECM_pseudobulk_logCPM_long.csv"

output_dir <- "results/ECM_signature_scoring_v2"
figure_dir <- "figures/ECM_signature_scoring_v2"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Load pseudobulk data
# ------------------------------------------------------------

pb <- read.csv(input_file, stringsAsFactors = FALSE)

required_cols <- c("sample", "cell_type", "gene", "logCPM")

missing_cols <- setdiff(required_cols, colnames(pb))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

pb <- pb %>%
  select(sample, cell_type, gene, logCPM) %>%
  mutate(
    sample = as.character(sample),
    cell_type = as.character(cell_type),
    gene = as.character(gene),
    logCPM = as.numeric(logCPM)
  )

# ------------------------------------------------------------
# 3. Integrity checks
# ------------------------------------------------------------

duplicate_records <- pb %>%
  count(sample, cell_type, gene) %>%
  filter(n > 1)

if (nrow(duplicate_records) > 0) {
  write.csv(
    duplicate_records,
    file.path(output_dir, "duplicate_records.csv"),
    row.names = FALSE
  )

  stop("Duplicate sample-cell_type-gene records detected.")
}

cat("Samples:", n_distinct(pb$sample), "\n")
cat("Cell types:", n_distinct(pb$cell_type), "\n")
cat("Genes:", n_distinct(pb$gene), "\n")

# ------------------------------------------------------------
# 4. Define ECM signatures
# ------------------------------------------------------------

signatures <- list(
  Core_stromal_ECM = c(
    "COL6A3",
    "AEBP1",
    "SULF1",
    "THBS2",
    "LUM",
    "POSTN",
    "CTSK"
  ),

  Vascular_ECM = c(
    "COL15A1",
    "CXCL12",
    "MMP2"
  ),

  Broad_ECM = c(
    "AEBP1",
    "COL14A1",
    "COL15A1",
    "COL6A3",
    "CTSK",
    "CXCL12",
    "DPT",
    "LUM",
    "MMP2",
    "POSTN",
    "SULF1",
    "THBS2",
    "THBS4"
  )
)

signature_table <- bind_rows(
  lapply(names(signatures), function(sig) {
    data.frame(
      signature = sig,
      gene = signatures[[sig]]
    )
  })
)

write.csv(
  signature_table,
  file.path(output_dir, "ECM_signature_gene_sets.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. Check signature gene availability
# ------------------------------------------------------------

gene_availability <- pb %>%
  distinct(cell_type, gene) %>%
  inner_join(signature_table, by = "gene") %>%
  count(cell_type, signature, name = "genes_available")

gene_availability <- gene_availability %>%
  left_join(
    signature_table %>%
      count(signature, name = "genes_in_signature"),
    by = "signature"
  ) %>%
  mutate(
    genes_missing = genes_in_signature - genes_available
  )

write.csv(
  gene_availability,
  file.path(output_dir, "ECM_signature_gene_availability.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Gene-wise Z-score within each cell type
#    IMPORTANT: calculate across samples, not across genes
# ------------------------------------------------------------

gene_z <- pb %>%
  group_by(cell_type, gene) %>%
  mutate(
    gene_mean = mean(logCPM, na.rm = TRUE),
    gene_sd = sd(logCPM, na.rm = TRUE),

    gene_z = ifelse(
      is.na(gene_sd) | gene_sd == 0,
      NA_real_,
      (logCPM - gene_mean) / gene_sd
    )
  ) %>%
  ungroup()

# Save gene-level Z scores for auditing
write.csv(
  gene_z,
  file.path(output_dir, "ECM_gene_level_Zscores.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Attach signature membership
# ------------------------------------------------------------

gene_z_sig <- gene_z %>%
  inner_join(signature_table, by = "gene")

# ------------------------------------------------------------
# 8. Calculate sample-level signature scores
#    Each score = mean gene-wise Z-score
#    for each sample x cell type x signature
# ------------------------------------------------------------

signature_scores <- gene_z_sig %>%
  group_by(sample, cell_type, signature) %>%
  summarise(
    signature_score = mean(gene_z, na.rm = TRUE),

    genes_expected = n_distinct(gene),

    genes_scored = sum(!is.na(gene_z)),

    .groups = "drop"
  ) %>%
  mutate(
    signature_score = ifelse(
      genes_scored == 0,
      NA_real_,
      signature_score
    )
  )

# ------------------------------------------------------------
# 9. CRITICAL VALIDATION
#    Scores should vary across samples where gene expression varies
# ------------------------------------------------------------

score_validation <- signature_scores %>%
  group_by(cell_type, signature) %>%
  summarise(
    n_samples = n(),
    n_unique_scores = n_distinct(round(signature_score, 8)),
    score_SD = sd(signature_score, na.rm = TRUE),
    score_min = min(signature_score, na.rm = TRUE),
    score_max = max(signature_score, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  score_validation,
  file.path(output_dir, "ECM_signature_score_validation.csv"),
  row.names = FALSE
)

print(score_validation)

# Warn if scores collapse to constants
collapsed_scores <- score_validation %>%
  filter(n_unique_scores <= 1 | is.na(score_SD) | score_SD == 0)

if (nrow(collapsed_scores) > 0) {
  warning(
    "Some signature scores remain constant. ",
    "Inspect gene-level Z scores and signature gene availability."
  )
}

# ------------------------------------------------------------
# 10. Summary statistics
# ------------------------------------------------------------

score_summary <- signature_scores %>%
  group_by(cell_type, signature) %>%
  summarise(
    n_samples = sum(!is.na(signature_score)),
    mean_score = mean(signature_score, na.rm = TRUE),
    median_score = median(signature_score, na.rm = TRUE),
    SD_score = sd(signature_score, na.rm = TRUE),
    min_score = min(signature_score, na.rm = TRUE),
    max_score = max(signature_score, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  score_summary,
  file.path(output_dir, "ECM_signature_score_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. Sample stratification: ECM-high vs ECM-low
#     Median split within each cell type/signature
# ------------------------------------------------------------

signature_stratification <- signature_scores %>%
  group_by(cell_type, signature) %>%
  mutate(
    cohort_median = median(signature_score, na.rm = TRUE),

    ECM_group = case_when(
      is.na(signature_score) ~ NA_character_,
      signature_score > cohort_median ~ "ECM_high",
      signature_score < cohort_median ~ "ECM_low",
      TRUE ~ "Median"
    )
  ) %>%
  ungroup()

write.csv(
  signature_stratification,
  file.path(output_dir, "ECM_high_low_sample_stratification.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 12. Stromal Core ECM sample ranking
# ------------------------------------------------------------

stromal_ranking <- signature_scores %>%
  filter(
    cell_type == "Stromal cells",
    signature == "Core_stromal_ECM"
  ) %>%
  arrange(desc(signature_score))

write.csv(
  stromal_ranking,
  file.path(output_dir, "Stromal_ECM_sample_ranking.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 13. Heatmap
# ------------------------------------------------------------

heatmap_data <- signature_scores %>%
  mutate(
    sample = factor(sample, levels = unique(sample))
  )

p_heatmap <- ggplot(
  heatmap_data,
  aes(x = signature, y = sample, fill = signature_score)
) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_wrap(~cell_type, scales = "free_y") +
  scale_fill_gradient2(
    low = "#C77C72",
    mid = "white",
    high = "#4B3CA8",
    midpoint = 0,
    name = "Score"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank()
  ) +
  labs(
    title = "ECM signature scores across BCBM samples",
    subtitle = "Gene-wise Z-scores calculated across samples within each cell type",
    x = "ECM signature",
    y = "Sample"
  )

ggsave(
  file.path(figure_dir, "ECM_signature_score_heatmap.png"),
  p_heatmap,
  width = 14,
  height = 11,
  dpi = 300
)

# ------------------------------------------------------------
# 14. Score distributions
# ------------------------------------------------------------

p_distribution <- ggplot(
  signature_scores,
  aes(x = signature, y = signature_score)
) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1) +
  facet_wrap(~cell_type, scales = "free_y") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "ECM signature score distributions",
    x = "Signature",
    y = "Mean gene-wise Z-score"
  )

ggsave(
  file.path(figure_dir, "ECM_signature_score_distributions.png"),
  p_distribution,
  width = 13,
  height = 9,
  dpi = 300
)

# ------------------------------------------------------------
# 15. Export final tables
# ------------------------------------------------------------

write.csv(
  signature_scores,
  file.path(output_dir, "ECM_signature_scores_long.csv"),
  row.names = FALSE
)

signature_scores_wide <- signature_scores %>%
  pivot_wider(
    names_from = signature,
    values_from = signature_score
  )

write.csv(
  signature_scores_wide,
  file.path(output_dir, "ECM_signature_scores_wide.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 16. Final report
# ------------------------------------------------------------

cat("\n====================================\n")
cat("ECM SIGNATURE SCORING COMPLETE\n")
cat("====================================\n")

cat("Output directory:", output_dir, "\n")
cat("Figure directory:", figure_dir, "\n")

cat("\nSamples with variable scores:\n")

score_validation %>%
  select(cell_type, signature, n_unique_scores, score_SD) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\nStromal Core ECM ranking preview:\n")

head(stromal_ranking, 10) %>%
  as.data.frame() %>%
  print(row.names = FALSE)

cat("\nCheck ECM_signature_score_validation.csv before proceeding.\n")
