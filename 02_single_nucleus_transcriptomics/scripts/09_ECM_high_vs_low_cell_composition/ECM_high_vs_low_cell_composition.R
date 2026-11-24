library(tidyverse)

# Load validated scores and composition
scores <- read.csv(
  "results/ECM_signature_scoring_v2/ECM_signature_scores_long.csv"
)

composition <- read.csv(
  "results/cellular_composition/sample_cell_proportions.csv"
) %>%
  rename(sample = X)

# Stromal Core ECM scores
stromal_core <- scores %>%
  filter(
    cell_type == "Stromal cells",
    signature == "Core_stromal_ECM"
  ) %>%
  select(sample, signature_score)

# Median-based grouping
median_score <- median(stromal_core$signature_score, na.rm = TRUE)

stromal_core <- stromal_core %>%
  mutate(
    ECM_group = if_else(
      signature_score >= median_score,
      "ECM_high",
      "ECM_low"
    )
  )

# Select cell populations
cell_data <- composition %>%
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

# Merge and compare
comparison <- stromal_core %>%
  select(sample, ECM_group) %>%
  inner_join(cell_data, by = "sample")

group_results <- comparison %>%
  group_by(cell_type) %>%
  summarise(
    n_low = sum(ECM_group == "ECM_low"),
    n_high = sum(ECM_group == "ECM_high"),
    median_low = median(proportion[ECM_group == "ECM_low"]),
    median_high = median(proportion[ECM_group == "ECM_high"]),
    p_value = wilcox.test(
      proportion ~ ECM_group,
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(FDR = p.adjust(p_value, method = "BH")) %>%
  arrange(FDR)

print(group_results)

# Save results
out_dir <- "results/ECM_composition_association"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  group_results,
  file.path(out_dir, "ECM_high_vs_low_cell_composition.csv"),
  row.names = FALSE
)


# 9B: Stromal ECM gene characterization

pb <- read.csv(
  "results/ECM_pseudobulk/ECM_pseudobulk_logCPM_long.csv"
)

ecm_genes <- c(
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

# Stromal-cell expression across samples
stromal_ecm <- pb %>%
  filter(
    cell_type == "Stromal cells",
    gene %in% ecm_genes
  )

gene_summary <- stromal_ecm %>%
  group_by(gene) %>%
  summarise(
    n_samples = n_distinct(sample),
    mean_logCPM = mean(logCPM, na.rm = TRUE),
    median_logCPM = median(logCPM, na.rm = TRUE),
    SD = sd(logCPM, na.rm = TRUE),
    min = min(logCPM, na.rm = TRUE),
    max = max(logCPM, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_logCPM))

print(gene_summary)

write.csv(
  gene_summary,
  "results/ECM_signature_scoring_v2/stromal_ECM_gene_summary.csv",
  row.names = FALSE
)


# Milestone 9C: Individual stromal ECM gene comparison

# 1. Load validated Stromal Core ECM scores
scores <- read.csv(
  "results/ECM_signature_scoring_v2/ECM_signature_scores_long.csv"
)

stromal_core <- scores %>%
  filter(
    cell_type == "Stromal cells",
    signature == "Core_stromal_ECM"
  ) %>%
  select(sample, signature_score)

# 2. Define ECM-high and ECM-low groups
median_score <- median(
  stromal_core$signature_score,
  na.rm = TRUE
)

stromal_core <- stromal_core %>%
  mutate(
    ECM_group = if_else(
      signature_score >= median_score,
      "ECM_high",
      "ECM_low"
    )
  )

print(table(stromal_core$ECM_group))

# 3. Load stromal pseudobulk expression
pb <- read.csv(
  "results/ECM_pseudobulk/ECM_pseudobulk_logCPM_long.csv"
)

ecm_genes <- c(
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

stromal_expression <- pb %>%
  filter(
    cell_type == "Stromal cells",
    gene %in% ecm_genes
  ) %>%
  select(sample, gene, logCPM)

# 4. Merge expression with group assignments
analysis_data <- stromal_expression %>%
  inner_join(
    stromal_core %>% select(sample, ECM_group),
    by = "sample"
  )

# Check sample coverage
print(
  analysis_data %>%
    distinct(sample, ECM_group) %>%
    count(ECM_group)
)

# 5. Compare each gene between groups
gene_comparison <- analysis_data %>%
  group_by(gene) %>%
  summarise(
    n_low = sum(ECM_group == "ECM_low" & !is.na(logCPM)),
    n_high = sum(ECM_group == "ECM_high" & !is.na(logCPM)),
    median_low = median(
      logCPM[ECM_group == "ECM_low"],
      na.rm = TRUE
    ),
    median_high = median(
      logCPM[ECM_group == "ECM_high"],
      na.rm = TRUE
    ),
    median_difference = median_high - median_low,
    p_value = wilcox.test(
      logCPM ~ ECM_group,
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    direction = case_when(
      median_difference > 0 ~ "Higher in ECM-high",
      median_difference < 0 ~ "Higher in ECM-low",
      TRUE ~ "No median difference"
    )
  ) %>%
  arrange(FDR)

print(gene_comparison, n = 13)

# 6. Save results
out_dir <- "results/ECM_composition_association"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  gene_comparison,
  file.path(
    out_dir,
    "stromal_ECM_genes_high_vs_low.csv"
  ),
  row.names = FALSE
)

# 7. Expression distribution plot
p <- ggplot(
  analysis_data,
  aes(x = ECM_group, y = logCPM, fill = ECM_group)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.65
  ) +
  geom_jitter(
    width = 0.12,
    size = 1.5,
    alpha = 0.7
  ) +
  facet_wrap(~gene, scales = "free_y", ncol = 4) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(
    title = "Stromal ECM gene expression: ECM-high versus ECM-low",
    x = NULL,
    y = "Pseudobulk logCPM"
  )

ggsave(
  file.path(
    out_dir,
    "stromal_ECM_genes_high_vs_low.png"
  ),
  p,
  width = 12,
  height = 10,
  dpi = 300
)
