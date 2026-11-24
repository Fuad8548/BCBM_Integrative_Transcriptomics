
library(dplyr)
library(ggplot2)

# =========================================================
# 1. LOAD EXISTING ECM SCORES
# =========================================================

scores <- read.csv(file.choose(), check.names = FALSE)

# =========================================================
# 2. EXTRACT STROMAL CORE ECM SCORE
# =========================================================

core_scores <- scores %>%
  filter(
    cell_type == "Stromal cells",
    signature == "Core_stromal_ECM"
  ) %>%
  select(sample, signature_score, genes_scored)

# Check for duplicated sample records
if (anyDuplicated(core_scores$sample)) {
  stop("Duplicate sample records found for Stromal Core ECM.")
}

# Check expected number of samples
cat("Number of samples:", nrow(core_scores), "\n")
cat("Genes scored:\n")
print(table(core_scores$genes_scored))

# =========================================================
# 3. ORDER SAMPLES BY ECM SCORE
# =========================================================

core_scores <- core_scores %>%
  arrange(signature_score) %>%
  mutate(
    sample = factor(sample, levels = sample)
  )

# =========================================================
# 4. PLOT
# =========================================================

p_5A <- ggplot(
  core_scores,
  aes(x = signature_score, y = sample)
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.6
  ) +
  geom_point(
    size = 3,
    color = "#2878B5",
    alpha = 0.85
  ) +
  labs(
    title = "Stromal Core ECM signature across brain metastases",
    subtitle = "Individual sample scores across 28 brain metastasis specimens",
    x = "Core stromal ECM signature score",
    y = "Brain metastasis sample"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 11)
  )

print(p_5A)

# =========================================================
# 5. EXPORT
# =========================================================

ggsave(
  "Figure5A_Stromal_Core_ECM_distribution.png",
  p_5A,
  width = 8,
  height = 9,
  dpi = 600,
  bg = "white"
)

ggsave(
  "Figure5A_Stromal_Core_ECM_distribution.pdf",
  p_5A,
  width = 8,
  height = 9,
  bg = "white"
)
