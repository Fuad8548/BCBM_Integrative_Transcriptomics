
library(dplyr)
library(ggplot2)

# =========================================================
# 1. LOAD EXISTING RESULTS
# =========================================================

scores <- read.csv(file.choose(), check.names = FALSE)
composition <- read.csv(file.choose(), check.names = FALSE)

# =========================================================
# 2. EXTRACT STROMAL CORE ECM SCORE
# =========================================================

core_scores <- scores %>%
  filter(
    cell_type == "Stromal cells",
    signature == "Core_stromal_ECM"
  ) %>%
  select(sample, signature_score)

# =========================================================
# 3. PREPARE CELL COMPOSITION DATA
# =========================================================

# The first column in the composition file contains sample IDs
names(composition)[1] <- "sample"

composition_stromal <- composition %>%
  select(sample, `Stromal cells`) %>%
  rename(stromal_proportion = `Stromal cells`)

# =========================================================
# 4. MERGE BY SAMPLE ID
# =========================================================

plot_data <- inner_join(
  core_scores,
  composition_stromal,
  by = "sample"
) %>%
  filter(
    !is.na(signature_score),
    !is.na(stromal_proportion)
  )

# Convert fraction to percentage for plotting
plot_data$stromal_percent <- plot_data$stromal_proportion * 100

cat("Matched samples:", nrow(plot_data), "\n")

# =========================================================
# 5. SPEARMAN CORRELATION
# =========================================================

cor_result <- cor.test(
  plot_data$signature_score,
  plot_data$stromal_percent,
  method = "spearman",
  exact = FALSE
)

rho <- unname(cor_result$estimate)
p_value <- cor_result$p.value

annotation <- paste0(
  "Spearman \u03c1 = ", round(rho, 3),
  "\nP = ", format.pval(p_value, digits = 3, eps = 0.001),
  "\nN = ", nrow(plot_data)
)

# =========================================================
# 6. FIGURE 5B
# =========================================================

p_5B <- ggplot(
  plot_data,
  aes(x = stromal_percent, y = signature_score)
) +
  geom_point(
    color = "#2878B5",
    size = 3,
    alpha = 0.8
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = "#C43C39",
    fill = "#E9A9A5",
    linewidth = 0.9
  ) +
  annotate(
    "text",
    x = Inf,
    y = Inf,
    label = annotation,
    hjust = 1.15,
    vjust = 1.2,
    size = 4.2
  ) +
  labs(
    title = "Stromal ECM signature is associated with stromal abundance",
    subtitle = "Core stromal ECM score across individual brain metastases",
    x = "Stromal cells (% of nuclei)",
    y = "Core stromal ECM signature score"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    axis.title = element_text(size = 11)
  )

print(p_5B)

# =========================================================
# 7. EXPORT
# =========================================================

ggsave(
  "Figure5B_ECM_vs_stromal_abundance.png",
  p_5B,
  width = 8,
  height = 6,
  dpi = 600,
  bg = "white"
)

ggsave(
  "Figure5B_ECM_vs_stromal_abundance.pdf",
  p_5B,
  width = 8,
  height = 6,
  bg = "white"
)
