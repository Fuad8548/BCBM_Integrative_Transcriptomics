
library(dplyr)
library(tidyr)
library(ggplot2)

# =========================================================
# 1. LOAD EXISTING DATA
# =========================================================

scores <- read.csv(file.choose(), check.names = FALSE)
composition <- read.csv(file.choose(), check.names = FALSE)

names(composition)[1] <- "sample"

# =========================================================
# 2. EXTRACT STROMAL CORE ECM SCORES
# =========================================================

core_scores <- scores %>%
  filter(
    cell_type == "Stromal cells",
    signature == "Core_stromal_ECM"
  ) %>%
  select(sample, signature_score)

# =========================================================
# 3. DEFINE ECM-HIGH AND ECM-LOW GROUPS
# =========================================================

score_median <- median(core_scores$signature_score, na.rm = TRUE)

core_scores <- core_scores %>%
  mutate(
    ECM_group = ifelse(
      signature_score >= score_median,
      "ECM-high",
      "ECM-low"
    )
  )

core_scores$ECM_group <- factor(
  core_scores$ECM_group,
  levels = c("ECM-low", "ECM-high")
)

print(table(core_scores$ECM_group))

# =========================================================
# 4. SELECT MAJOR CELL POPULATIONS
# =========================================================

cell_types_to_plot <- c(
  "Stromal cells",
  "Pericytes",
  "Endothelial cells",
  "Tumor cells",
  "Tumor cells (prol.)",
  "Myeloid cells"
)

composition_long <- composition %>%
  select(sample, all_of(cell_types_to_plot)) %>%
  pivot_longer(
    cols = all_of(cell_types_to_plot),
    names_to = "cell_type",
    values_to = "proportion"
  ) %>%
  mutate(
    proportion_percent = proportion * 100
  )

# =========================================================
# 5. MERGE ECM GROUPS AND CELL COMPOSITION
# =========================================================

plot_data <- inner_join(
  core_scores %>% select(sample, ECM_group, signature_score),
  composition_long,
  by = "sample"
)

# =========================================================
# 6. STATISTICAL COMPARISON
# =========================================================

composition_stats <- plot_data %>%
  group_by(cell_type) %>%
  summarise(
    median_ECM_low = median(
      proportion_percent[ECM_group == "ECM-low"],
      na.rm = TRUE
    ),
    median_ECM_high = median(
      proportion_percent[ECM_group == "ECM-high"],
      na.rm = TRUE
    ),
    difference = median_ECM_high - median_ECM_low,
    p_value = wilcox.test(
      proportion_percent ~ ECM_group
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(FDR)

print(composition_stats)

# =========================================================
# 7. FIGURE 5C
# =========================================================

p_5C <- ggplot(
  plot_data,
  aes(x = ECM_group, y = proportion_percent, fill = ECM_group)
) +
  geom_boxplot(
    width = 0.55,
    alpha = 0.65,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 1.8,
    alpha = 0.7
  ) +
  facet_wrap(
    ~cell_type,
    scales = "free_y",
    ncol = 3
  ) +
  scale_fill_manual(
    values = c(
      "ECM-low" = "#2878B5",
      "ECM-high" = "#C43C39"
    )
  ) +
  labs(
    title = "Cellular composition of ECM-high and ECM-low metastases",
    subtitle = "Exploratory median split of stromal Core ECM scores",
    x = NULL,
    y = "Cell population (% of nuclei)",
    fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(size = 9),
    legend.position = "bottom"
  )

print(p_5C)

# =========================================================
# 8. EXPORT
# =========================================================

ggsave(
  "Figure5C_ECM_group_cell_composition.png",
  p_5C,
  width = 11,
  height = 8,
  dpi = 600,
  bg = "white"
)

ggsave(
  "Figure5C_ECM_group_cell_composition.pdf",
  p_5C,
  width = 11,
  height = 8,
  bg = "white"
)

write.csv(
  composition_stats,
  "Figure5C_ECM_group_composition_statistics.csv",
  row.names = FALSE
)
