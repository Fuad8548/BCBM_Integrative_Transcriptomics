library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# =========================================================
# 1. LOAD ECM PSEUDOBULK DATA
# =========================================================

ecm <- read.csv(file.choose(), check.names = FALSE)

# =========================================================
# 2. DEFINE GENES AND ACTUAL CELL TYPES
# =========================================================

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

celltype_order <- c(
  "Tumor cells",
  "Tumor cells (prol.)",
  "Stromal cells",
  "Pericytes",
  "Endothelial cells",
  "Astrocytes"
)

# =========================================================
# 3. FILTER AND SUMMARIZE
# =========================================================

ecm_mean <- ecm %>%
  filter(
    gene %in% ecm_genes,
    cell_type %in% celltype_order
  ) %>%
  group_by(gene, cell_type) %>%
  summarise(
    mean_logCPM = mean(logCPM, na.rm = TRUE),
    n_samples = sum(!is.na(logCPM)),
    .groups = "drop"
  )

# =========================================================
# 4. COMPLETE ONLY THE SIX REAL CELL TYPES
# =========================================================

ecm_mean_complete <- ecm_mean %>%
  complete(
    gene = ecm_genes,
    cell_type = celltype_order
  )

# =========================================================
# 5. GENE-WISE Z-SCORE ACROSS CELL TYPES
# =========================================================

ecm_heatmap_data <- ecm_mean_complete %>%
  group_by(gene) %>%
  mutate(
    z_score = if (
      sum(!is.na(mean_logCPM)) > 1 &&
        sd(mean_logCPM, na.rm = TRUE) > 0
    ) {
      as.numeric(scale(mean_logCPM))
    } else {
      0
    }
  ) %>%
  ungroup()

# =========================================================
# 6. ORDER ROWS AND COLUMNS
# =========================================================

ecm_heatmap_data$gene <- factor(
  ecm_heatmap_data$gene,
  levels = rev(ecm_genes)
)

ecm_heatmap_data$cell_type <- factor(
  ecm_heatmap_data$cell_type,
  levels = celltype_order
)

# =========================================================
# 7. PLOT
# =========================================================

p_4C_corrected <- ggplot(
  ecm_heatmap_data,
  aes(x = cell_type, y = gene, fill = z_score)
) +
  geom_tile(
    color = "white",
    linewidth = 0.4
  ) +
  scale_fill_gradient2(
    low = "#2878B5",
    mid = "white",
    high = "#C43C39",
    midpoint = 0,
    limits = c(-2, 2),
    oob = scales::squish,
    na.value = "grey85",
    name = "Gene-wise\nZ-score"
  ) +
  labs(
    title = "Cellular distribution of ECM consensus genes",
    subtitle = "Mean pseudobulk expression across available samples within each cell population",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 10
    ),
    axis.text.y = element_text(
      face = "italic",
      size = 10
    ),
    panel.grid = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(10, 15, 10, 10)
  )

print(p_4C_corrected)

# =========================================================
# 8. EXPORT
# =========================================================

ggsave(
  "Figure4C_ECM_celltype_heatmap_corrected.png",
  p_4C_corrected,
  width = 9,
  height = 7,
  dpi = 600,
  bg = "white"
)

ggsave(
  "Figure4C_ECM_celltype_heatmap_corrected.pdf",
  p_4C_corrected,
  width = 9,
  height = 7,
  bg = "white"
)
