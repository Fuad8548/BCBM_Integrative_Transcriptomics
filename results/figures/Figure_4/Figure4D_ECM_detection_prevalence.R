
library(dplyr)
library(ggplot2)

# =========================================================
# 1. LOAD EXISTING ECM DETECTION RESULTS
# =========================================================

det <- read.csv(file.choose(), check.names = FALSE)

required_cols <- c(
  "cell_type", "gene", "n_nuclei",
  "n_detected", "pct_detected"
)

if (!all(required_cols %in% names(det))) {
  stop(
    "Missing required columns. Found: ",
    paste(names(det), collapse = ", ")
  )
}

# =========================================================
# 2. DEFINE ECM GENES AND ACTUAL CELL POPULATIONS
# =========================================================

ecm_genes <- c(
  "AEBP1", "COL14A1", "COL15A1", "COL6A3",
  "CTSK", "CXCL12", "DPT", "LUM", "MMP2",
  "POSTN", "SULF1", "THBS2", "THBS4"
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
# 3. FILTER DATA
# =========================================================

det_plot <- det %>%
  filter(
    gene %in% ecm_genes,
    cell_type %in% celltype_order
  )

# =========================================================
# 4. VALIDATE COVERAGE
# =========================================================

cat("Genes represented:", length(unique(det_plot$gene)), "\n")
cat("Cell populations represented:",
    length(unique(det_plot$cell_type)), "\n")

missing_combinations <- expand.grid(
  gene = ecm_genes,
  cell_type = celltype_order
) %>%
  anti_join(
    det_plot %>% distinct(gene, cell_type),
    by = c("gene", "cell_type")
  )

if (nrow(missing_combinations) > 0) {
  print(missing_combinations)
  warning("Some gene-cell type combinations are missing.")
}

# =========================================================
# 5. ORDER GENES AND CELL POPULATIONS
# =========================================================

det_plot$gene <- factor(
  det_plot$gene,
  levels = rev(ecm_genes)
)

det_plot$cell_type <- factor(
  det_plot$cell_type,
  levels = celltype_order
)

# =========================================================
# 6. PLOT DETECTION PREVALENCE
# =========================================================

p_4D <- ggplot(
  det_plot,
  aes(x = cell_type, y = gene, fill = pct_detected)
) +
  geom_tile(
    color = "white",
    linewidth = 0.4
  ) +
  geom_text(
    aes(label = sprintf("%.1f", pct_detected)),
    size = 3.0,
    color = "black"
  ) +
  scale_fill_gradient(
    low = "#F7F7F7",
    high = "#2166AC",
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    name = "Nuclei with\ngene detected (%)"
  ) +
  labs(
    title = "Prevalence of ECM gene expression across cell populations",
    subtitle = "Percentage of nuclei with detectable expression",
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

print(p_4D)

# =========================================================
# 7. EXPORT
# =========================================================

ggsave(
  "Figure4D_ECM_detection_prevalence.png",
  p_4D,
  width = 10,
  height = 7,
  dpi = 600,
  bg = "white"
)

ggsave(
  "Figure4D_ECM_detection_prevalence.pdf",
  p_4D,
  width = 10,
  height = 7,
  bg = "white"
)
