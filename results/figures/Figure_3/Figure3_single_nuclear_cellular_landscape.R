# ============================================================
# FIGURE 3: FINAL ASSEMBLY
# 3A Mean cellular composition
# 3B Sample-level stacked composition
# 3C Sample-level composition heatmap
# ============================================================

project_dir <- "E:/computational-biology/BCBM_Integrative_Transcriptomics"

figure_dir <- file.path(
  project_dir,
  "results",
  "figures",
  "Figure_3"
)

composition_dir <- file.path(
  project_dir,
  "02_single_nucleus_transcriptomics",
  "results",
  "cellular_composition"
)

existing_plot_dir <- file.path(
  project_dir,
  "02_single_nucleus_transcriptomics",
  "figures",
  "cellular_composition"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

library(ggplot2)
library(patchwork)
library(grid)
library(png)

# ------------------------------------------------------------
# 1. LOAD SUMMARY DATA
# ------------------------------------------------------------

summary_file <- file.path(
  composition_dir,
  "celltype_composition_summary.csv"
)

summary_df <- read.csv(
  summary_file,
  stringsAsFactors = FALSE
)

# Sort from lowest to highest abundance
summary_df$cell_type <- factor(
  summary_df$cell_type,
  levels = summary_df$cell_type[
    order(summary_df$mean_percent)
  ]
)

# ------------------------------------------------------------
# 2. PANEL 3A
# ------------------------------------------------------------

p3A <- ggplot(
  summary_df,
  aes(
    x = mean_percent,
    y = cell_type
  )
) +
  geom_col(
    fill = "#4776A8",
    width = 0.72,
    color = NA
  ) +
  geom_errorbar(
    aes(
      xmin = pmax(mean_percent - SD_percent, 0),
      xmax = mean_percent + SD_percent
    ),
    height = 0.18,
    linewidth = 0.45
  ) +
  labs(
    title = "Mean cellular composition across samples",
    x = "Mean proportion of nuclei (%)",
    y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 12
    ),
    axis.text.y = element_text(
      size = 9,
      color = "black"
    ),
    axis.text.x = element_text(
      size = 9,
      color = "black"
    ),
    axis.title.x = element_text(size = 10),
    plot.margin = margin(8, 12, 8, 12)
  )

# Save standalone panel A
ggsave(
  file.path(
    figure_dir,
    "Figure3A_mean_cellular_abundance.png"
  ),
  p3A,
  width = 10,
  height = 6.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    figure_dir,
    "Figure3A_mean_cellular_abundance.pdf"
  ),
  p3A,
  width = 10,
  height = 6.5
)

# ------------------------------------------------------------
# 3. LOAD EXISTING PANELS B AND C
# ------------------------------------------------------------

plot3B_file <- file.path(
  existing_plot_dir,
  "sample_cellular_composition_stacked_bar.png"
)

plot3C_file <- file.path(
  existing_plot_dir,
  "sample_cellular_composition_heatmap.png"
)

stopifnot(file.exists(plot3B_file))
stopifnot(file.exists(plot3C_file))

img3B <- png::readPNG(plot3B_file)
img3C <- png::readPNG(plot3C_file)

g3B <- grid::rasterGrob(
  img3B,
  interpolate = TRUE
)

g3C <- grid::rasterGrob(
  img3C,
  interpolate = TRUE
)

# Convert images to patchwork-compatible panels
panel3B <- patchwork::wrap_elements(
  full = g3B
)

panel3C <- patchwork::wrap_elements(
  full = g3C
)

# ------------------------------------------------------------
# 4. ADD PANEL LABELS
# ------------------------------------------------------------

panel3A <- p3A +
  labs(tag = "A") +
  theme(
    plot.tag = element_text(
      size = 16,
      face = "bold"
    ),
    plot.tag.position = "topleft"
  )

panel3B <- panel3B +
  labs(tag = "B") +
  theme(
    plot.tag = element_text(
      size = 16,
      face = "bold"
    ),
    plot.tag.position = "topleft"
  )

panel3C <- panel3C +
  labs(tag = "C") +
  theme(
    plot.tag = element_text(
      size = 16,
      face = "bold"
    ),
    plot.tag.position = "topleft"
  )

# ------------------------------------------------------------
# 5. ASSEMBLE PANELS
# ------------------------------------------------------------

Figure3 <- (panel3A /
  (panel3B | panel3C)) +
  plot_layout(
    heights = c(1.15, 1.35)
  ) +
  plot_annotation(
    title = "Single-nucleus cellular landscape of breast cancer brain metastases",
    theme = theme(
      plot.title = element_text(
        size = 15,
        face = "bold",
        hjust = 0.5
      )
    )
  )

# ------------------------------------------------------------
# 6. EXPORT FINAL FIGURE
# ------------------------------------------------------------

ggsave(
  file.path(
    figure_dir,
    "Figure3_Single_nucleus_cellular_landscape.png"
  ),
  Figure3,
  width = 16,
  height = 13,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    figure_dir,
    "Figure3_Single_nucleus_cellular_landscape.pdf"
  ),
  Figure3,
  width = 16,
  height = 13,
  device = cairo_pdf,
  bg = "white"
)

# ------------------------------------------------------------
# 7. CONFIRM OUTPUTS
# ------------------------------------------------------------

cat("\n========== FIGURE 3 OUTPUTS ==========\n")

print(
  list.files(
    figure_dir,
    pattern = "Figure3",
    full.names = TRUE
  )
)

cat("\nFigure 3 panels:\n")
cat("A: Mean cellular abundance across samples\n")
cat("B: Sample-level stacked composition\n")
cat("C: Sample-level composition heatmap\n")

cat("\nAssembly completed.\n")
