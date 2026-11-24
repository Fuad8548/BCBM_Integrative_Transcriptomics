# ============================================================
# FIGURE 5: SAMPLE-LEVEL ECM HETEROGENEITY
# ASSEMBLE PANELS 5A-5D
# ============================================================

project_dir <- "E:/computational-biology/BCBM_Integrative_Transcriptomics"

figure5_dir <- file.path(
  project_dir,
  "results",
  "figures",
  "Figure_5"
)

# ------------------------------------------------------------
# 1. Locate existing panels
# ------------------------------------------------------------

panel_files <- c(
  "5A" = file.path(
    figure5_dir,
    "Figure5A_Stromal_Core_ECM_distribution.png"
  ),
  "5B" = file.path(
    figure5_dir,
    "Figure5B_ECM_vs_stromal_abundance.png"
  ),
  "5C" = file.path(
    figure5_dir,
    "Figure5C_ECM_group_cell_composition.png"
  ),
  "5D" = file.path(
    figure5_dir,
    "Figure5D_ECM_gene_comparison.png"
  )
)

# Verify all panels exist
missing_files <- panel_files[!file.exists(panel_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing panel files:\n",
    paste(missing_files, collapse = "\n")
  )
}

cat("All four panels found.\n")

# ------------------------------------------------------------
# 2. Load packages
# ------------------------------------------------------------

library(png)
library(grid)
library(patchwork)

# ------------------------------------------------------------
# 3. Import panels
# ------------------------------------------------------------

panel_grobs <- lapply(
  panel_files,
  function(f) {
    img <- png::readPNG(f)

    grid::rasterGrob(
      img,
      interpolate = TRUE
    )
  }
)

# Convert image grobs into patchwork-compatible objects
# Convert image grobs into patchwork-compatible objects
panel_plots <- lapply(
  panel_grobs,
  patchwork::wrap_elements # or wrap_elements(panel = g)
)

# ------------------------------------------------------------
# 4. Assemble 2 x 2 layout
# ------------------------------------------------------------

Figure5 <- (panel_plots[[1]] | panel_plots[[2]]) /
  (panel_plots[[3]] | panel_plots[[4]])

Figure5 <- Figure5 +
  plot_layout(
    widths = c(1, 1),
    heights = c(1, 1)
  )

# ------------------------------------------------------------
# 5. Export
# ------------------------------------------------------------

ggsave(
  file.path(
    figure5_dir,
    "Figure5_ECM_sample_heterogeneity.png"
  ),
  plot = Figure5,
  width = 16,
  height = 12,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(
    figure5_dir,
    "Figure5_ECM_sample_heterogeneity.pdf"
  ),
  plot = Figure5,
  width = 16,
  height = 12,
  bg = "white"
)

cat("\n========================================\n")
cat("FIGURE 5 ASSEMBLED SUCCESSFULLY\n")
cat("========================================\n")

cat(
  "\nPNG:\n",
  file.path(
    figure5_dir,
    "Figure5_ECM_sample_heterogeneity.png"
  ),
  "\n\nPDF:\n",
  file.path(
    figure5_dir,
    "Figure5_ECM_sample_heterogeneity.pdf"
  ),
  "\n"
)
