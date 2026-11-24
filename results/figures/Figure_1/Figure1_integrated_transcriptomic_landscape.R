library(ggplot2)
library(patchwork)
library(png)
library(grid)

# =========================================================
# PANEL A: STUDY DESIGN WORKFLOW
# =========================================================

workflow_nodes <- data.frame(
  x = c(1.5, 4.1, 4.1, 7.0, 9.6, 12.2),
  y = c(2.0, 2.8, 1.2, 2.0, 2.0, 2.0),
  label = c(
    "Discovery cohort\nGSE184869\n45 paired samples",
    "Model A\nPaired limma\nPatient-adjusted",
    "Model B\nPurity-adjusted limma\nPatient + purity",
    "Independent validation\nGSE125989\n16 paired samples",
    "Three-way consensus\nDirection consistency\nacross analyses",
    "Biological interpretation\nECM remodelling\nMetabolic adaptation"
  ),
  fill = c(
    "#DCEAF5",
    "#DCEAF5",
    "#DDF0E8",
    "#E8E8F2",
    "#FBE8D6",
    "#D9EEEA"
  ),
  border = c(
    "#2878B5",
    "#2878B5",
    "#31916F",
    "#65758A",
    "#D68132",
    "#16877C"
  ),
  stringsAsFactors = FALSE
)

workflow_edges <- data.frame(
  x = c(2.5, 2.5, 5.2, 5.2, 8.1, 10.7),
  y = c(2.0, 2.0, 2.8, 1.2, 2.0, 2.0),
  xend = c(3.1, 3.1, 6.0, 6.0, 8.5, 11.1),
  yend = c(2.8, 1.2, 2.0, 2.0, 2.0, 2.0)
)

p_A <- ggplot() +
  geom_segment(
    data = workflow_edges,
    aes(x = x, y = y, xend = xend, yend = yend),
    arrow = arrow(length = unit(0.16, "cm"), type = "closed"),
    linewidth = 0.75,
    color = "#344E68"
  ) +
  geom_label(
    data = workflow_nodes,
    aes(x = x, y = y, label = label, fill = fill, color = border),
    text.colour = "#172B3A",
    linewidth = 0.75,
    label.r = unit(0.15, "cm"),
    label.padding = unit(0.30, "cm"),
    size = 3.4,
    fontface = "bold",
    lineheight = 1.15,
    show.legend = FALSE
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_x_continuous(limits = c(0.2, 13.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.45, 3.35), expand = c(0, 0)) +
  labs(title = "A  Study design and integrative analysis") +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0),
    plot.margin = margin(6, 10, 6, 10)
  )

# Export standalone Panel A
ggsave(
  "Figure1A_study_design.png",
  p_A,
  width = 16,
  height = 4.5,
  dpi = 600,
  bg = "white"
)

# =========================================================
# LOAD & TAG PANELS B, C, D
# =========================================================

panel_files <- c(
  "Figure1B_discovery_volcano.png",
  "Figure1C_purity_adjusted_volcano.png",
  "Figure1D_three_way_consensus_effects.png"
)

missing_files <- panel_files[!file.exists(panel_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing panel files:\n",
    paste(missing_files, collapse = "\n")
  )
}

make_tagged_panel <- function(file, tag_label) {
  img <- png::readPNG(file)
  grob <- grid::rasterGrob(img, interpolate = FALSE)

  ggplot() +
    annotation_custom(grob, xmin = 0, xmax = 1, ymin = 0, ymax = 1) +
    annotate(
      "text",
      x = 0.01,
      y = 0.99,
      label = tag_label,
      hjust = 0,
      vjust = 1,
      fontface = "bold",
      size = 5,
      color = "#111111"
    ) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(2, 4, 2, 4))
}

p_B <- make_tagged_panel(panel_files[1], "B")
p_C <- make_tagged_panel(panel_files[2], "C")
p_D <- make_tagged_panel(panel_files[3], "D")

# =========================================================
# ASSEMBLE FIGURE 1
# =========================================================

bottom_row <- (p_B | p_C | p_D) +
  plot_layout(widths = c(1, 1, 1))

Figure1 <- (p_A / bottom_row) +
  plot_layout(heights = c(1.0, 2.15)) +
  plot_annotation(
    theme = theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(6, 8, 6, 8)
    )
  )

print(Figure1)

# =========================================================
# EXPORT FINAL FIGURE
# =========================================================

ggsave(
  "Figure1_integrated_transcriptomic_landscape.png",
  Figure1,
  width = 18,
  height = 11.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  "Figure1_integrated_transcriptomic_landscape.pdf",
  Figure1,
  width = 18,
  height = 11.5,
  device = cairo_pdf,
  bg = "white"
)
