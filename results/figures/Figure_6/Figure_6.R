# ============================================================
# FIGURE 6: INTEGRATIVE BIOLOGICAL MODEL
# BCBM TRANSCRIPTOMIC REMODELING
# ============================================================

library(ggplot2)
library(grid)

project_dir <- "E:/computational-biology/BCBM_Integrative_Transcriptomics"

output_dir <- file.path(
  project_dir,
  "results",
  "figures",
  "Figure_6"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 1. Define main section positions
# ------------------------------------------------------------

sections <- data.frame(
  x = c(1.7, 5.0, 8.3),
  title = c(
    "A. BULK TRANSCRIPTOMIC REMODELING",
    "B. CELLULAR CONTEXT",
    "C. ECM HETEROGENEITY"
  ),
  subtitle = c(
    "Primary tumour vs brain metastasis",
    "Single-nucleus profiling",
    "Across metastatic samples"
  )
)

# ------------------------------------------------------------
# 2. Define biological content
# ------------------------------------------------------------

boxes <- data.frame(
  x = c(
    1.7,
    1.7,
    1.7,
    5.0,
    5.0,
    5.0,
    5.0,
    8.3,
    8.3,
    8.3
  ),

  y = c(
    5.0,
    3.5,
    2.0,
    5.3,
    4.0,
    2.7,
    1.4,
    4.7,
    3.2,
    1.7
  ),

  label = c(
    # A
    "↑ Oxidative phosphorylation\n↑ Mitochondrial respiration",
    "↓ ECM and collagen programs\n↓ Matrix organization",
    "↓ Immune and cytokine signaling",

    # B
    "Tumour compartment\nMalignant + proliferating cells",
    "Stromal ECM hub\nBroad ECM-associated expression",
    "Vascular niches\nEndothelial cells + pericytes",
    "Neural-associated populations\nAstrocytes + other brain-associated cells",

    # C
    "ECM-high samples\nHigher ECM signature scores",
    "Continuous variation\nIntermediate ECM expression",
    "ECM-low samples\nLower ECM signature scores"
  ),

  group = c(
    "A",
    "A",
    "A",
    "B",
    "B",
    "B",
    "B",
    "C",
    "C",
    "C"
  )
)

# ------------------------------------------------------------
# 3. Define conceptual connections
# ------------------------------------------------------------

arrows <- data.frame(
  x = c(2.85, 6.15),
  xend = c(3.55, 6.85),
  y = c(3.5, 3.5),
  yend = c(3.5, 3.5)
)

# ------------------------------------------------------------
# 4. Define section colors
# ------------------------------------------------------------

section_colors <- c(
  A = "#DCEBDD",
  B = "#DDEAF5",
  C = "#EDE3F2"
)

# ------------------------------------------------------------
# 5. Plot
# ------------------------------------------------------------

p6 <- ggplot() +

  # Main section backgrounds
  annotate(
    "rect",
    xmin = 0.35,
    xmax = 3.05,
    ymin = 0.65,
    ymax = 6.35,
    fill = "#F5F8F5",
    color = "#BBC8BC",
    linewidth = 0.7,
    radius = unit(0.12, "in")
  ) +

  annotate(
    "rect",
    xmin = 3.65,
    xmax = 6.35,
    ymin = 0.65,
    ymax = 6.35,
    fill = "#F4F7FA",
    color = "#BCCBD8",
    linewidth = 0.7
  ) +

  annotate(
    "rect",
    xmin = 6.95,
    xmax = 9.65,
    ymin = 0.65,
    ymax = 6.35,
    fill = "#F8F5F9",
    color = "#CBBCCF",
    linewidth = 0.7
  ) +

  # Section titles
  geom_text(
    data = sections,
    aes(x = x, y = 6.05, label = title),
    fontface = "bold",
    size = 3.7
  ) +

  geom_text(
    data = sections,
    aes(x = x, y = 5.72, label = subtitle),
    size = 3.2,
    color = "grey35"
  ) +

  # Biological content boxes
  geom_label(
    data = boxes,
    aes(
      x = x,
      y = y,
      label = label,
      fill = group
    ),
    size = 3.6,
    fontface = "plain",
    lineheight = 1.15,
    label.size = 0.35,
    label.padding = unit(0.32, "lines"),
    label.r = unit(0.12, "lines"),
    color = "#202020",
    show.legend = FALSE
  ) +

  scale_fill_manual(
    values = section_colors
  ) +

  # Conceptual connecting arrows
  geom_segment(
    data = arrows,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend
    ),
    arrow = arrow(
      length = unit(0.18, "in"),
      type = "closed"
    ),
    linewidth = 0.8,
    color = "grey40"
  ) +

  # Footer interpretation
  annotate(
    "text",
    x = 5,
    y = 0.25,
    label = paste(
      "Integrative interpretation of bulk and single-nucleus transcriptomic findings;",
      "associations do not establish causality."
    ),
    size = 3.1,
    color = "grey35"
  ) +

  coord_cartesian(
    xlim = c(0, 10),
    ylim = c(0, 7),
    clip = "off"
  ) +

  theme_void() +

  theme(
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    plot.margin = margin(
      15,
      15,
      15,
      15
    )
  )

# ------------------------------------------------------------
# 6. Export
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "Figure6_Integrative_Biological_Model.png"
  ),
  plot = p6,
  width = 15,
  height = 9,
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(
    output_dir,
    "Figure6_Integrative_Biological_Model.pdf"
  ),
  plot = p6,
  width = 15,
  height = 9,
  bg = "white"
)

# Display in RStudio
print(p6)

cat("\nFigure 6 generated successfully.\n")
