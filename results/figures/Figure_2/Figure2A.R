# ============================================================
# FIGURE 2: FUNCTIONAL PATHWAY LANDSCAPE
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

project_dir <- "E:/computational-biology/BCBM_Integrative_Transcriptomics"

pathway_dir <- file.path(
  project_dir,
  "01_bulk_transcriptomics",
  "results",
  "tables",
  "consensus_pathways"
)

figure_dir <- file.path(
  project_dir,
  "results",
  "figures",
  "Figure_2"
)

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

go <- read.csv(
  file.path(pathway_dir, "GO_robust_consensus.csv"),
  check.names = FALSE
)

kegg <- read.csv(
  file.path(pathway_dir, "KEGG_robust_consensus.csv"),
  check.names = FALSE
)

# ------------------------------------------------------------
# 1. Define selected pathways
# ------------------------------------------------------------

metabolic_terms <- c(
  "Oxidative phosphorylation",
  "aerobic respiration",
  "respiratory electron transport chain"
)

ecm_terms <- c(
  "extracellular matrix organization",
  "collagen fibril organization",
  "collagen metabolic process",
  "extracellular matrix assembly"
)

immune_terms <- c(
  "Cytokine-cytokine receptor interaction",
  "Integrin signaling",
  "negative regulation of leukocyte migration"
)

# ------------------------------------------------------------
# 2. Combine GO and KEGG
# ------------------------------------------------------------

pathway_data <- bind_rows(
  go %>% mutate(source = "GO"),
  kegg %>% mutate(source = "KEGG")
)

# ------------------------------------------------------------
# 3. Convert NES columns into long format
# ------------------------------------------------------------

nes_long <- pathway_data %>%
  select(
    Description,
    NES_discovery,
    NES_purity,
    NES_validation
  ) %>%
  pivot_longer(
    cols = c(
      NES_discovery,
      NES_purity,
      NES_validation
    ),
    names_to = "Dataset",
    values_to = "NES"
  ) %>%
  mutate(
    Dataset = recode(
      Dataset,
      NES_discovery = "Discovery",
      NES_purity = "Purity-adjusted",
      NES_validation = "Validation"
    ),
    Dataset = factor(
      Dataset,
      levels = c(
        "Discovery",
        "Purity-adjusted",
        "Validation"
      )
    )
  )

# ------------------------------------------------------------
# 4. Plotting function
# ------------------------------------------------------------

plot_pathways <- function(terms, title, subtitle = NULL) {
  plot_data <- nes_long %>%
    filter(Description %in% terms) %>%
    mutate(
      Description = factor(
        Description,
        levels = rev(terms)
      )
    )

  ggplot(
    plot_data,
    aes(
      x = NES,
      y = Description,
      color = Dataset,
      group = Dataset
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey50"
    ) +
    geom_line(
      aes(group = Description),
      color = "grey75",
      linewidth = 0.4
    ) +
    geom_point(
      size = 2.6
    ) +
    scale_color_manual(
      values = c(
        "Discovery" = "#2166AC",
        "Purity-adjusted" = "#D97721",
        "Validation" = "#27934A"
      )
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Normalized enrichment score (NES)",
      y = NULL,
      color = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 12
      ),
      plot.subtitle = element_text(size = 9),
      axis.text.y = element_text(
        size = 8.5,
        color = "black"
      ),
      axis.text.x = element_text(size = 8.5),
      legend.position = "bottom",
      panel.grid.major.x = element_line(
        color = "grey90",
        linewidth = 0.3
      )
    )
}

# ------------------------------------------------------------
# 5. Panel 2B: Metabolic reprogramming
# ------------------------------------------------------------

p2B <- plot_pathways(
  metabolic_terms,
  "2B | Mitochondrial energy metabolism",
  "Consistent enrichment of respiratory pathways"
)

# ------------------------------------------------------------
# 6. Panel 2C: ECM remodeling
# ------------------------------------------------------------

p2C <- plot_pathways(
  ecm_terms,
  "2C | Extracellular matrix and collagen remodeling",
  "Consistent depletion of ECM-associated transcriptional programs"
)

# ------------------------------------------------------------
# 7. Panel 2D: Immune and microenvironment signaling
# ------------------------------------------------------------

p2D <- plot_pathways(
  immune_terms,
  "2D | Microenvironment-associated signaling",
  "Cytokine, integrin, and leukocyte-associated pathways"
)

# ------------------------------------------------------------
# 8. Panel 2A: Global pathway overview
# ------------------------------------------------------------

global_terms <- c(
  metabolic_terms,
  ecm_terms[1:2],
  "angiogenesis",
  immune_terms
)

global_data <- nes_long %>%
  filter(Description %in% global_terms) %>%
  mutate(
    Description = factor(
      Description,
      levels = rev(global_terms)
    )
  )

p2A <- ggplot(
  global_data,
  aes(
    x = Dataset,
    y = Description,
    fill = NES
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.5
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0
  ) +
  labs(
    title = "2A | Global pathway landscape",
    subtitle = "Direction and consistency across datasets",
    x = NULL,
    y = NULL,
    fill = "NES"
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 12
    ),
    axis.text.y = element_text(
      size = 8,
      color = "black"
    ),
    axis.text.x = element_text(
      size = 8.5,
      angle = 30,
      hjust = 1
    ),
    legend.position = "right"
  )

# ------------------------------------------------------------
# 9. Export individual panels
# ------------------------------------------------------------

ggsave(
  file.path(figure_dir, "Figure2A_Global_pathway_landscape.pdf"),
  p2A,
  width = 8,
  height = 5.5
)

ggsave(
  file.path(figure_dir, "Figure2B_Metabolic_reprogramming.pdf"),
  p2B,
  width = 7,
  height = 3.5
)

ggsave(
  file.path(figure_dir, "Figure2C_ECM_remodeling.pdf"),
  p2C,
  width = 8,
  height = 4.5
)

ggsave(
  file.path(figure_dir, "Figure2D_Immune_microenvironment.pdf"),
  p2D,
  width = 7.5,
  height = 3.5
)

# ------------------------------------------------------------
# 10. Assemble combined Figure 2
# ------------------------------------------------------------

Figure2 <- (p2A |
  (p2B / p2C / p2D)) +
  plot_layout(
    widths = c(1.05, 1.5),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom"
  )

ggsave(
  file.path(
    figure_dir,
    "Figure2_Functional_pathway_landscape.pdf"
  ),
  Figure2,
  width = 15,
  height = 11,
  limitsize = FALSE
)

ggsave(
  file.path(
    figure_dir,
    "Figure2_Functional_pathway_landscape.png"
  ),
  Figure2,
  width = 15,
  height = 11,
  dpi = 400,
  limitsize = FALSE
)

# ------------------------------------------------------------
# 11. Confirm output
# ------------------------------------------------------------

cat("\nFigure 2 saved in:\n", figure_dir, "\n")

list.files(
  figure_dir,
  pattern = "Figure2",
  full.names = TRUE
)
