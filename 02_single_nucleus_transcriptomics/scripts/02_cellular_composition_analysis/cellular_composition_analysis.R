# ============================================================
# bc_bm Single-Nucleus RNA-seq
# Milestone 2: Cellular Composition Analysis
# Dataset: GSE324453
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

input_file <- "data/raw/GSE324453_seurat_filtered_annotated.rds"

output_dir <- "results/cellular_composition"
figure_dir <- "figures/cellular_composition"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Load original Seurat object
# ------------------------------------------------------------

bc_bm <- readRDS(input_file)

meta <- bc_bm@meta.data

sample_col <- "sampleDemuxIDstrict"
celltype_col <- "common_cell_clusters"

# Check required metadata
stopifnot(sample_col %in% colnames(meta))
stopifnot(celltype_col %in% colnames(meta))

meta <- meta %>%
  select(all_of(c(sample_col, celltype_col))) %>%
  filter(
    !is.na(.data[[sample_col]]),
    !is.na(.data[[celltype_col]])
  )

colnames(meta) <- c("sample", "cell_type")

# ------------------------------------------------------------
# 3. Cell counts per sample
# ------------------------------------------------------------

cell_counts <- table(meta$sample, meta$cell_type)

cell_counts_df <- as.data.frame.matrix(cell_counts)

write.csv(
  cell_counts_df,
  file.path(output_dir, "sample_cell_counts.csv")
)

# ------------------------------------------------------------
# 4. Cell proportions per sample
# ------------------------------------------------------------

cell_proportions <- prop.table(cell_counts, margin = 1)

cell_proportions_df <- as.data.frame.matrix(cell_proportions)

write.csv(
  cell_proportions_df,
  file.path(output_dir, "sample_cell_proportions.csv")
)

# ------------------------------------------------------------
# 5. Long-format composition table
# ------------------------------------------------------------

composition_long <- as.data.frame(cell_counts) %>%
  rename(
    sample = Var1,
    cell_type = Var2,
    cell_count = Freq
  ) %>%
  left_join(
    as.data.frame(cell_proportions) %>%
      rename(
        sample = Var1,
        cell_type = Var2,
        proportion = Freq
      ),
    by = c("sample", "cell_type")
  ) %>%
  mutate(percent = proportion * 100)

write.csv(
  composition_long,
  file.path(output_dir, "composition_long.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Sample-level total nuclei
# ------------------------------------------------------------

sample_sizes <- composition_long %>%
  group_by(sample) %>%
  summarise(
    total_nuclei = sum(cell_count),
    .groups = "drop"
  )

write.csv(
  sample_sizes,
  file.path(output_dir, "sample_total_nuclei.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Cell-type summary across samples
# ------------------------------------------------------------

celltype_summary <- composition_long %>%
  group_by(cell_type) %>%
  summarise(
    total_nuclei = sum(cell_count),
    mean_percent = mean(percent),
    median_percent = median(percent),
    SD_percent = sd(percent),
    min_percent = min(percent),
    max_percent = max(percent),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_percent))

write.csv(
  celltype_summary,
  file.path(output_dir, "celltype_composition_summary.csv"),
  row.names = FALSE
)

print(celltype_summary)

# ------------------------------------------------------------
# 8. Order samples by total tumor-cell proportion
# ------------------------------------------------------------

tumor_types <- c(
  "Tumor cells",
  "Tumor cells (prol.)"
)

tumor_fraction <- composition_long %>%
  filter(cell_type %in% tumor_types) %>%
  group_by(sample) %>%
  summarise(
    tumor_percent = sum(percent),
    .groups = "drop"
  ) %>%
  arrange(desc(tumor_percent))

sample_order <- tumor_fraction$sample

composition_long$sample <- factor(
  composition_long$sample,
  levels = sample_order
)

# ------------------------------------------------------------
# 9. 100% stacked bar plot
# ------------------------------------------------------------

p_composition <- ggplot(
  composition_long,
  aes(x = sample, y = percent, fill = cell_type)
) +
  geom_col(width = 0.85) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0),
    breaks = seq(0, 100, 20)
  ) +
  labs(
    title = "Cellular Composition of bc_bm Samples",
    x = "Sample",
    y = "Cell proportion (%)",
    fill = "Cell population"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    panel.grid.major.x = element_blank()
  )

ggsave(
  file.path(figure_dir, "sample_cellular_composition_stacked_bar.png"),
  p_composition,
  width = 14,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# 10. Cellular composition heatmap
# ------------------------------------------------------------

heatmap_data <- composition_long

heatmap_data$sample <- factor(
  heatmap_data$sample,
  levels = sample_order
)

p_heatmap <- ggplot(
  heatmap_data,
  aes(x = cell_type, y = sample, fill = percent)
) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_gradient(
    low = "white",
    high = "steelblue",
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    title = "Cellular Composition Across bc_bm Samples",
    x = "Cell population",
    y = "Sample",
    fill = "Proportion"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    panel.grid = element_blank()
  )

ggsave(
  file.path(figure_dir, "sample_cellular_composition_heatmap.png"),
  p_heatmap,
  width = 12,
  height = 10,
  dpi = 300
)

# ------------------------------------------------------------
# 11. Save sample tumor fraction ranking
# ------------------------------------------------------------

write.csv(
  tumor_fraction,
  file.path(output_dir, "sample_tumor_fraction.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 12. Final summary
# ------------------------------------------------------------

cat("\n========================================\n")
cat("CELLULAR COMPOSITION ANALYSIS COMPLETE\n")
cat("========================================\n")

cat("Samples:", nrow(cell_counts), "\n")
cat("Cell populations:", ncol(cell_counts), "\n")
cat("Total nuclei:", sum(cell_counts), "\n")

cat("\nOutput directory:", output_dir, "\n")
cat("Figure directory:", figure_dir, "\n")

print(tumor_fraction)
