# ============================================================
# BCBM snRNA-seq
# Milestone 3: Cellular Origin of Robust ECM Genes
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

output_dir <- "results/ECM_cellular_origin"
figure_dir <- "figures/ECM_cellular_origin"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Load object
# ------------------------------------------------------------

bcbm <- readRDS(input_file)

DefaultAssay(bcbm) <- "RNA"

# ------------------------------------------------------------
# 3. Robust ECM-associated bulk candidates
# ------------------------------------------------------------

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

available_genes <- intersect(
  ecm_genes,
  rownames(bcbm)
)

missing_genes <- setdiff(
  ecm_genes,
  available_genes
)

cat("Genes available:", length(available_genes), "\n")
cat("Missing genes:", missing_genes, "\n")

if (length(available_genes) == 0) {
  stop("None of the ECM genes were found in the Seurat object.")
}

# ------------------------------------------------------------
# 4. Normalize expression
# ------------------------------------------------------------

bcbm <- NormalizeData(
  bcbm,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# ------------------------------------------------------------
# 5. Cellular population expression overview
# ------------------------------------------------------------

p_dot <- DotPlot(
  bcbm,
  features = available_genes,
  group.by = "common_cell_clusters",
  assay = "RNA"
) +
  RotatedAxis() +
  labs(
    title = "Expression of Robust ECM Genes Across BCBM Cell Populations",
    x = "ECM-associated gene",
    y = "Cell population",
    size = "Detection (%)",
    color = "Average expression"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(figure_dir, "ECM_genes_celltype_dotplot.png"),
  p_dot,
  width = 13,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# 6. Extract normalized expression
# ------------------------------------------------------------

expression_data <- FetchData(
  bcbm,
  vars = c(
    "common_cell_clusters",
    "sampleDemuxIDstrict",
    available_genes
  ),
  assay = "RNA",
  layer = "data"
)

expression_data <- expression_data %>%
  rename(
    cell_type = common_cell_clusters,
    sample = sampleDemuxIDstrict
  )

# ------------------------------------------------------------
# 7. Long-format expression table
# ------------------------------------------------------------

expression_long <- expression_data %>%
  pivot_longer(
    cols = all_of(available_genes),
    names_to = "gene",
    values_to = "expression"
  )

write.csv(
  expression_long,
  file.path(output_dir, "ECM_expression_per_nucleus.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Expression by cell population
# ------------------------------------------------------------

celltype_expression <- expression_long %>%
  group_by(cell_type, gene) %>%
  summarise(
    n_nuclei = n(),
    mean_expression = mean(expression),
    median_expression = median(expression),
    percent_detected = mean(expression > 0) * 100,
    .groups = "drop"
  )

write.csv(
  celltype_expression,
  file.path(output_dir, "ECM_expression_by_celltype.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. Expression by sample and cell population
# ------------------------------------------------------------

sample_celltype_expression <- expression_long %>%
  group_by(sample, cell_type, gene) %>%
  summarise(
    n_nuclei = n(),
    mean_expression = mean(expression),
    percent_detected = mean(expression > 0) * 100,
    .groups = "drop"
  )

write.csv(
  sample_celltype_expression,
  file.path(output_dir, "ECM_expression_sample_celltype.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. Cell population heatmap: mean expression
# ------------------------------------------------------------

p_expression_heatmap <- ggplot(
  celltype_expression,
  aes(x = gene, y = cell_type, fill = mean_expression)
) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient(
    low = "white",
    high = "firebrick"
  ) +
  labs(
    title = "Mean ECM Gene Expression Across Cell Populations",
    x = "Gene",
    y = "Cell population",
    fill = "Mean expression"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank()
  )

ggsave(
  file.path(figure_dir, "ECM_celltype_expression_heatmap.png"),
  p_expression_heatmap,
  width = 12,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# 11. Detection frequency heatmap
# ------------------------------------------------------------

p_detection_heatmap <- ggplot(
  celltype_expression,
  aes(x = gene, y = cell_type, fill = percent_detected)
) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient(
    low = "white",
    high = "darkgreen"
  ) +
  labs(
    title = "ECM Gene Detection Across Cell Populations",
    x = "Gene",
    y = "Cell population",
    fill = "Detected (%)"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank()
  )

ggsave(
  file.path(figure_dir, "ECM_celltype_detection_heatmap.png"),
  p_detection_heatmap,
  width = 12,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# 12. Rank cell populations for each ECM gene
# ------------------------------------------------------------

celltype_ranking <- celltype_expression %>%
  group_by(gene) %>%
  arrange(desc(mean_expression), .by_group = TRUE) %>%
  mutate(expression_rank = row_number()) %>%
  ungroup()

write.csv(
  celltype_ranking,
  file.path(output_dir, "ECM_celltype_expression_ranking.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 13. Summary
# ------------------------------------------------------------

cat("\n========================================\n")
cat("ECM CELLULAR ORIGIN ANALYSIS COMPLETE\n")
cat("========================================\n")

cat("Genes analyzed:", length(available_genes), "\n")
cat("Cell populations:", length(unique(expression_data$cell_type)), "\n")
cat("Samples:", length(unique(expression_data$sample)), "\n")

print(
  celltype_expression %>%
    arrange(gene, desc(mean_expression))
)
