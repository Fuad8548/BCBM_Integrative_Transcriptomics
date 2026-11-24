# ============================================================
# BCBM snRNA-seq
# Milestone 4: ECM Expression vs Cellular Composition
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

output_dir <- "results/ECM_composition_association"
figure_dir <- "figures/ECM_composition_association"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. Load object and metadata
# ------------------------------------------------------------

bcbm <- readRDS(input_file)

meta <- bcbm@meta.data %>%
  select(
    sampleDemuxIDstrict,
    common_cell_clusters
  ) %>%
  filter(
    !is.na(sampleDemuxIDstrict),
    !is.na(common_cell_clusters)
  )

meta <- meta %>%
  rename(
    sample = sampleDemuxIDstrict,
    cell_type = common_cell_clusters
  )

# ------------------------------------------------------------
# 3. ECM genes
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

# ------------------------------------------------------------
# 4. Cell populations to evaluate
# ------------------------------------------------------------

target_celltypes <- c(
  "Stromal cells",
  "Pericytes",
  "Endothelial cells",
  "Astrocytes",
  "Tumor cells",
  "Tumor cells (prol.)"
)

# ------------------------------------------------------------
# 5. Calculate sample-level cell proportions
# ------------------------------------------------------------

cell_counts <- table(meta$sample, meta$cell_type)

cell_proportions <- prop.table(
  cell_counts,
  margin = 1
)

composition_df <- as.data.frame(cell_proportions) %>%
  rename(
    sample = Var1,
    cell_type = Var2,
    proportion = Freq
  )

composition_df <- composition_df %>%
  filter(cell_type %in% target_celltypes)

write.csv(
  composition_df,
  file.path(output_dir, "target_celltype_proportions.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Normalize RNA expression
# ------------------------------------------------------------

DefaultAssay(bcbm) <- "RNA"

bcbm <- NormalizeData(
  bcbm,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# ------------------------------------------------------------
# 7. Extract ECM expression
# ------------------------------------------------------------

expression_df <- FetchData(
  bcbm,
  vars = c(
    "sampleDemuxIDstrict",
    "common_cell_clusters",
    available_genes
  ),
  assay = "RNA",
  layer = "data"
) %>%
  rename(
    sample = sampleDemuxIDstrict,
    cell_type = common_cell_clusters
  ) %>%
  pivot_longer(
    cols = all_of(available_genes),
    names_to = "gene",
    values_to = "expression"
  )

# ------------------------------------------------------------
# 8. Global mean expression per sample
# ------------------------------------------------------------

global_expression <- expression_df %>%
  group_by(sample, gene) %>%
  summarise(
    mean_expression = mean(expression),
    percent_detected = mean(expression > 0) * 100,
    n_nuclei = n(),
    .groups = "drop"
  )

write.csv(
  global_expression,
  file.path(output_dir, "global_ECM_expression_per_sample.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. Mean expression per sample within each cell population
# ------------------------------------------------------------

celltype_expression <- expression_df %>%
  group_by(sample, cell_type, gene) %>%
  summarise(
    mean_expression = mean(expression),
    percent_detected = mean(expression > 0) * 100,
    n_nuclei = n(),
    .groups = "drop"
  )

write.csv(
  celltype_expression,
  file.path(output_dir, "celltype_ECM_expression_per_sample.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. Merge global expression with cellular proportions
# ------------------------------------------------------------

global_association_data <- global_expression %>%
  left_join(
    composition_df %>%
      pivot_wider(
        names_from = cell_type,
        values_from = proportion
      ),
    by = "sample"
  )

write.csv(
  global_association_data,
  file.path(output_dir, "global_expression_composition_merged.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. Spearman correlation function
# ------------------------------------------------------------

run_spearman <- function(data, expression_col, proportion_col) {
  x <- data[[expression_col]]
  y <- data[[proportion_col]]

  valid <- complete.cases(x, y)

  x <- x[valid]
  y <- y[valid]

  if (
    length(x) < 5 ||
      length(unique(x)) < 2 ||
      length(unique(y)) < 2
  ) {
    return(
      data.frame(
        n = length(x),
        rho = NA_real_,
        p_value = NA_real_
      )
    )
  }

  test <- suppressWarnings(
    cor.test(
      x,
      y,
      method = "spearman",
      exact = FALSE
    )
  )

  data.frame(
    n = length(x),
    rho = unname(test$estimate),
    p_value = test$p.value
  )
}

# ------------------------------------------------------------
# 12. Correlation: global gene expression vs cell proportions
# ------------------------------------------------------------

proportion_celltypes <- target_celltypes

global_cor_results <- list()

for (gene in available_genes) {
  gene_data <- global_association_data %>%
    filter(gene == !!gene)

  for (celltype in proportion_celltypes) {
    result <- run_spearman(
      gene_data,
      "mean_expression",
      celltype
    )

    result$gene <- gene
    result$cell_type_proportion <- celltype

    global_cor_results[[length(global_cor_results) + 1]] <-
      result
  }
}

global_cor_results <- bind_rows(global_cor_results) %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(FDR)

write.csv(
  global_cor_results,
  file.path(output_dir, "global_expression_vs_composition_correlations.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 13. Correlation: within-cell-type expression vs own abundance
# ------------------------------------------------------------

within_cor_results <- list()

for (celltype in target_celltypes) {
  for (gene in available_genes) {
    expression_subset <- celltype_expression %>%
      filter(
        cell_type == celltype,
        gene == !!gene
      )

    proportion_subset <- composition_df %>%
      filter(cell_type == celltype) %>%
      select(sample, proportion)

    merged <- expression_subset %>%
      inner_join(
        proportion_subset,
        by = "sample"
      )

    result <- run_spearman(
      merged,
      "mean_expression",
      "proportion"
    )

    result$gene <- gene
    result$cell_type <- celltype

    within_cor_results[[length(within_cor_results) + 1]] <-
      result
  }
}

within_cor_results <- bind_rows(within_cor_results) %>%
  group_by(cell_type) %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH")
  ) %>%
  ungroup() %>%
  arrange(FDR)

write.csv(
  within_cor_results,
  file.path(output_dir, "within_celltype_expression_vs_abundance.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 14. Global correlation heatmap
# ------------------------------------------------------------

p_cor_heatmap <- ggplot(
  global_cor_results,
  aes(
    x = cell_type_proportion,
    y = gene,
    fill = rho
  )
) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    midpoint = 0,
    low = "blue",
    mid = "white",
    high = "red",
    na.value = "grey85"
  ) +
  labs(
    title = "ECM Expression vs Cellular Composition",
    x = "Cell population proportion",
    y = "ECM gene",
    fill = "Spearman rho"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank()
  )

ggsave(
  file.path(figure_dir, "ECM_global_expression_composition_correlation.png"),
  p_cor_heatmap,
  width = 11,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# 15. Within-cell-type correlation heatmap
# ------------------------------------------------------------

p_within_heatmap <- ggplot(
  within_cor_results,
  aes(
    x = gene,
    y = cell_type,
    fill = rho
  )
) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    midpoint = 0,
    low = "blue",
    mid = "white",
    high = "red",
    na.value = "grey85"
  ) +
  labs(
    title = "Cell-Type-Specific ECM Expression vs Cell Abundance",
    x = "ECM gene",
    y = "Cell population",
    fill = "Spearman rho"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank()
  )

ggsave(
  file.path(figure_dir, "ECM_within_celltype_correlation_heatmap.png"),
  p_within_heatmap,
  width = 12,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------
# 16. Export significant global associations
# ------------------------------------------------------------

significant_global <- global_cor_results %>%
  filter(!is.na(FDR), FDR < 0.05)

write.csv(
  significant_global,
  file.path(output_dir, "significant_global_correlations.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 17. Export significant within-cell-type associations
# ------------------------------------------------------------

significant_within <- within_cor_results %>%
  filter(!is.na(FDR), FDR < 0.05)

write.csv(
  significant_within,
  file.path(output_dir, "significant_within_celltype_correlations.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 18. Summary
# ------------------------------------------------------------

cat("\n============================================\n")
cat("ECM-COMPOSITION ASSOCIATION ANALYSIS COMPLETE\n")
cat("============================================\n")

cat("Samples:", length(unique(meta$sample)), "\n")
cat("Genes:", length(available_genes), "\n")
cat("Cell populations tested:", length(target_celltypes), "\n")

cat("\nSignificant global associations:", nrow(significant_global), "\n")
cat(
  "Significant within-cell-type associations:",
  nrow(significant_within),
  "\n"
)

print(significant_global)
print(significant_within)
