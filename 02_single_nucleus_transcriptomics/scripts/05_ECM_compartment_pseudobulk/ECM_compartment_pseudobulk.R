# ============================================================
# BCBM snRNA-seq | GSE324453
# Milestone 5: RAM-efficient ECM pseudobulk
# Extract only 13 genes and 6 cell populations
# ============================================================

library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

input_file <- "data/raw/GSE324453_seurat_filtered_annotated.rds"

output_dir <- "results/ECM_pseudobulk"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Load object
# ------------------------------------------------------------

bcbm <- readRDS(input_file)

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

target_celltypes <- c(
  "Stromal cells",
  "Pericytes",
  "Endothelial cells",
  "Astrocytes",
  "Tumor cells",
  "Tumor cells (prol.)"
)

# ------------------------------------------------------------
# 3. Metadata first
# ------------------------------------------------------------

meta <- bcbm@meta.data

meta$cell <- rownames(meta)

meta_selected <- meta %>%
  filter(
    !is.na(sampleDemuxIDstrict),
    !is.na(common_cell_clusters),
    common_cell_clusters %in% target_celltypes
  ) %>%
  mutate(
    sample = as.character(sampleDemuxIDstrict),
    cell_type = as.character(common_cell_clusters),
    group = paste(sample, cell_type, sep = "__")
  )

selected_cells <- meta_selected$cell

# ------------------------------------------------------------
# 4. Identify available ECM genes
# ------------------------------------------------------------

available_genes <- intersect(
  ecm_genes,
  rownames(bcbm)
)

missing_genes <- setdiff(
  ecm_genes,
  available_genes
)

cat("ECM genes available:", length(available_genes), "\n")
cat("Selected nuclei:", length(selected_cells), "\n")
cat("Missing genes:", missing_genes, "\n")

# ------------------------------------------------------------
# 5. Extract only target genes and selected nuclei
# ------------------------------------------------------------

ecm_counts <- LayerData(
  bcbm,
  assay = "RNA",
  layer = "counts",
  features = available_genes,
  cells = selected_cells
)

# Align metadata to extracted matrix columns
meta_selected <- meta_selected[
  match(colnames(ecm_counts), meta_selected$cell),
]

stopifnot(
  identical(meta_selected$cell, colnames(ecm_counts))
)

cat(
  "Extracted matrix dimensions:",
  nrow(ecm_counts),
  "genes x",
  ncol(ecm_counts),
  "nuclei\n"
)

# ------------------------------------------------------------
# 6. Aggregate raw counts by sample x cell type
# ------------------------------------------------------------

group_levels <- unique(meta_selected$group)

group_factor <- factor(
  meta_selected$group,
  levels = group_levels
)

aggregation_matrix <- Matrix::sparse.model.matrix(
  ~ 0 + group_factor
)

colnames(aggregation_matrix) <- group_levels

pseudobulk_counts <- ecm_counts %*% aggregation_matrix

# ------------------------------------------------------------
# 7. Pseudobulk metadata
# ------------------------------------------------------------

group_info <- meta_selected %>%
  group_by(group, sample, cell_type) %>%
  summarise(
    n_nuclei = n(),
    .groups = "drop"
  )

group_info <- group_info[
  match(colnames(pseudobulk_counts), group_info$group),
]

stopifnot(
  identical(group_info$group, colnames(pseudobulk_counts))
)

write.csv(
  group_info,
  file.path(output_dir, "pseudobulk_group_metadata.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Export raw pseudobulk counts
# ------------------------------------------------------------

raw_counts_df <- as.data.frame(
  as.matrix(pseudobulk_counts)
)

raw_counts_df$gene <- rownames(raw_counts_df)

write.csv(
  raw_counts_df,
  file.path(output_dir, "ECM_pseudobulk_raw_counts.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. Normalize pseudobulk expression
# ------------------------------------------------------------

# Library size uses all genes in each selected cell population,
# not only the 13 ECM genes.

# Use full layer column sums without loading the full matrix.
# This obtains total RNA counts per selected nucleus.

selected_library_sizes <- bcbm$nCount_RNA[
  match(colnames(ecm_counts), colnames(bcbm))
]

stopifnot(
  length(selected_library_sizes) == ncol(ecm_counts),
  !anyNA(selected_library_sizes)
)

pseudobulk_library_sizes <- as.numeric(
  Matrix::crossprod(
    aggregation_matrix,
    selected_library_sizes
  )
)

names(pseudobulk_library_sizes) <- group_levels

# CPM
cpm_matrix <- sweep(
  as.matrix(pseudobulk_counts),
  2,
  pseudobulk_library_sizes,
  "/"
) *
  1e6

# log2 CPM
logCPM_matrix <- log2(cpm_matrix + 0.5)

# ------------------------------------------------------------
# 10. Export expression matrix
# ------------------------------------------------------------

logCPM_df <- as.data.frame(logCPM_matrix)

logCPM_df$gene <- rownames(logCPM_df)

write.csv(
  logCPM_df,
  file.path(output_dir, "ECM_pseudobulk_logCPM.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 11. Long-format expression
# ------------------------------------------------------------

logCPM_long <- as.data.frame(as.table(logCPM_matrix)) %>%
  rename(
    gene = Var1,
    group = Var2,
    logCPM = Freq
  ) %>%
  left_join(
    group_info,
    by = "group"
  )

write.csv(
  logCPM_long,
  file.path(output_dir, "ECM_pseudobulk_logCPM_long.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 12. Expression summary by cell population
# ------------------------------------------------------------

expression_summary <- logCPM_long %>%
  group_by(cell_type, gene) %>%
  summarise(
    n_samples = n(),
    mean_logCPM = mean(logCPM),
    median_logCPM = median(logCPM),
    SD_logCPM = sd(logCPM),
    min_logCPM = min(logCPM),
    max_logCPM = max(logCPM),
    .groups = "drop"
  ) %>%
  arrange(cell_type, desc(mean_logCPM))

write.csv(
  expression_summary,
  file.path(output_dir, "ECM_pseudobulk_expression_summary.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 13. Export individual compartment matrices
# ------------------------------------------------------------

for (celltype in target_celltypes) {
  celltype_data <- logCPM_long %>%
    filter(cell_type == celltype) %>%
    select(sample, gene, logCPM) %>%
    pivot_wider(
      names_from = sample,
      values_from = logCPM
    )

  safe_name <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    celltype
  )

  write.csv(
    celltype_data,
    file.path(
      output_dir,
      paste0("ECM_logCPM_", safe_name, ".csv")
    ),
    row.names = FALSE
  )
}

# ------------------------------------------------------------
# 14. Cell population sample summary
# ------------------------------------------------------------

celltype_sample_summary <- group_info %>%
  group_by(cell_type) %>%
  summarise(
    n_samples = n(),
    total_nuclei = sum(n_nuclei),
    median_nuclei_per_sample = median(n_nuclei),
    min_nuclei_per_sample = min(n_nuclei),
    max_nuclei_per_sample = max(n_nuclei),
    .groups = "drop"
  )

write.csv(
  celltype_sample_summary,
  file.path(
    output_dir,
    "pseudobulk_celltype_sample_summary.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 15. Final summary
# ------------------------------------------------------------

cat("\n========================================\n")
cat("RAM-EFFICIENT PSEUDOBULK COMPLETE\n")
cat("========================================\n")

cat("Samples:", length(unique(meta_selected$sample)), "\n")
cat("Cell populations:", length(unique(meta_selected$cell_type)), "\n")
cat("ECM genes:", length(available_genes), "\n")
cat("Pseudobulk groups:", ncol(pseudobulk_counts), "\n")

print(celltype_sample_summary)
print(expression_summary)
