# ============================================================
# BCBM snRNA-seq | GSE324453
# Milestone 6: ECM detection prevalence per nucleus
# RAM-efficient: 13 genes only
# ============================================================

library(Seurat)
library(dplyr)

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

input_file <- "data/raw/GSE324453_seurat_filtered_annotated.rds"

output_dir <- "results/ECM_detection_prevalence"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Load Seurat object
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
# 3. Prepare metadata
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
    cell_type = as.character(common_cell_clusters)
  )

selected_cells <- meta_selected$cell

# ------------------------------------------------------------
# 4. Extract only 13 genes
# ------------------------------------------------------------

available_genes <- intersect(
  ecm_genes,
  rownames(bcbm)
)

ecm_counts <- LayerData(
  bcbm,
  assay = "RNA",
  layer = "counts",
  features = available_genes,
  cells = selected_cells
)

# Align metadata to matrix columns
meta_selected <- meta_selected[
  match(colnames(ecm_counts), meta_selected$cell),
]

stopifnot(
  identical(meta_selected$cell, colnames(ecm_counts))
)

# ------------------------------------------------------------
# 5. Overall detection by cell type
# ------------------------------------------------------------

overall_results <- list()

for (gene in available_genes) {
  gene_counts <- as.numeric(
    ecm_counts[gene, ]
  )

  gene_meta <- meta_selected

  gene_meta$count <- gene_counts
  gene_meta$detected <- gene_counts > 0

  result <- gene_meta %>%
    group_by(cell_type) %>%
    summarise(
      gene = gene,
      n_nuclei = n(),
      n_detected = sum(detected),
      pct_detected = mean(detected) * 100,
      mean_count_all_nuclei = mean(count),
      mean_count_positive_nuclei = ifelse(
        sum(detected) > 0,
        mean(count[detected]),
        NA_real_
      ),
      median_count_positive_nuclei = ifelse(
        sum(detected) > 0,
        median(count[detected]),
        NA_real_
      ),
      .groups = "drop"
    )

  overall_results[[gene]] <- result
}

overall_detection <- bind_rows(overall_results)

write.csv(
  overall_detection,
  file.path(output_dir, "ECM_detection_by_celltype.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 6. Sample-level detection
# ------------------------------------------------------------

sample_results <- list()

for (gene in available_genes) {
  gene_counts <- as.numeric(
    ecm_counts[gene, ]
  )

  gene_meta <- meta_selected

  gene_meta$count <- gene_counts
  gene_meta$detected <- gene_counts > 0

  result <- gene_meta %>%
    group_by(sample, cell_type) %>%
    summarise(
      gene = gene,
      n_nuclei = n(),
      n_detected = sum(detected),
      pct_detected = mean(detected) * 100,
      mean_count_all_nuclei = mean(count),
      mean_count_positive_nuclei = ifelse(
        sum(detected) > 0,
        mean(count[detected]),
        NA_real_
      ),
      .groups = "drop"
    )

  sample_results[[gene]] <- result
}

sample_detection <- bind_rows(sample_results)

write.csv(
  sample_detection,
  file.path(output_dir, "ECM_detection_by_sample_celltype.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 7. Focused malignant-cell results
# ------------------------------------------------------------

tumor_detection <- overall_detection %>%
  filter(
    cell_type %in%
      c(
        "Tumor cells",
        "Tumor cells (prol.)"
      )
  ) %>%
  arrange(gene, cell_type)

write.csv(
  tumor_detection,
  file.path(output_dir, "ECM_tumor_cell_detection.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 8. Rank cellular sources by detection prevalence
# ------------------------------------------------------------

source_ranking <- overall_detection %>%
  group_by(gene) %>%
  arrange(desc(pct_detected), .by_group = TRUE) %>%
  mutate(source_rank = row_number()) %>%
  ungroup()

write.csv(
  source_ranking,
  file.path(output_dir, "ECM_cellular_source_ranking.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. Final output
# ------------------------------------------------------------

cat("\n========================================\n")
cat("ECM DETECTION ANALYSIS COMPLETE\n")
cat("========================================\n")

print(
  overall_detection %>%
    arrange(gene, desc(pct_detected))
)

cat("\nTumour cell detection summary:\n")

print(tumor_detection)

cat("\nOutput directory:", output_dir, "\n")
