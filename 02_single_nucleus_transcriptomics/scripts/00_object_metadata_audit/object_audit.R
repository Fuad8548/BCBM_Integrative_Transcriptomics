# ============================================================
# GSE324453 | Step 1: Seurat object inspection
# ============================================================

library(Seurat)
library(dplyr)

# 1. Load object
rds_file <- "data/raw/GSE324453_seurat_filtered_annotated.rds"

bc_bm <- readRDS(rds_file)

# 2. Basic object information
cat("\n--- OBJECT CLASS ---\n")
print(class(bc_bm))

cat("\n--- OBJECT SUMMARY ---\n")
print(bc_bm)

cat("\n--- NUMBER OF CELLS ---\n")
print(ncol(bc_bm))

cat("\n--- NUMBER OF FEATURES ---\n")
print(nrow(bc_bm))

# 3. Assays and available data
cat("\n--- ASSAYS ---\n")
print(Assays(bc_bm))

for (assay_name in Assays(bc_bm)) {
  cat("\n--- ASSAY:", assay_name, "---\n")

  assay_obj <- bc_bm[[assay_name]]

  # Seurat v5
  if (inherits(assay_obj, "Assay5")) {
    print(Layers(assay_obj))

    for (layer_name in Layers(assay_obj)) {
      cat(
        layer_name,
        ":",
        paste(
          dim(LayerData(bc_bm, assay = assay_name, layer = layer_name)),
          collapse = " x "
        ),
        "\n"
      )
    }
  } else {
    # Seurat v4
    print(dim(GetAssayData(
      bc_bm,
      assay = assay_name,
      slot = "counts"
    )))
    print(dim(GetAssayData(
      bc_bm,
      assay = assay_name,
      slot = "data"
    )))
  }
}

# 4. Metadata inspection
cat("\n--- METADATA COLUMNS ---\n")
print(colnames(bc_bm@meta.data))

cat("\n--- FIRST 6 METADATA ROWS ---\n")
print(head(bc_bm@meta.data))

# 5. Sample distribution
metadata <- bc_bm@meta.data

cat("\n--- POSSIBLE SAMPLE IDENTIFIERS ---\n")
sample_columns <- grep(
  "sample|patient|orig.ident|donor",
  colnames(metadata),
  ignore.case = TRUE,
  value = TRUE
)
print(sample_columns)

for (col in sample_columns) {
  cat("\nDistribution for:", col, "\n")
  print(table(metadata[[col]], useNA = "ifany"))
}

# 6. Cell annotation distribution
cat("\n--- POSSIBLE CELL ANNOTATION COLUMNS ---\n")
annotation_columns <- grep(
  "cell|cluster|type|annotation|class|identity",
  colnames(metadata),
  ignore.case = TRUE,
  value = TRUE
)
print(annotation_columns)

for (col in annotation_columns) {
  cat("\nDistribution for:", col, "\n")
  print(head(
    sort(table(metadata[[col]], useNA = "ifany"), decreasing = TRUE),
    30
  ))
}

# 7. Existing dimensional reductions
cat("\n--- REDUCTIONS ---\n")
print(Reductions(bc_bm))

# 8. Save audit report
sink("GSE324453_object_audit.txt")

cat("GSE324453 Seurat Object Audit\n")
cat("Cells:", ncol(bc_bm), "\n")
cat("Features:", nrow(bc_bm), "\n")
cat("Assays:\n")
print(Assays(bc_bm))
cat("\nMetadata columns:\n")
print(colnames(bc_bm@meta.data))
cat("\nSample identifier columns:\n")
print(sample_columns)
cat("\nAnnotation columns:\n")
print(annotation_columns)
cat("\nReductions:\n")
print(Reductions(bc_bm))
cat("\nObject summary:\n")
print(bc_bm)

sink()

cat("\nAudit completed. Report saved as GSE324453_object_audit.txt\n")


# Patient-level QC
patient_qc <- meta %>%
  group_by(sampleDemuxIDstrict) %>%
  summarise(
    nuclei = n(),
    median_UMI = median(nCount_RNA),
    median_genes = median(nFeature_RNA),
    median_mt = median(percent.mt),
    .groups = "drop"
  ) %>%
  arrange(desc(nuclei))

print(patient_qc)

# Library-level QC
library_qc <- meta %>%
  group_by(orig.ident) %>%
  summarise(
    nuclei = n(),
    median_UMI = median(nCount_RNA),
    median_genes = median(nFeature_RNA),
    median_mt = median(percent.mt),
    .groups = "drop"
  ) %>%
  arrange(desc(nuclei))

print(library_qc)

# Save
write.csv(
  patient_qc,
  "results/patient_level_QC.csv",
  row.names = FALSE
)

write.csv(
  library_qc,
  "results/library_level_QC.csv",
  row.names = FALSE
)
