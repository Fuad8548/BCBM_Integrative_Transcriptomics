# ============================================
# BCBM PROJECT
# 03_paired_differential_expression.R
# ============================================

# 1. Install and load limma
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("limma", quietly = TRUE)) {
  BiocManager::install("limma")
}

library(limma)
library(readxl)

# 2. File paths
project_dir <- "BCBM-Transcriptomic-Analysis"

expression_file <- file.path(
  project_dir,
  "data",
  "raw",
  "GSE184869",
  "GSE184869_rna_seq_batch_corrected_log2_TMM_normalised_CPM_protein_coding_genes.xlsx"
)

metadata_file <- file.path(
  project_dir,
  "data",
  "raw",
  "GSE184869",
  "41467_2022_27987_MOESM4_ESM.xlsx"
)

results_dir <- file.path(project_dir, "results", "differential_expression")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# 3. Load expression matrix
expr <- read_excel(expression_file)
expr <- as.data.frame(expr)

rownames(expr) <- expr$ensembl_gene_id
expr$ensembl_gene_id <- NULL

expr <- as.matrix(expr)
storage.mode(expr) <- "numeric"

# 4. Load metadata
meta <- read_excel(metadata_file, skip = 2)
meta <- as.data.frame(meta)

meta$sampleID <- as.character(meta$sampleID)
meta$patientID <- as.character(meta$patientID)

# Standardise sample IDs
colnames(expr) <- gsub("_RCS$", "", colnames(expr))

# 5. Match metadata to expression matrix
stopifnot(setequal(colnames(expr), meta$sampleID))

meta <- meta[match(colnames(expr), meta$sampleID), ]

# 6. Define paired design
meta$patientID <- factor(meta$patientID)
meta$tumourID <- factor(meta$tumourID, levels = c("P", "M"))

design <- model.matrix(
  ~ patientID + tumourID,
  data = meta
)

# 7. Fit limma model
fit <- lmFit(expr, design)
fit <- eBayes(fit)

# 8. Extract brain metastasis vs primary tumour results
coef_name <- "tumourIDM"

deg <- topTable(
  fit,
  coef = coef_name,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

deg$ensembl_gene_id <- rownames(deg)

# 9. Classify differential expression
deg$regulation <- ifelse(
  deg$adj.P.Val < 0.05 & deg$logFC > 0,
  "Up in brain metastasis",
  ifelse(
    deg$adj.P.Val < 0.05 & deg$logFC < 0,
    "Down in brain metastasis",
    "Not significant"
  )
)

# 10. Save results
write.csv(
  deg,
  file.path(results_dir, "GSE184869_paired_limma_results.csv"),
  row.names = FALSE
)

# 11. Summary
cat("Total genes analysed:", nrow(deg), "\n")
cat(
  "Significant genes (FDR < 0.05):",
  sum(deg$adj.P.Val < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Upregulated in brain metastasis:",
  sum(deg$adj.P.Val < 0.05 & deg$logFC > 0, na.rm = TRUE),
  "\n"
)

cat(
  "Downregulated in brain metastasis:",
  sum(deg$adj.P.Val < 0.05 & deg$logFC < 0, na.rm = TRUE),
  "\n"
)

cat("\nTop 20 differential expression results:\n")
print(
  head(
    deg[, c("ensembl_gene_id", "logFC", "P.Value", "adj.P.Val")],
    20
  )
)

cat("\nMILESTONE 3 COMPLETED\n")
