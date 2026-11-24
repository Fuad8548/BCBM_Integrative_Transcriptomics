# ============================================================
# 06_purity_sensitivity.R
# Compare paired differential expression models:
# A: Patient-adjusted
# B: Patient + tumour-purity-adjusted
# ============================================================

required_packages <- c("readxl", "limma")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(readxl)
library(limma)

# 1. Paths
project_dir <- "BCBM-Transcriptomic-Analysis"

expression_file <- file.path(
  project_dir,
  "data/raw/GSE184869",
  "GSE184869_rna_seq_batch_corrected_log2_TMM_normalised_CPM_protein_coding_genes.xlsx"
)

metadata_file <- file.path(
  project_dir,
  "data/raw/GSE184869",
  "41467_2022_27987_MOESM4_ESM.xlsx"
)

output_dir <- file.path(
  project_dir,
  "results/purity_sensitivity"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Load expression data
expr_df <- read_excel(expression_file)

expr_df <- as.data.frame(expr_df)

gene_ids <- expr_df[[1]]
expr_df <- expr_df[, -1, drop = FALSE]

expr <- as.matrix(
  data.frame(lapply(expr_df, as.numeric), check.names = FALSE)
)

rownames(expr) <- gene_ids
colnames(expr) <- sub("_RCS$", "", colnames(expr))

# 3. Load metadata
meta <- read_excel(metadata_file, skip = 2)
meta <- as.data.frame(meta)

# Check required metadata columns
required_cols <- c("sampleID", "patientID", "tumourID")

missing_cols <- setdiff(required_cols, colnames(meta))

if (length(missing_cols) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# Identify tumour-purity column
purity_col <- grep(
  "purity",
  colnames(meta),
  ignore.case = TRUE,
  value = TRUE
)

if (length(purity_col) != 1) {
  stop(
    "Could not identify exactly one tumour-purity column. Available columns:\n",
    paste(colnames(meta), collapse = ", ")
  )
}

cat("Purity column identified:", purity_col, "\n")

# 4. Standardise sample identifiers
meta$sampleID <- sub("_RCS$", "", as.character(meta$sampleID))

meta$patientID <- as.factor(meta$patientID)

meta$tumourID <- factor(
  meta$tumourID,
  levels = c("P", "M")
)

# Convert purity to numeric
meta$purity <- as.numeric(
  gsub("%", "", as.character(meta[[purity_col]]))
)

# 5. Match expression and metadata
common_samples <- intersect(colnames(expr), meta$sampleID)

if (length(common_samples) == 0) {
  stop("No matching sample IDs between expression and metadata.")
}

expr <- expr[, common_samples, drop = FALSE]

meta <- meta[
  match(common_samples, meta$sampleID),
  ,
  drop = FALSE
]

stopifnot(all(colnames(expr) == meta$sampleID))

# 6. Retain samples with complete purity information
meta <- meta[!is.na(meta$purity), , drop = FALSE]

expr <- expr[, meta$sampleID, drop = FALSE]

# Retain only patients with both primary and metastatic samples
pair_table <- table(meta$patientID, meta$tumourID)

complete_patients <- rownames(pair_table)[
  pair_table[, "P"] == 1 &
    pair_table[, "M"] == 1
]

meta <- meta[meta$patientID %in% complete_patients, , drop = FALSE]

meta$patientID <- droplevels(meta$patientID)

expr <- expr[, meta$sampleID, drop = FALSE]

cat("\nSamples retained:", ncol(expr), "\n")
cat("Complete patient pairs:", length(complete_patients), "\n")
cat("Purity summary:\n")
print(summary(meta$purity))

if (length(complete_patients) < 10) {
  stop("Too few complete patient pairs for a reliable sensitivity analysis.")
}

# 7. Model A: paired model without purity
design_A <- model.matrix(
  ~ patientID + tumourID,
  data = meta
)

if (qr(design_A)$rank != ncol(design_A)) {
  stop("Model A design matrix is not full rank.")
}

fit_A <- lmFit(expr, design_A)
fit_A <- eBayes(fit_A)

results_A <- topTable(
  fit_A,
  coef = "tumourIDM",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)

results_A$ensembl_gene_id <- rownames(results_A)

# 8. Model B: paired model adjusted for purity
design_B <- model.matrix(
  ~ patientID + tumourID + purity,
  data = meta
)

if (qr(design_B)$rank != ncol(design_B)) {
  stop(
    "Model B is not full rank. Tumour purity may be collinear ",
    "with patient or metastatic status."
  )
}

fit_B <- lmFit(expr, design_B)
fit_B <- eBayes(fit_B)

results_B <- topTable(
  fit_B,
  coef = "tumourIDM",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)

results_B$ensembl_gene_id <- rownames(results_B)

# 9. Save individual model results
write.csv(
  results_A,
  file.path(output_dir, "Model_A_paired_DE.csv"),
  row.names = FALSE
)

write.csv(
  results_B,
  file.path(output_dir, "Model_B_purity_adjusted_DE.csv"),
  row.names = FALSE
)

# 10. Compare models
comparison <- merge(
  results_A[, c(
    "ensembl_gene_id",
    "logFC",
    "P.Value",
    "adj.P.Val"
  )],
  results_B[, c(
    "ensembl_gene_id",
    "logFC",
    "P.Value",
    "adj.P.Val"
  )],
  by = "ensembl_gene_id",
  suffixes = c("_A", "_B")
)

comparison$logFC_change <- comparison$logFC_B -
  comparison$logFC_A

comparison$direction_consistent <- sign(comparison$logFC_A) ==
  sign(comparison$logFC_B)

comparison$significant_A <- comparison$adj.P.Val_A < 0.05
comparison$significant_B <- comparison$adj.P.Val_B < 0.05

write.csv(
  comparison,
  file.path(output_dir, "Model_A_vs_Model_B_comparison.csv"),
  row.names = FALSE
)

# 11. Summary
cat("\n========== MODEL COMPARISON ==========\n")

cat(
  "Model A significant genes:",
  sum(results_A$adj.P.Val < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Model B significant genes:",
  sum(results_B$adj.P.Val < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Genes significant in both models:",
  sum(comparison$significant_A & comparison$significant_B),
  "\n"
)

cat(
  "Genes retaining direction:",
  sum(comparison$direction_consistent, na.rm = TRUE),
  "of",
  nrow(comparison),
  "\n"
)

cat(
  "LogFC correlation:",
  cor(comparison$logFC_A, comparison$logFC_B, use = "complete.obs"),
  "\n"
)

# 12. Candidate-gene comparison
candidate_genes <- c(
  "AGR2",
  "ESR1",
  "HDAC2",
  "MMP9",
  "SCN5A",
  "ERBB2",
  "PGR",
  "REST"
)

# Ensembl-to-symbol mapping
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
  BiocManager::install("AnnotationDbi", ask = FALSE)
}

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db", ask = FALSE)
}

library(AnnotationDbi)
library(org.Hs.eg.db)

gene_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = sub("\\.[0-9]+$", "", comparison$ensembl_gene_id),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

comparison$gene_symbol <- unname(gene_map)

candidate_results <- comparison[
  comparison$gene_symbol %in% candidate_genes,
  c(
    "gene_symbol",
    "logFC_A",
    "adj.P.Val_A",
    "logFC_B",
    "adj.P.Val_B",
    "direction_consistent"
  )
]

write.csv(
  candidate_results,
  file.path(output_dir, "Candidate_gene_purity_comparison.csv"),
  row.names = FALSE
)

cat("\n========== CANDIDATE GENES ==========\n")
print(candidate_results)

cat("\nAnalysis completed.\n")
