# ============================================
# BCBM PROJECT
# 02_data_qc.R
# ============================================

library(readxl)

# 1. File paths
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

# 2. Load expression matrix
expr <- read_excel(expression_file)
expr <- as.data.frame(expr)

rownames(expr) <- expr$ensembl_gene_id
expr$ensembl_gene_id <- NULL

# 3. Load metadata
meta <- read_excel(metadata_file, skip = 2)
meta <- as.data.frame(meta)

# 4. Standardise sample IDs
colnames(expr) <- gsub("_RCS$", "", colnames(expr))

meta$sampleID <- as.character(meta$sampleID)
meta$patientID <- as.character(meta$patientID)

# 5. Basic QC
cat("Expression matrix:", dim(expr), "\n")
cat("Metadata:", dim(meta), "\n")

cat("Duplicate gene IDs:", sum(duplicated(rownames(expr))), "\n")
cat("Duplicate expression sample IDs:", sum(duplicated(colnames(expr))), "\n")
cat("Duplicate metadata sample IDs:", sum(duplicated(meta$sampleID)), "\n")

cat("Missing expression values:", sum(is.na(expr)), "\n")
cat("Infinite expression values:", sum(!is.finite(as.matrix(expr))), "\n")

# 6. Sample matching
cat("\nExpression samples:", length(colnames(expr)), "\n")
cat("Metadata samples:", length(meta$sampleID), "\n")

cat("Missing from metadata:\n")
print(setdiff(colnames(expr), meta$sampleID))

cat("Missing from expression:\n")
print(setdiff(meta$sampleID, colnames(expr)))

# 7. Reorder metadata
stopifnot(setequal(colnames(expr), meta$sampleID))

meta <- meta[match(colnames(expr), meta$sampleID), ]

# 8. Patient pairing
pair_table <- reshape(
  meta[, c("patientID", "tumourID", "sampleID")],
  idvar = "patientID",
  timevar = "tumourID",
  direction = "wide"
)

cat("\nTotal patients:", nrow(pair_table), "\n")

complete_pairs <- sum(
  !is.na(pair_table$sampleID.P) &
    !is.na(pair_table$sampleID.M)
)

cat("Complete matched pairs:", complete_pairs, "\n")

# 9. Subtype distribution
cat("\nPrimary tumour subtype distribution:\n")
print(table(meta$IHC_subtype[meta$tumourID == "P"]))

# 10. Final validation
stopifnot(nrow(expr) == 16714)
stopifnot(ncol(expr) == 90)
stopifnot(nrow(pair_table) == 45)
stopifnot(complete_pairs == 45)
stopifnot(sum(is.na(expr)) == 0)
stopifnot(setequal(colnames(expr), meta$sampleID))

cat("\nMILESTONE 2 PASSED\n")
