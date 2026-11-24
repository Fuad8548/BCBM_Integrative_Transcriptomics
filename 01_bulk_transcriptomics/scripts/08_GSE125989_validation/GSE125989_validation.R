
# ============================================================
# 08_GSE125989_validation.R
# Independent paired validation of BCBM transcriptomics
# ============================================================

# 1. Install and load packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

packages <- c(
  "GEOquery",
  "limma",
  "AnnotationDbi",
  "hgu133a2.db"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(GEOquery)
library(limma)
library(AnnotationDbi)
library(hgu133a2.db)

# 2. Project paths
project_dir <- "BCBM-Transcriptomic-Analysis"

raw_dir <- file.path(project_dir, "data/raw/GSE125989")
output_dir <- file.path(
  project_dir, "results/GSE125989_validation"
)

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 3. Download GEO series matrix
cat("Downloading GSE125989...\n")

gse <- getGEO(
  "GSE125989",
  GSEMatrix = TRUE,
  getGPL = FALSE
)

if (length(gse) != 1) {
  stop("Expected one ExpressionSet for GSE125989.")
}

eset <- gse[[1]]

# Save downloaded object
saveRDS(
  eset,
  file.path(raw_dir, "GSE125989_ExpressionSet.rds")
)

# 4. Extract expression and metadata
expr <- exprs(eset)
meta <- pData(eset)

cat("\nExpression dimensions:\n")
print(dim(expr))

cat("\nMetadata columns:\n")
print(colnames(meta))

cat("\nSample titles:\n")
print(meta$title)

# 5. Identify sample numbers and tissue groups
sample_number <- as.integer(
  sub(
    ".*EA1690_([0-9]+).*",
    "\\1",
    meta$title
  )
)

if (any(is.na(sample_number)) ||
    any(sample_number < 1 | sample_number > 32)) {
  stop("Could not reliably extract all 32 sample numbers.")
}

group <- ifelse(
  grepl("primary breast cancer", meta$title, ignore.case = TRUE),
  "P",
  ifelse(
    grepl("brain metastas", meta$title, ignore.case = TRUE),
    "M",
    NA
  )
)

if (any(is.na(group))) {
  stop("Could not identify tissue type for every sample.")
}

# GEO sample numbering:
# EA1690_01 to EA1690_16 = primary tumours
# EA1690_17 to EA1690_32 = brain metastases
#
# Pairing assumption:
# _01 with _17, _02 with _18, etc.
# Verify this numbering-based pairing before interpreting results.

patient_number <- ifelse(
  sample_number <= 16,
  sample_number,
  sample_number - 16
)

meta$sample_number <- sample_number
meta$group <- factor(group, levels = c("P", "M"))
meta$patient_id <- factor(patient_number)

meta$sampleID <- rownames(meta)

# 6. Check pairing
pair_table <- table(meta$patient_id, meta$group)

cat("\nPairing table:\n")
print(pair_table)

if (nrow(pair_table) != 16 ||
    ncol(pair_table) != 2 ||
    any(pair_table != 1)) {
  stop("Pairing check failed. Do not proceed with paired analysis.")
}

if (!all(
  meta$sample_number[meta$group == "P"] ==
    meta$sample_number[meta$group == "M"] - 16
)) {
  stop("Sample-number pairing does not match the expected pattern.")
}

cat("\nAll 16 numbered pairs passed structural checks.\n")

# 7. Expression QC
if (anyNA(expr) || any(!is.finite(expr))) {
  stop("Expression matrix contains missing or non-finite values.")
}

cat("\nExpression summary:\n")
print(summary(as.vector(expr)))

cat("\nExpression range:\n")
print(range(expr))

# GEO reports log2 MAS5 values.
# Do not apply additional normalisation.

# 8. Paired limma analysis
design <- model.matrix(
  ~ patient_id + group,
  data = meta
)

if (qr(design)$rank != ncol(design)) {
  stop("Design matrix is not full rank.")
}

fit <- lmFit(expr, design)
fit <- eBayes(fit)

results <- topTable(
  fit,
  coef = "groupM",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

results$probe_id <- rownames(results)

# 9. Probe annotation
annotation <- AnnotationDbi::select(
  hgu133a2.db,
  keys = unique(results$probe_id),
  columns = c("SYMBOL", "GENENAME"),
  keytype = "PROBEID"
)

annotation <- annotation[
  !duplicated(annotation$PROBEID),
]

results_annotated <- merge(
  results,
  annotation,
  by.x = "probe_id",
  by.y = "PROBEID",
  all.x = TRUE,
  sort = FALSE
)

# 10. Save results
write.csv(
  results_annotated,
  file.path(
    output_dir,
    "GSE125989_paired_limma_annotated.csv"
  ),
  row.names = FALSE
)

write.csv(
  meta,
  file.path(output_dir, "GSE125989_sample_metadata.csv"),
  row.names = FALSE
)

# 11. Candidate genes
candidate_genes <- c(
  "AGR2", "ESR1", "HDAC2", "MMP9",
  "SCN5A", "ERBB2", "PGR", "REST"
)

candidate_results <- results_annotated[
  results_annotated$SYMBOL %in% candidate_genes,
  c(
    "probe_id", "SYMBOL", "GENENAME",
    "logFC", "P.Value", "adj.P.Val"
  )
]

write.csv(
  candidate_results,
  file.path(
    output_dir,
    "GSE125989_candidate_gene_results.csv"
  ),
  row.names = FALSE
)

# 12. Summary
cat("\n========== GSE125989 VALIDATION ==========\n")

cat("Genes/probes analysed:", nrow(results), "\n")

cat(
  "Significant probes (FDR < 0.05):",
  sum(results$adj.P.Val < 0.05),
  "\n"
)

cat("\nTop 20 differential probes:\n")

print(
  head(
    results_annotated[
      order(results_annotated$adj.P.Val),
      c("probe_id", "SYMBOL", "logFC", "P.Value", "adj.P.Val")
    ],
    20
  )
)

cat("\nCandidate genes:\n")
print(candidate_results)

cat("\nValidation analysis completed.\n")
