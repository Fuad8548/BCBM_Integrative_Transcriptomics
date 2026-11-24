# ============================================
# BCBM PROJECT
# 04_gene_annotation.R
# ============================================

# 1. Install required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in c("AnnotationDbi", "org.Hs.eg.db")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}

library(AnnotationDbi)
library(org.Hs.eg.db)

# 2. File paths
project_dir <- "BCBM-Transcriptomic-Analysis"

input_file <- file.path(
  project_dir,
  "results",
  "differential_expression",
  "GSE184869_paired_limma_results.csv"
)

output_dir <- file.path(
  project_dir,
  "results",
  "gene_annotation"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 3. Load differential expression results
deg <- read.csv(input_file, stringsAsFactors = FALSE)

# 4. Remove Ensembl version suffix if present
deg$ensembl_gene_id <- sub(
  "\\.[0-9]+$",
  "",
  deg$ensembl_gene_id
)

# 5. Map Ensembl IDs to gene symbols and gene names
gene_ids <- unique(deg$ensembl_gene_id)

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = gene_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

gene_names <- mapIds(
  org.Hs.eg.db,
  keys = gene_ids,
  column = "GENENAME",
  keytype = "ENSEMBL",
  multiVals = "first"
)

# 6. Add annotations
deg$gene_symbol <- unname(
  gene_symbols[deg$ensembl_gene_id]
)

deg$gene_name <- unname(
  gene_names[deg$ensembl_gene_id]
)

# 7. Annotation summary
cat("Total genes:", nrow(deg), "\n")

cat(
  "Genes with symbols:",
  sum(!is.na(deg$gene_symbol)),
  "\n"
)

cat(
  "Genes without symbols:",
  sum(is.na(deg$gene_symbol)),
  "\n"
)

# 8. Define candidate genes
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

candidate_results <- deg[
  deg$gene_symbol %in% candidate_genes,
  c(
    "ensembl_gene_id",
    "gene_symbol",
    "gene_name",
    "logFC",
    "P.Value",
    "adj.P.Val",
    "regulation"
  )
]

candidate_results <- candidate_results[
  match(candidate_genes, candidate_results$gene_symbol),
]

# 9. Save complete annotated results
write.csv(
  deg,
  file.path(output_dir, "GSE184869_annotated_DEG_results.csv"),
  row.names = FALSE
)

# 10. Save candidate-gene results
write.csv(
  candidate_results,
  file.path(output_dir, "BCBM_candidate_gene_results.csv"),
  row.names = FALSE
)

# 11. Display candidate results
cat("\nCANDIDATE GENE RESULTS\n")
print(candidate_results)

cat("\nTop 20 annotated genes:\n")

top20 <- deg[
  order(deg$adj.P.Val),
  c(
    "gene_symbol",
    "gene_name",
    "logFC",
    "P.Value",
    "adj.P.Val"
  )
]

print(head(top20, 20))

cat("\nMILESTONE 4 COMPLETED\n")
