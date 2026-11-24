# ============================================================
# 10_consensus_pathways.R
# Identify reproducible pathways across BCBM cohorts
# ============================================================

project_dir <- "BCBM-Transcriptomic-Analysis"

# 1. Input files
files <- list(
  GO_discovery = file.path(
    project_dir,
    "results/pathway_enrichment",
    "GSEA_GO_Biological_Process.csv"
  ),
  GO_purity = file.path(
    project_dir,
    "results/purity_adjusted_GSEA",
    "Purity_adjusted_GSEA_GO_BP.csv"
  ),
  GO_validation = file.path(
    project_dir,
    "results/GSE125989_pathway_validation",
    "GSE125989_GSEA_GO_BP.csv"
  ),
  KEGG_discovery = file.path(
    project_dir,
    "results/pathway_enrichment",
    "GSEA_KEGG.csv"
  ),
  KEGG_purity = file.path(
    project_dir,
    "results/purity_adjusted_GSEA",
    "Purity_adjusted_GSEA_KEGG.csv"
  ),
  KEGG_validation = file.path(
    project_dir,
    "results/GSE125989_pathway_validation",
    "GSE125989_GSEA_KEGG.csv"
  )
)

# 2. Check files
missing_files <- files[!file.exists(unlist(files))]

if (length(missing_files) > 0) {
  stop(
    "Missing input files:\n",
    paste(unlist(missing_files), collapse = "\n")
  )
}

# 3. Read pathway results
read_pathway <- function(path, label) {
  x <- read.csv(path, check.names = FALSE)

  required <- c("ID", "Description", "NES", "p.adjust")

  if (!all(required %in% colnames(x))) {
    stop("Missing required columns in: ", path)
  }

  x <- x[, required]

  colnames(x)[3:4] <- c(
    paste0("NES_", label),
    paste0("FDR_", label)
  )

  x
}

go_A <- read_pathway(files$GO_discovery, "discovery")
go_B <- read_pathway(files$GO_purity, "purity")
go_C <- read_pathway(files$GO_validation, "validation")

kegg_A <- read_pathway(files$KEGG_discovery, "discovery")
kegg_B <- read_pathway(files$KEGG_purity, "purity")
kegg_C <- read_pathway(files$KEGG_validation, "validation")

# 4. Merge pathways by pathway ID
merge_consensus <- function(A, B, C) {
  x <- merge(A, B, by = c("ID", "Description"))
  x <- merge(x, C, by = c("ID", "Description"))

  x$significant_all_three <-
    x$FDR_discovery < 0.05 &
    x$FDR_purity < 0.05 &
    x$FDR_validation < 0.05

  x$direction_consistent <-
    sign(x$NES_discovery) == sign(x$NES_purity) &
    sign(x$NES_purity) == sign(x$NES_validation)

  x$robust_consensus <-
    x$significant_all_three &
    x$direction_consistent

  x
}

GO_consensus <- merge_consensus(go_A, go_B, go_C)

KEGG_consensus <- merge_consensus(
  kegg_A,
  kegg_B,
  kegg_C
)

# 5. Save complete comparisons
output_dir <- file.path(
  project_dir,
  "results/consensus_pathways"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  GO_consensus,
  file.path(output_dir, "GO_consensus_all.csv"),
  row.names = FALSE
)

write.csv(
  KEGG_consensus,
  file.path(output_dir, "KEGG_consensus_all.csv"),
  row.names = FALSE
)

# 6. Extract robust consensus pathways
GO_robust <- GO_consensus[
  GO_consensus$robust_consensus,
]

KEGG_robust <- KEGG_consensus[
  KEGG_consensus$robust_consensus,
]

GO_robust <- GO_robust[
  order(GO_robust$FDR_validation),
]

KEGG_robust <- KEGG_robust[
  order(KEGG_robust$FDR_validation),
]

write.csv(
  GO_robust,
  file.path(output_dir, "GO_robust_consensus.csv"),
  row.names = FALSE
)

write.csv(
  KEGG_robust,
  file.path(output_dir, "KEGG_robust_consensus.csv"),
  row.names = FALSE
)

# 7. Summary
cat("\n========== CONSENSUS PATHWAY ANALYSIS ==========\n")

cat("GO pathways significant and direction-consistent:", nrow(GO_robust), "\n")

cat(
  "KEGG pathways significant and direction-consistent:",
  nrow(KEGG_robust),
  "\n"
)

cat("\nTop GO consensus pathways:\n")

if (nrow(GO_robust) > 0) {
  print(
    head(
      GO_robust[, c(
        "ID",
        "Description",
        "NES_discovery",
        "NES_purity",
        "NES_validation",
        "FDR_discovery",
        "FDR_purity",
        "FDR_validation"
      )],
      20
    )
  )
}

cat("\nTop KEGG consensus pathways:\n")

if (nrow(KEGG_robust) > 0) {
  print(
    head(
      KEGG_robust[, c(
        "ID",
        "Description",
        "NES_discovery",
        "NES_purity",
        "NES_validation",
        "FDR_discovery",
        "FDR_purity",
        "FDR_validation"
      )],
      20
    )
  )
}

cat("\nAnalysis completed.\n")
