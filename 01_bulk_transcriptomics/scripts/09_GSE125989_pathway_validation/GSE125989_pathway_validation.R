# ============================================================
# 09_GSE125989_pathway_validation.R
# Independent pathway validation using GSE125989
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

packages <- c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

# 1. Paths
project_dir <- "BCBM-Transcriptomic-Analysis"

input_file <- file.path(
  project_dir,
  "results/GSE125989_validation",
  "GSE125989_paired_limma_annotated.csv"
)

output_dir <- file.path(
  project_dir,
  "results/GSE125989_pathway_validation"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# 2. Load differential expression results
deg <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c("SYMBOL", "t")

if (!all(required_columns %in% colnames(deg))) {
  stop("Input file must contain SYMBOL and t columns.")
}

deg <- deg[
  !is.na(deg$SYMBOL) &
    deg$SYMBOL != "" &
    !is.na(deg$t) &
    is.finite(deg$t),
]

# 3. Collapse multiple probes per gene
# Retain the probe with the strongest absolute moderated t-statistic.

deg <- deg[
  order(abs(deg$t), decreasing = TRUE),
]

gene_level <- deg[
  !duplicated(deg$SYMBOL),
  c("SYMBOL", "t", "logFC", "P.Value", "adj.P.Val")
]

cat("Unique genes after probe collapsing:", nrow(gene_level), "\n")

# 4. Map gene symbols to Entrez IDs
gene_map <- bitr(
  gene_level$SYMBOL,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

mapped <- merge(
  gene_level,
  gene_map,
  by = "SYMBOL"
)

# Retain one entry per Entrez ID
mapped <- mapped[
  order(abs(mapped$t), decreasing = TRUE),
]

mapped <- mapped[
  !duplicated(mapped$ENTREZID),
]

gene_list <- mapped$t
names(gene_list) <- mapped$ENTREZID

gene_list <- sort(
  gene_list,
  decreasing = TRUE
)

cat("Genes mapped to Entrez IDs:", length(gene_list), "\n")

if (length(gene_list) < 1000) {
  stop("Too few mapped genes. Check annotation and input data.")
}

# 5. GO Biological Process GSEA
gsea_go <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE
)

# 6. KEGG GSEA
gsea_kegg <- gseKEGG(
  geneList = gene_list,
  organism = "hsa",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE
)

# 7. Convert results to data frames
go_df <- as.data.frame(gsea_go)
kegg_df <- as.data.frame(gsea_kegg)

# 8. Save complete pathway results
write.csv(
  go_df,
  file.path(output_dir, "GSE125989_GSEA_GO_BP.csv"),
  row.names = FALSE
)

write.csv(
  kegg_df,
  file.path(output_dir, "GSE125989_GSEA_KEGG.csv"),
  row.names = FALSE
)

# 9. Extract pathways of primary interest
target_terms <- c(
  "oxidative phosphorylation",
  "extracellular matrix organization",
  "extracellular structure organization",
  "external encapsulating structure organization",
  "adaptive immune response",
  "immune response",
  "cytokine-mediated signaling pathway",
  "leukocyte activation",
  "antigen receptor-mediated signaling pathway"
)

go_targets <- go_df[
  grepl(
    paste(target_terms, collapse = "|"),
    go_df$Description,
    ignore.case = TRUE
  ),
  c("ID", "Description", "NES", "pvalue", "p.adjust")
]

kegg_targets <- kegg_df[
  grepl(
    "oxidative phosphorylation|cytokine-cytokine receptor interaction|hematopoietic cell lineage",
    kegg_df$Description,
    ignore.case = TRUE
  ),
  c("ID", "Description", "NES", "pvalue", "p.adjust")
]

write.csv(
  go_targets,
  file.path(output_dir, "GSE125989_target_GO_pathways.csv"),
  row.names = FALSE
)

write.csv(
  kegg_targets,
  file.path(output_dir, "GSE125989_target_KEGG_pathways.csv"),
  row.names = FALSE
)

# 10. Summary
cat("\n========== GSE125989 PATHWAY VALIDATION ==========\n")

cat(
  "Significant GO pathways (FDR < 0.05):",
  sum(go_df$p.adjust < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Significant KEGG pathways (FDR < 0.05):",
  sum(kegg_df$p.adjust < 0.05, na.rm = TRUE),
  "\n"
)

cat("\nTarget GO pathways:\n")

if (nrow(go_targets) > 0) {
  print(
    go_targets[
      order(go_targets$p.adjust),
    ]
  )
} else {
  cat("No target GO pathways returned.\n")
}

cat("\nTarget KEGG pathways:\n")

if (nrow(kegg_targets) > 0) {
  print(
    kegg_targets[
      order(kegg_targets$p.adjust),
    ]
  )
} else {
  cat("No target KEGG pathways returned.\n")
}

cat("\nAnalysis completed.\n")
