# ============================================================
# 07_purity_adjusted_GSEA.R
# GSEA using purity-adjusted differential expression results
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

packages <- c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "enrichplot"
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
  "results/purity_sensitivity",
  "Model_B_purity_adjusted_DE.csv"
)

output_dir <- file.path(
  project_dir,
  "results/purity_adjusted_GSEA"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 2. Load results
deg <- read.csv(input_file, check.names = FALSE)

deg$ensembl_gene_id <- sub(
  "\\.[0-9]+$",
  "",
  deg$ensembl_gene_id
)

deg <- deg[
  !is.na(deg$t) & is.finite(deg$t),
]

# 3. Map Ensembl IDs to Entrez IDs
gene_map <- bitr(
  unique(deg$ensembl_gene_id),
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

mapped_deg <- merge(
  deg,
  gene_map,
  by.x = "ensembl_gene_id",
  by.y = "ENSEMBL"
)

# Keep one gene per Entrez ID
mapped_deg <- mapped_deg[
  order(abs(mapped_deg$t), decreasing = TRUE),
]

mapped_deg <- mapped_deg[
  !duplicated(mapped_deg$ENTREZID),
]

gene_list <- mapped_deg$t
names(gene_list) <- mapped_deg$ENTREZID

gene_list <- sort(gene_list, decreasing = TRUE)

cat("Mapped genes:", length(gene_list), "\n")

# 4. GO Biological Process GSEA
gsea_go <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

# 5. KEGG GSEA
gsea_kegg <- gseKEGG(
  geneList = gene_list,
  organism = "hsa",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

# 6. Save results
go_df <- as.data.frame(gsea_go)
kegg_df <- as.data.frame(gsea_kegg)

write.csv(
  go_df,
  file.path(output_dir, "Purity_adjusted_GSEA_GO_BP.csv"),
  row.names = FALSE
)

write.csv(
  kegg_df,
  file.path(output_dir, "Purity_adjusted_GSEA_KEGG.csv"),
  row.names = FALSE
)

# 7. Summary
cat("\n========== PURITY-ADJUSTED GSEA ==========\n")

cat(
  "Significant GO pathways:",
  sum(go_df$p.adjust < 0.05, na.rm = TRUE),
  "\n"
)

cat(
  "Significant KEGG pathways:",
  sum(kegg_df$p.adjust < 0.05, na.rm = TRUE),
  "\n"
)

cat("\nTop 15 GO pathways:\n")

if (nrow(go_df) > 0) {
  print(
    head(
      go_df[
        order(go_df$p.adjust),
        c("ID", "Description", "NES", "p.adjust")
      ],
      15
    )
  )
}

cat("\nTop 15 KEGG pathways:\n")

if (nrow(kegg_df) > 0) {
  print(
    head(
      kegg_df[
        order(kegg_df$p.adjust),
        c("ID", "Description", "NES", "p.adjust")
      ],
      15
    )
  )
}

cat("\nAnalysis completed.\n")
