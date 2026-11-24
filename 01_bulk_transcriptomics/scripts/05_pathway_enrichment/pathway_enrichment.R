# ============================================
# BCBM PROJECT
# 05_pathway_enrichment.R
# ============================================

# 1. Install packages
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
    BiocManager::install(pkg, ask = FALSE)
  }
}

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)

# 2. File paths
project_dir <- "BCBM-Transcriptomic-Analysis"

input_file <- file.path(
  project_dir,
  "results",
  "gene_annotation",
  "GSE184869_annotated_DEG_results.csv"
)

output_dir <- file.path(
  project_dir,
  "results",
  "pathway_enrichment"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 3. Load annotated differential expression results
deg <- read.csv(input_file, stringsAsFactors = FALSE)

# 4. Prepare Ensembl IDs
deg$ensembl_gene_id <- sub(
  "\\.[0-9]+$",
  "",
  deg$ensembl_gene_id
)

deg <- deg[
  !is.na(deg$ensembl_gene_id) &
    !is.na(deg$t),
]

# 5. Convert Ensembl IDs to Entrez IDs
gene_mapping <- bitr(
  unique(deg$ensembl_gene_id),
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# 6. Merge mapping with differential expression results
mapped_deg <- merge(
  deg,
  gene_mapping,
  by.x = "ensembl_gene_id",
  by.y = "ENSEMBL"
)

# 7. Remove duplicate Entrez IDs
# Retain the gene with the strongest absolute t-statistic
mapped_deg <- mapped_deg[
  order(abs(mapped_deg$t), decreasing = TRUE),
]

mapped_deg <- mapped_deg[
  !duplicated(mapped_deg$ENTREZID),
]

# 8. Create ranked gene list
gene_list <- mapped_deg$t
names(gene_list) <- mapped_deg$ENTREZID

gene_list <- sort(gene_list, decreasing = TRUE)

cat("Genes entering GSEA:", length(gene_list), "\n")

# 9. Gene Ontology Biological Process GSEA
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

# 10. KEGG pathway GSEA
gsea_kegg <- gseKEGG(
  geneList = gene_list,
  organism = "hsa",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = FALSE
)

# 11. Convert Entrez IDs to readable gene symbols
gsea_go_readable <- setReadable(
  gsea_go,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

gsea_kegg_readable <- setReadable(
  gsea_kegg,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

# 12. Save results
write.csv(
  as.data.frame(gsea_go_readable),
  file.path(output_dir, "GSEA_GO_Biological_Process.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(gsea_kegg_readable),
  file.path(output_dir, "GSEA_KEGG.csv"),
  row.names = FALSE
)

# 13. Display summaries
cat("\nGO Biological Process results:\n")
cat("Significant pathways:", nrow(as.data.frame(gsea_go)), "\n")

if (nrow(as.data.frame(gsea_go)) > 0) {
  print(
    head(
      as.data.frame(gsea_go_readable)[,
        c("ID", "Description", "NES", "p.adjust")
      ],
      15
    )
  )
}

cat("\nKEGG results:\n")
cat("Significant pathways:", nrow(as.data.frame(gsea_kegg)), "\n")

if (nrow(as.data.frame(gsea_kegg)) > 0) {
  print(
    head(
      as.data.frame(gsea_kegg_readable)[,
        c("ID", "Description", "NES", "p.adjust")
      ],
      15
    )
  )
}

cat("\nMILESTONE 5 COMPLETED\n")
