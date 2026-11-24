# ============================================================
# 12_consensus_gene_summary.R
# ============================================================

library(dplyr)

input_file <- "results/leading_edge_analysis/pathway_gene_support.csv"

outdir <- "results/leading_edge_analysis"

genes <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Three-way shared genes, nonredundant by symbol
three_way <- genes %>%
  filter(datasets_supported == 3) %>%
  group_by(gene_symbol) %>%
  summarise(
    gene_id = first(na.omit(gene_id)),
    pathways = paste(sort(unique(pathway)), collapse = "; "),
    pathway_count = n_distinct(pathway),
    .groups = "drop"
  )

# At least two datasets, nonredundant by symbol
two_way <- genes %>%
  filter(datasets_supported >= 2) %>%
  group_by(gene_symbol) %>%
  summarise(
    gene_id = first(na.omit(gene_id)),
    pathways = paste(sort(unique(pathway)), collapse = "; "),
    pathway_count = n_distinct(pathway),
    .groups = "drop"
  )

write.csv(
  three_way,
  file.path(outdir, "unique_three_way_consensus_genes.csv"),
  row.names = FALSE
)

write.csv(
  two_way,
  file.path(outdir, "unique_two_way_consensus_genes.csv"),
  row.names = FALSE
)

cat("Unique three-way consensus genes:", nrow(three_way), "\n")
cat("Unique two-way consensus genes:", nrow(two_way), "\n")

cat("\nThree-way genes by pathway membership:\n")
print(table(three_way$pathway_count))

cat("\nTop genes appearing in multiple pathways:\n")
print(
  three_way %>%
    filter(pathway_count >= 2) %>%
    arrange(desc(pathway_count))
)
