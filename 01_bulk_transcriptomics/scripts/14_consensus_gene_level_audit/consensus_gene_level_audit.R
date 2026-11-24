# ============================================================
# 14_consensus_gene_level_audit.R
# Audit duplicate gene records and validation probe agreement
# ============================================================

library(dplyr)

consensus_file <- "results/leading_edge_analysis/unique_three_way_consensus_genes.csv"
modelA_file <- "results/gene_annotation/GSE184869_annotated_DEG_results.csv"
modelB_file <- "results/purity_sensitivity/Model_B_purity_adjusted_DE.csv"
validation_file <- "results/GSE125989_validation/GSE125989_paired_limma_annotated.csv"

outdir <- "results/consensus_DE_integration"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# 1. Consensus genes
consensus <- read.csv(consensus_file, check.names = FALSE) %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  distinct(gene_symbol, .keep_all = TRUE)

# 2. Model A: audit multiple Ensembl records per symbol
A_raw <- read.csv(modelA_file, check.names = FALSE)

A_audit <- A_raw %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  group_by(gene_symbol) %>%
  summarise(
    A_records = n(),
    A_direction_count = n_distinct(sign(logFC)),
    .groups = "drop"
  )

A_best <- A_raw %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  group_by(gene_symbol) %>%
  slice_min(adj.P.Val, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    gene_symbol,
    logFC_A = logFC,
    P_A = P.Value,
    FDR_A = adj.P.Val
  )

A <- A_best %>%
  left_join(A_audit, by = "gene_symbol")

# 3. Model B: map Ensembl IDs to gene symbols
annotation <- A_raw %>%
  filter(
    !is.na(ensembl_gene_id),
    !is.na(gene_symbol),
    gene_symbol != ""
  ) %>%
  select(ensembl_gene_id, gene_symbol) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

B_raw <- read.csv(modelB_file, check.names = FALSE) %>%
  left_join(annotation, by = "ensembl_gene_id") %>%
  filter(!is.na(gene_symbol), gene_symbol != "")

B_audit <- B_raw %>%
  group_by(gene_symbol) %>%
  summarise(
    B_records = n(),
    B_direction_count = n_distinct(sign(logFC)),
    .groups = "drop"
  )

B_best <- B_raw %>%
  group_by(gene_symbol) %>%
  slice_min(adj.P.Val, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    gene_symbol,
    logFC_B = logFC,
    P_B = P.Value,
    FDR_B = adj.P.Val
  )

B <- B_best %>%
  left_join(B_audit, by = "gene_symbol")

# 4. Validation: audit all probes per gene
V_raw <- read.csv(validation_file, check.names = FALSE) %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  mutate(gene_symbol = trimws(SYMBOL))

V_audit <- V_raw %>%
  group_by(gene_symbol) %>%
  summarise(
    validation_n_probes = n(),
    validation_direction_count = n_distinct(sign(logFC)),
    validation_positive_probes = sum(logFC > 0, na.rm = TRUE),
    validation_negative_probes = sum(logFC < 0, na.rm = TRUE),
    .groups = "drop"
  )

V_best <- V_raw %>%
  group_by(gene_symbol) %>%
  slice_min(adj.P.Val, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    gene_symbol,
    validation_probe = probe_id,
    logFC_validation = logFC,
    P_validation = P.Value,
    FDR_validation = adj.P.Val
  )

V <- V_best %>%
  left_join(V_audit, by = "gene_symbol")

# 5. Integrate and apply strict consistency criteria
result <- consensus %>%
  left_join(A, by = "gene_symbol") %>%
  left_join(B, by = "gene_symbol") %>%
  left_join(V, by = "gene_symbol") %>%
  mutate(
    direction_A = sign(logFC_A),
    direction_B = sign(logFC_B),
    direction_validation = sign(logFC_validation),

    representative_direction_consistent = !is.na(direction_A) &
      !is.na(direction_B) &
      !is.na(direction_validation) &
      direction_A == direction_B &
      direction_A == direction_validation,

    all_records_direction_consistent = A_direction_count == 1 &
      B_direction_count == 1 &
      validation_direction_count == 1,

    strict_direction_consistent = representative_direction_consistent &
      all_records_direction_consistent,

    significant_all_three = !is.na(FDR_A) &
      FDR_A < 0.05 &
      !is.na(FDR_B) &
      FDR_B < 0.05 &
      !is.na(FDR_validation) &
      FDR_validation < 0.05,

    strict_robust_candidate = strict_direction_consistent &
      significant_all_three
  )

# 6. Save outputs
write.csv(
  result,
  file.path(outdir, "consensus_gene_level_audit.csv"),
  row.names = FALSE
)

write.csv(
  result %>% filter(strict_direction_consistent),
  file.path(outdir, "strict_direction_consistent_genes.csv"),
  row.names = FALSE
)

write.csv(
  result %>% filter(strict_robust_candidate),
  file.path(outdir, "strict_robust_candidates.csv"),
  row.names = FALSE
)

write.csv(
  result %>% filter(validation_direction_count > 1),
  file.path(outdir, "validation_probe_discordant_genes.csv"),
  row.names = FALSE
)

# 7. Print diagnostics
cat("\n===== GENE-LEVEL AUDIT =====\n")
cat("Consensus genes:", nrow(result), "\n")
cat(
  "Genes matched in all three analyses:",
  sum(complete.cases(result[, c(
    "logFC_A",
    "logFC_B",
    "logFC_validation"
  )])),
  "\n"
)

cat(
  "Representative directions consistent:",
  sum(result$representative_direction_consistent, na.rm = TRUE),
  "\n"
)

cat(
  "Strict direction consistency across all records:",
  sum(result$strict_direction_consistent, na.rm = TRUE),
  "\n"
)

cat(
  "Significant in all three analyses:",
  sum(result$significant_all_three, na.rm = TRUE),
  "\n"
)

cat(
  "Strict robust candidates:",
  sum(result$strict_robust_candidate, na.rm = TRUE),
  "\n"
)

cat("\nValidation genes with discordant probe directions:\n")
print(
  result %>%
    filter(validation_direction_count > 1) %>%
    select(
      gene_symbol,
      validation_n_probes,
      validation_positive_probes,
      validation_negative_probes
    )
)

cat("\nStrict robust candidates:\n")
print(
  result %>%
    filter(strict_robust_candidate) %>%
    select(
      gene_symbol,
      logFC_A,
      FDR_A,
      logFC_B,
      FDR_B,
      logFC_validation,
      FDR_validation,
      validation_n_probes
    )
)

cat("\nOutputs saved in:", outdir, "\n")
