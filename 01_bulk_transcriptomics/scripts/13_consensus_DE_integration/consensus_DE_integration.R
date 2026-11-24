# ============================================================
# 13_consensus_DE_integration.R
# Integrate three-way consensus leading-edge genes with DE
# ============================================================

library(dplyr)

# ----------------------------
# 1. Input and output paths
# ----------------------------

consensus_file <- "results/leading_edge_analysis/unique_three_way_consensus_genes.csv"

modelA_file <- "results/gene_annotation/GSE184869_annotated_DEG_results.csv"

modelB_file <- "results/purity_sensitivity/Model_B_purity_adjusted_DE.csv"

validation_file <- "results/GSE125989_validation/GSE125989_paired_limma_annotated.csv"

outdir <- "results/consensus_DE_integration"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)


# ----------------------------
# 2. Read consensus genes
# ----------------------------

consensus <- read.csv(
  consensus_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot("gene_symbol" %in% names(consensus))

consensus <- consensus %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  mutate(gene_symbol = trimws(gene_symbol)) %>%
  distinct(gene_symbol, .keep_all = TRUE)


# ----------------------------
# 3. GSE184869 Model A
# ----------------------------

modelA <- read.csv(
  modelA_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

modelA <- modelA %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  transmute(
    gene_symbol = trimws(gene_symbol),
    logFC_A = logFC,
    P_A = P.Value,
    FDR_A = adj.P.Val
  ) %>%
  distinct(gene_symbol, .keep_all = TRUE)


# ----------------------------
# 4. GSE184869 Model B
# Map Ensembl IDs to symbols using Model A annotation
# ----------------------------

annotation <- read.csv(
  modelA_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  select(ensembl_gene_id, gene_symbol) %>%
  filter(
    !is.na(ensembl_gene_id),
    !is.na(gene_symbol),
    gene_symbol != ""
  ) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

modelB_raw <- read.csv(
  modelB_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

modelB <- modelB_raw %>%
  left_join(annotation, by = "ensembl_gene_id") %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  transmute(
    gene_symbol = trimws(gene_symbol),
    logFC_B = logFC,
    P_B = P.Value,
    FDR_B = adj.P.Val
  ) %>%
  distinct(gene_symbol, .keep_all = TRUE)


# ----------------------------
# 5. GSE125989 validation
# Handle multiple probes per gene
# ----------------------------

validation_raw <- read.csv(
  validation_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

validation_probes <- validation_raw %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  mutate(gene_symbol = trimws(SYMBOL)) %>%
  group_by(gene_symbol) %>%
  mutate(
    n_probes = n(),
    probe_directions = n_distinct(sign(logFC))
  ) %>%
  arrange(adj.P.Val, P.Value, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    gene_symbol,
    validation_probe = probe_id,
    logFC_validation = logFC,
    P_validation = P.Value,
    FDR_validation = adj.P.Val,
    validation_n_probes = n_probes,
    validation_probe_direction_count = probe_directions
  )


# ----------------------------
# 6. Integrate all analyses
# ----------------------------

integrated <- consensus %>%
  left_join(modelA, by = "gene_symbol") %>%
  left_join(modelB, by = "gene_symbol") %>%
  left_join(validation_probes, by = "gene_symbol") %>%
  mutate(
    direction_A = sign(logFC_A),
    direction_B = sign(logFC_B),
    direction_validation = sign(logFC_validation),

    direction_consistent_all_three = !is.na(direction_A) &
      !is.na(direction_B) &
      !is.na(direction_validation) &
      direction_A == direction_B &
      direction_A == direction_validation,

    significant_A = !is.na(FDR_A) & FDR_A < 0.05,
    significant_B = !is.na(FDR_B) & FDR_B < 0.05,
    significant_validation = !is.na(FDR_validation) & FDR_validation < 0.05,

    significant_in_all_three = significant_A &
      significant_B &
      significant_validation,

    support_category = case_when(
      significant_in_all_three & direction_consistent_all_three ~
        "Significant and direction-consistent in all three",

      direction_consistent_all_three ~
        "Direction-consistent in all three",

      significant_A & significant_validation ~
        "Significant in discovery Model A and validation",

      significant_A ~
        "Significant in discovery Model A only",

      TRUE ~
        "Other / incomplete support"
    )
  ) %>%
  arrange(
    desc(direction_consistent_all_three),
    desc(significant_in_all_three),
    FDR_validation,
    FDR_A
  )


# ----------------------------
# 7. Save outputs
# ----------------------------

write.csv(
  integrated,
  file.path(outdir, "three_way_consensus_DE_integration.csv"),
  row.names = FALSE
)

write.csv(
  integrated %>% filter(direction_consistent_all_three),
  file.path(outdir, "direction_consistent_consensus_genes.csv"),
  row.names = FALSE
)

write.csv(
  integrated %>% filter(significant_in_all_three),
  file.path(outdir, "significant_all_three_consensus_genes.csv"),
  row.names = FALSE
)

write.csv(
  validation_probes,
  file.path(outdir, "validation_gene_probe_summary.csv"),
  row.names = FALSE
)


# ----------------------------
# 8. Diagnostic summary
# ----------------------------

cat("\n====================================\n")
cat("CONSENSUS-DE INTEGRATION SUMMARY\n")
cat("====================================\n")

cat("Consensus genes:", nrow(consensus), "\n")

cat("Genes matched in Model A:", sum(!is.na(integrated$logFC_A)), "\n")

cat("Genes matched in Model B:", sum(!is.na(integrated$logFC_B)), "\n")

cat(
  "Genes matched in validation:",
  sum(!is.na(integrated$logFC_validation)),
  "\n"
)

cat(
  "\nDirection-consistent across all three:",
  sum(integrated$direction_consistent_all_three),
  "\n"
)

cat("Significant in all three:", sum(integrated$significant_in_all_three), "\n")

cat("\nSupport categories:\n")
print(table(integrated$support_category))

cat("\nTop direction-consistent genes:\n")

print(
  integrated %>%
    filter(direction_consistent_all_three) %>%
    select(
      gene_symbol,
      logFC_A,
      FDR_A,
      logFC_B,
      FDR_B,
      logFC_validation,
      FDR_validation,
      validation_n_probes,
      validation_probe_direction_count
    ) %>%
    head(30)
)

cat("\nOutput directory:", outdir, "\n")
