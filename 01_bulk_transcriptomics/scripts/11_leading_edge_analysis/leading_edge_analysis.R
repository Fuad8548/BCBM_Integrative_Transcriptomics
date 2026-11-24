# ============================================================
# 11_leading_edge_analysis.R
# Identifier-aware leading-edge consensus analysis
# ============================================================

# 1. Packages -------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in c("AnnotationDbi", "org.Hs.eg.db")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

library(AnnotationDbi)
library(org.Hs.eg.db)

# 2. Input files ----------------------------------------------

files <- list(
  discovery_GO = "results/pathway_enrichment/GSEA_GO_Biological_Process.csv",
  purity_GO = "results/purity_adjusted_GSEA/Purity_adjusted_GSEA_GO_BP.csv",
  validation_GO = "results/GSE125989_pathway_validation/GSE125989_GSEA_GO_BP.csv",
  discovery_KEGG = "results/pathway_enrichment/GSEA_KEGG.csv",
  purity_KEGG = "results/purity_adjusted_GSEA/Purity_adjusted_GSEA_KEGG.csv",
  validation_KEGG = "results/GSE125989_pathway_validation/GSE125989_GSEA_KEGG.csv"
)

missing_files <- files[!file.exists(unlist(files))]

if (length(missing_files) > 0) {
  stop("Missing files:\n", paste(unlist(missing_files), collapse = "\n"))
}

outdir <- "results/leading_edge_analysis"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# 3. Selected pathways ----------------------------------------

pathways <- data.frame(
  ID = c(
    "hsa00190",
    "GO:0009060",
    "GO:0030198",
    "GO:0030199",
    "hsa04512",
    "hsa04518",
    "hsa04060"
  ),
  Description = c(
    "Oxidative phosphorylation",
    "Aerobic respiration",
    "Extracellular matrix organization",
    "Collagen fibril organization",
    "ECM-receptor interaction",
    "Integrin signaling",
    "Cytokine-cytokine receptor interaction"
  ),
  Type = c("KEGG", "GO", "GO", "GO", "KEGG", "KEGG", "KEGG"),
  stringsAsFactors = FALSE
)

# 4. Read GSEA files ------------------------------------------

gsea_results <- lapply(files, function(f) {
  read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
})

required <- c("ID", "Description", "NES", "p.adjust", "core_enrichment")

for (nm in names(gsea_results)) {
  missing <- setdiff(required, colnames(gsea_results[[nm]]))
  if (length(missing) > 0) {
    stop("Missing columns in ", nm, ": ", paste(missing, collapse = ", "))
  }
}

# 5. Identifier conversion ------------------------------------
# Numeric IDs are treated as Entrez IDs.
# Non-numeric IDs are treated as gene symbols.
# Unmapped values remain explicitly labelled RAW.

known_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
known_aliases <- keys(org.Hs.eg.db, keytype = "ALIAS")

convert_gene_ids <- function(ids) {
  ids <- unique(trimws(as.character(ids)))
  ids <- ids[!is.na(ids) & nzchar(ids)]

  result <- data.frame(
    original_id = ids,
    gene_symbol = NA_character_,
    gene_id = NA_character_,
    mapping_type = NA_character_,
    stringsAsFactors = FALSE
  )

  numeric_idx <- grepl("^[0-9]+$", ids)
  symbol_idx <- !numeric_idx

  # Entrez IDs
  if (any(numeric_idx)) {
    numeric_ids <- ids[numeric_idx]

    mapped <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = numeric_ids,
      column = "SYMBOL",
      keytype = "ENTREZID",
      multiVals = "first"
    )

    result$gene_symbol[numeric_idx] <- unname(mapped[numeric_ids])
    result$gene_id[numeric_idx] <- numeric_ids
    result$mapping_type[numeric_idx] <- "EntrezID"
  }

  # Gene symbols and possible aliases
  if (any(symbol_idx)) {
    symbol_ids <- ids[symbol_idx]

    direct_match <- symbol_ids %in% known_symbols

    result$gene_symbol[symbol_idx] <- symbol_ids
    result$mapping_type[symbol_idx] <- "GeneSymbol"

    # Resolve known aliases only when no direct symbol exists
    alias_idx <- !direct_match & symbol_ids %in% known_aliases

    if (any(alias_idx)) {
      alias_values <- symbol_ids[alias_idx]

      alias_map <- AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys = alias_values,
        column = "SYMBOL",
        keytype = "ALIAS",
        multiVals = "first"
      )

      positions <- which(symbol_idx)[alias_idx]

      result$gene_symbol[positions] <- unname(alias_map[alias_values])
      result$mapping_type[positions] <- "Alias"
    }

    # Retrieve Entrez IDs for recognised symbols
    recognised <- which(symbol_idx & !is.na(result$gene_symbol))

    if (length(recognised) > 0) {
      symbol_map <- AnnotationDbi::mapIds(
        org.Hs.eg.db,
        keys = unique(result$gene_symbol[recognised]),
        column = "ENTREZID",
        keytype = "SYMBOL",
        multiVals = "first"
      )

      result$gene_id[recognised] <- unname(
        symbol_map[result$gene_symbol[recognised]]
      )
    }
  }

  result$gene_key <- ifelse(
    !is.na(result$gene_symbol),
    paste0("SYMBOL:", result$gene_symbol),
    paste0("RAW:", result$original_id)
  )

  result
}

# 6. Extract leading-edge genes -------------------------------

extract_leading_edge <- function(data, pathway_id, dataset_name) {
  selected <- data[
    trimws(as.character(data$ID)) == pathway_id,
    ,
    drop = FALSE
  ]

  if (nrow(selected) == 0) {
    warning(pathway_id, " not found in ", dataset_name)
    return(NULL)
  }

  selected <- selected[1, , drop = FALSE]

  core <- as.character(selected$core_enrichment)

  if (is.na(core) || !nzchar(trimws(core))) {
    warning("Empty core_enrichment: ", pathway_id, " / ", dataset_name)
    return(NULL)
  }

  ids <- unique(trimws(unlist(strsplit(core, "/", fixed = TRUE))))
  ids <- ids[!is.na(ids) & nzchar(ids)]

  mapping <- convert_gene_ids(ids)

  data.frame(
    pathway_id = pathway_id,
    pathway = as.character(selected$Description),
    dataset = dataset_name,
    original_id = mapping$original_id,
    gene_id = mapping$gene_id,
    gene_symbol = mapping$gene_symbol,
    gene_key = mapping$gene_key,
    mapping_type = mapping$mapping_type,
    NES = as.numeric(selected$NES),
    FDR = as.numeric(selected$p.adjust),
    stringsAsFactors = FALSE
  )
}

# 7. Extract all selected pathways ----------------------------

leading_edge_list <- list()

for (i in seq_len(nrow(pathways))) {
  pathway_id <- pathways$ID[i]

  datasets <- if (pathways$Type[i] == "GO") {
    c("discovery_GO", "purity_GO", "validation_GO")
  } else {
    c("discovery_KEGG", "purity_KEGG", "validation_KEGG")
  }

  for (ds in datasets) {
    x <- extract_leading_edge(
      gsea_results[[ds]],
      pathway_id,
      ds
    )

    if (!is.null(x)) {
      leading_edge_list[[length(leading_edge_list) + 1]] <- x
    }
  }
}

if (length(leading_edge_list) == 0) {
  stop("No leading-edge genes extracted.")
}

leading_edge <- do.call(rbind, leading_edge_list)

# 8. Remove duplicate records ---------------------------------

leading_edge <- unique(leading_edge)

write.csv(
  leading_edge,
  file.path(outdir, "all_leading_edge_genes.csv"),
  row.names = FALSE
)

# 9. Pathway-level gene support -------------------------------

support_list <- list()
summary_list <- list()

for (i in seq_len(nrow(pathways))) {
  pid <- pathways$ID[i]

  pathway_data <- leading_edge[
    leading_edge$pathway_id == pid,
    ,
    drop = FALSE
  ]

  expected <- if (pathways$Type[i] == "GO") {
    c("discovery_GO", "purity_GO", "validation_GO")
  } else {
    c("discovery_KEGG", "purity_KEGG", "validation_KEGG")
  }

  # Use only successfully mapped symbols for consensus.
  mapped_data <- pathway_data[
    !is.na(pathway_data$gene_symbol),
    ,
    drop = FALSE
  ]

  gene_sets <- lapply(expected, function(ds) {
    unique(mapped_data$gene_key[mapped_data$dataset == ds])
  })

  names(gene_sets) <- c("discovery", "purity", "validation")

  all_keys <- sort(unique(unlist(gene_sets)))

  if (length(all_keys) > 0) {
    support_matrix <- sapply(
      gene_sets,
      function(x) all_keys %in% x
    )

    support_matrix <- as.matrix(support_matrix)
    colnames(support_matrix) <- names(gene_sets)

    counts <- rowSums(support_matrix)

    annotation <- mapped_data[
      match(all_keys, mapped_data$gene_key),
      c("gene_id", "gene_symbol"),
      drop = FALSE
    ]

    support_df <- data.frame(
      pathway_id = pid,
      pathway = pathways$Description[i],
      gene_key = all_keys,
      discovery = support_matrix[, "discovery"],
      purity = support_matrix[, "purity"],
      validation = support_matrix[, "validation"],
      datasets_supported = counts,
      gene_id = annotation$gene_id,
      gene_symbol = annotation$gene_symbol,
      stringsAsFactors = FALSE
    )

    support_df$support_category <- ifelse(
      counts == 3,
      "Shared_all_three",
      ifelse(counts == 2, "Shared_two", "Dataset_specific")
    )

    support_list[[length(support_list) + 1]] <- support_df
  } else {
    counts <- integer(0)
  }

  summary_list[[length(summary_list) + 1]] <- data.frame(
    pathway_id = pid,
    pathway = pathways$Description[i],
    discovery_genes = length(gene_sets$discovery),
    purity_genes = length(gene_sets$purity),
    validation_genes = length(gene_sets$validation),
    union_genes = length(all_keys),
    shared_all_three = sum(counts == 3),
    shared_at_least_two = sum(counts >= 2),
    stringsAsFactors = FALSE
  )
}

support_genes <- if (length(support_list) > 0) {
  do.call(rbind, support_list)
} else {
  data.frame()
}

summary_table <- do.call(rbind, summary_list)

# 10. Export support and consensus ----------------------------

write.csv(
  support_genes,
  file.path(outdir, "pathway_gene_support.csv"),
  row.names = FALSE
)

consensus_genes <- if (nrow(support_genes) > 0) {
  support_genes[support_genes$datasets_supported >= 2, , drop = FALSE]
} else {
  data.frame()
}

write.csv(
  consensus_genes,
  file.path(outdir, "consensus_leading_edge_genes.csv"),
  row.names = FALSE
)

write.csv(
  summary_table,
  file.path(outdir, "leading_edge_pathway_summary.csv"),
  row.names = FALSE
)

# 11. Unique gene lists ---------------------------------------

unique_all <- unique(
  leading_edge[
    !is.na(leading_edge$gene_symbol),
    c("gene_id", "gene_symbol"),
    drop = FALSE
  ]
)

unique_all <- unique_all[order(unique_all$gene_symbol), ]

write.csv(
  unique_all,
  file.path(outdir, "unique_leading_edge_gene_list.csv"),
  row.names = FALSE
)

unique_consensus <- if (nrow(consensus_genes) > 0) {
  unique(consensus_genes[, c("gene_id", "gene_symbol")])
} else {
  data.frame(
    gene_id = character(),
    gene_symbol = character()
  )
}

unique_consensus <- unique_consensus[
  order(unique_consensus$gene_symbol),
  ,
  drop = FALSE
]

write.csv(
  unique_consensus,
  file.path(outdir, "unique_consensus_gene_list.csv"),
  row.names = FALSE
)

# 12. Identifier diagnostics ----------------------------------

mapping_summary <- as.data.frame(
  table(
    leading_edge$dataset,
    leading_edge$mapping_type,
    useNA = "ifany"
  )
)

write.csv(
  mapping_summary,
  file.path(outdir, "identifier_mapping_diagnostics.csv"),
  row.names = FALSE
)

unmapped <- leading_edge[
  is.na(leading_edge$gene_symbol),
  ,
  drop = FALSE
]

write.csv(
  unmapped,
  file.path(outdir, "unmapped_leading_edge_identifiers.csv"),
  row.names = FALSE
)

# 13. Print results -------------------------------------------

cat("\n========================================\n")
cat("IDENTIFIER-AWARE LEADING-EDGE ANALYSIS\n")
cat("========================================\n")

cat("\nTotal leading-edge records:", nrow(leading_edge), "\n")

cat(
  "Unique mapped symbols:",
  length(unique(na.omit(leading_edge$gene_symbol))),
  "\n"
)

cat(
  "Unique genes across selected pathways:",
  nrow(unique_all),
  "\n"
)

cat(
  "Unique genes supported by >=2 analyses:",
  nrow(unique_consensus),
  "\n"
)

cat("\nIdentifier mapping diagnostics:\n")
print(mapping_summary)

cat("\nPathway summary:\n")
print(summary_table)

cat("\nConsensus support categories:\n")
if (nrow(support_genes) > 0) {
  print(table(support_genes$support_category))
}

cat("\nUnmapped identifier records:", nrow(unmapped), "\n")

cat("\nOutput directory:", outdir, "\n")
