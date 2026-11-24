# ==========================================
# FIGURE 1B | DISCOVERY VOLCANO
# ==========================================

results_A <- read.csv(file.choose(), check.names = FALSE)

results_A$FDR <- results_A$adj.P.Val
results_A$gene <- results_A$ensembl_gene_id
results_A$plot_FDR <- pmax(results_A$FDR, 1e-300)

# Ensembl-to-symbol mapping for selected genes
gene_map <- c(
  "ENSG00000131095" = "GFAP",
  "ENSG00000170801" = "HTRA3",
  "ENSG00000109061" = "MYH1",
  "ENSG00000264424" = "MYH4",
  "ENSG00000143387" = "CTSK",
  "ENSG00000006788" = "MYH13",
  "ENSG00000125414" = "MYH2",
  "ENSG00000133020" = "MYH8",
  "ENSG00000174348" = "PODN",
  "ENSG00000172061" = "LRRC15"
)

results_A$gene_label <- ifelse(
  results_A$gene %in% names(gene_map),
  gene_map[results_A$gene],
  NA
)

# Significance classification
results_A$status <- "Not significant"
results_A$status[results_A$FDR < 0.05 & results_A$logFC > 0] <- "Upregulated"
results_A$status[results_A$FDR < 0.05 & results_A$logFC < 0] <- "Downregulated"

# Only the 10 biologically selected genes are labeled
label_data <- subset(
  results_A,
  gene %in% names(gene_map) & FDR < 0.05 & !is.na(gene_label)
)

# Volcano plot
p_B <- ggplot(
  results_A,
  aes(x = logFC, y = -log10(plot_FDR))
) +
  geom_point(
    aes(color = status),
    size = 1.2,
    alpha = 0.65
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_text_repel(
    data = label_data,
    aes(label = gene_label),
    size = 3.5,
    fontface = "italic",
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = Inf,
    min.segment.length = 0,
    seed = 42
  ) +
  scale_color_manual(
    values = c(
      "Not significant" = "grey75",
      "Upregulated" = "#C43C39",
      "Downregulated" = "#2878B5"
    )
  ) +
  labs(
    title = "Discovery differential expression",
    subtitle = "GSE184869 | Paired primary breast tumour vs brain metastasis",
    x = "log2 fold-change (brain metastasis vs primary tumour)",
    y = expression(-log[10](FDR)),
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    axis.title = element_text(size = 11),
    legend.position = "bottom"
  )

print(p_B)

ggsave(
  "Figure1B_discovery_volcano.png",
  p_B,
  width = 9,
  height = 7,
  dpi = 600
)
