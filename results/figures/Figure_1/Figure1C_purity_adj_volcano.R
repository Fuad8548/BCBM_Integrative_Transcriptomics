# ==========================================
# FIGURE 1C | PURITY-ADJUSTED VOLCANO
# ==========================================

library(ggplot2)
library(ggrepel)

results_B <- read.csv(file.choose(), check.names = FALSE)

results_B$FDR <- results_B$adj.P.Val
results_B$gene <- results_B$ensembl_gene_id
results_B$plot_FDR <- pmax(results_B$FDR, 1e-300)

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

results_B$gene_label <- gene_map[results_B$gene]

results_B$status <- "Not significant"
results_B$status[results_B$FDR < 0.05 & results_B$logFC > 0] <- "Upregulated"
results_B$status[results_B$FDR < 0.05 & results_B$logFC < 0] <- "Downregulated"

label_data_B <- subset(
  results_B,
  gene %in% names(gene_map) & FDR < 0.05 & !is.na(gene_label)
)

p_C <- ggplot(
  results_B,
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
    data = label_data_B,
    aes(label = gene_label),
    size = 3.5,
    fontface = "italic",
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
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
    title = "Purity-adjusted differential expression",
    subtitle = "GSE184869 | Adjusted for tumour purity",
    x = "log2 fold-change (brain metastasis vs primary tumour)",
    y = expression(-log[10](FDR)),
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    legend.position = "bottom"
  )

print(p_C)

ggsave(
  "Figure1C_purity_adjusted_volcano.png",
  p_C,
  width = 9,
  height = 7,
  dpi = 600
)
