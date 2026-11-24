library(ggplot2)
library(ggrepel)

# Load validation differential-expression results
val <- read.csv(file.choose(), check.names = FALSE)

# Inspect column names
print(names(val))

# --------------------------------------------------
# Identify essential columns
# --------------------------------------------------

get_col <- function(data, candidates) {
  hit <- candidates[candidates %in% names(data)]
  if (length(hit) == 0) {
    stop(
      "Cannot find required column. Available columns: ",
      paste(names(data), collapse = ", ")
    )
  }
  hit[1]
}

fc_col <- get_col(val, c("logFC", "log2FC", "log2FoldChange"))
fdr_col <- get_col(val, c("adj.P.Val", "FDR", "padj", "adjusted_p_value"))

gene_col <- intersect(
  c("SYMBOL", "gene_symbol", "gene", "symbol", "GeneSymbol", "gene_name"),
  names(val)
)

# --------------------------------------------------
# Prepare data
# --------------------------------------------------

val$logFC_plot <- val[[fc_col]]
val$FDR_plot <- val[[fdr_col]]

if (length(gene_col) > 0) {
  val$gene_label <- as.character(val[[gene_col[1]]])
} else {
  val$gene_label <- NA_character_
}

val$status <- "Not significant"

val$status[
  val$FDR_plot < 0.05 &
    val$logFC_plot > 1
] <- "Upregulated"

val$status[
  val$FDR_plot < 0.05 &
    val$logFC_plot < -1
] <- "Downregulated"

val$plot_FDR <- pmax(val$FDR_plot, 1e-300)

# --------------------------------------------------
# Select labels (one top probe per gene by lowest FDR)
# --------------------------------------------------

consensus_labels <- c(
  "AEBP1",
  "COL14A1",
  "COL15A1",
  "COL6A3",
  "CTSK",
  "CXCL12",
  "DPT",
  "IGF1",
  "LUM",
  "MMP2",
  "PDGFRA",
  "POSTN",
  "SULF1",
  "THBS2",
  "THBS4"
)

label_candidates <- val[
  val$status != "Not significant" &
    !is.na(val$gene_label) &
    val$gene_label %in% consensus_labels,
]

# Keep the most significant probe per gene symbol
label_candidates <- label_candidates[
  order(label_candidates$FDR_plot, label_candidates$logFC_plot),
]
label_data <- label_candidates[!duplicated(label_candidates$gene_label), ]

# --------------------------------------------------
# Volcano plot
# --------------------------------------------------

p_S1 <- ggplot(
  val,
  aes(x = logFC_plot, y = -log10(plot_FDR))
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
    size = 3.2,
    fontface = "italic",
    box.padding = 0.45,
    point.padding = 0.25,
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
    title = "Independent validation differential expression",
    subtitle = "GSE125989 | Paired primary breast tumour vs brain metastasis",
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

print(p_S1)

# Export
ggsave(
  "Supplementary_Figure_S1_validation_volcano.png",
  p_S1,
  width = 9,
  height = 7,
  dpi = 600
)
