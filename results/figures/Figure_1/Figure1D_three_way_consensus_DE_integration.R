library(ggplot2)

consensus <- read.csv(file.choose(), check.names = FALSE)

# Keep the 15 genes significant in all three analyses
consensus <- consensus[consensus$significant_in_all_three == TRUE, ]

# Order genes by Model A effect
consensus$gene_symbol <- factor(
  consensus$gene_symbol,
  levels = consensus$gene_symbol[order(consensus$logFC_A)]
)

# Convert to long format
plot_data <- data.frame(
  gene_symbol = rep(consensus$gene_symbol, 3),
  Analysis = rep(
    c("Model A", "Model B", "Validation"),
    each = nrow(consensus)
  ),
  logFC = c(
    consensus$logFC_A,
    consensus$logFC_B,
    consensus$logFC_validation
  )
)

plot_data$Analysis <- factor(
  plot_data$Analysis,
  levels = c("Model A", "Model B", "Validation")
)

# Panel D
p_D <- ggplot(plot_data, aes(x = logFC, y = gene_symbol, color = Analysis)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(
    position = position_dodge(width = 0.55),
    size = 2.8
  ) +
  scale_color_manual(
    values = c(
      "Model A" = "#2878B5",
      "Model B" = "#D98B32",
      "Validation" = "#38966B"
    )
  ) +
  labs(
    title = "Cross-analysis consistency of consensus genes",
    subtitle = "15 genes significant in discovery, purity-adjusted analysis, and validation",
    x = "log2 fold-change (brain metastasis vs primary tumour)",
    y = NULL,
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    axis.title.x = element_text(size = 11),
    axis.text.y = element_text(face = "italic", size = 10),
    legend.position = "bottom"
  )

print(p_D)

ggsave(
  "Figure1D_three_way_consensus_effects.png",
  p_D,
  width = 9,
  height = 7,
  dpi = 600
)
