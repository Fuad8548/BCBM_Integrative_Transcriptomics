library(dplyr)
library(ggplot2)
library(tidyr)

# ============================================================
# 1. Load existing ECM score file
# ============================================================

scores <- read.csv(file.choose(), check.names = FALSE)


# 1. Rebuild stromal ECM-high/low groups
core_scores <- scores %>%
  filter(
    cell_type == "Stromal cells",
    signature == "Core_stromal_ECM"
  ) %>%
  select(sample, signature_score)

cutoff <- median(core_scores$signature_score, na.rm = TRUE)

core_scores <- core_scores %>%
  mutate(group = ifelse(signature_score > cutoff, "High", "Low"))

print(table(core_scores$group))

# 2. ECM genes and core signature
ecm_pb <- read.csv(file.choose(), check.names = FALSE)

ecm_genes <- c(
  "AEBP1",
  "COL14A1",
  "COL15A1",
  "COL6A3",
  "CTSK",
  "CXCL12",
  "DPT",
  "LUM",
  "MMP2",
  "POSTN",
  "SULF1",
  "THBS2",
  "THBS4"
)

core_genes <- c(
  "COL6A3",
  "AEBP1",
  "SULF1",
  "THBS2",
  "LUM",
  "POSTN",
  "CTSK"
)

# 3. Extract stromal pseudobulk expression
stromal_ecm <- ecm_pb %>%
  filter(
    cell_type == "Stromal cells",
    gene %in% ecm_genes
  ) %>%
  inner_join(core_scores, by = "sample")


gene_summary <- stromal_ecm %>%
  group_by(gene, ECM_group) %>%
  summarise(
    median_logCPM = median(logCPM, na.rm = TRUE),
    mean_logCPM = mean(logCPM, na.rm = TRUE),
    n_samples = n_distinct(sample),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = ECM_group,
    values_from = c(median_logCPM, mean_logCPM, n_samples),
    names_glue = "{.value}_{ECM_group}"
  )

gene_tests <- stromal_ecm %>%
  group_by(gene) %>%
  summarise(
    p_value = {
      df_g <- cur_data()
      suppressWarnings(
        wilcox.test(
          df_g$logCPM[df_g$ECM_group == "ECM-high"],
          df_g$logCPM[df_g$ECM_group == "ECM-low"],
          exact = FALSE
        )$p.value
      )
    },
    .groups = "drop"
  ) %>%
  mutate(FDR = p.adjust(p_value, method = "BH"))

gene_results <- gene_summary %>%
  left_join(gene_tests, by = "gene") %>%
  mutate(
    difference = `median_logCPM_ECM-high` - `median_logCPM_ECM-low`,
    gene_class = ifelse(
      gene %in% core_genes,
      "Core score gene",
      "Non-core gene"
    ),
    significance = ifelse(
      FDR < 0.05,
      "FDR < 0.05",
      "Not significant"
    )
  ) %>%
  arrange(desc(difference))

p5D <- ggplot(
  gene_results,
  aes(x = difference, y = gene)
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_point(
    aes(
      color = significance,
      shape = gene_class
    ),
    size = 3.5
  ) +
  scale_color_manual(
    values = c(
      "FDR < 0.05" = "#C43C39",
      "Not significant" = "grey60"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Core score gene" = 16,
      "Non-core gene" = 17
    )
  ) +
  labs(
    title = "ECM-associated gene expression in ECM-high versus ECM-low metastases",
    subtitle = "Stromal pseudobulk; difference in median logCPM",
    x = "Median logCPM difference (ECM-high minus ECM-low)",
    y = NULL,
    color = "Statistical significance",
    shape = "Gene category"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(face = "italic"),
    legend.position = "right"
  )

print(p5D)


ggsave(
  "Figure5D_ECM_gene_comparison.png",
  p5D,
  width = 8,
  height = 6,
  dpi = 300
)
