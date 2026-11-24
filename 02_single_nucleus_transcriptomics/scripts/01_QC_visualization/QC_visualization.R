# ============================================================
# GSE324453 | Step 3: QC visualization
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

bc_bm <- readRDS(
  "data/raw/GSE324453_seurat_filtered_annotated.rds"
)

meta <- bc_bm@meta.data

dir.create("results/quality_control", recursive = TRUE, showWarnings = FALSE)

dir.create("figures/quality_control", recursive = TRUE, showWarnings = FALSE)

# 1. Global QC distributions
p1 <- VlnPlot(
  bc_bm,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

ggsave(
  "figures/quality_control/global_QC_distribution.png",
  p1,
  width = 14,
  height = 5,
  dpi = 300
)

# 2. QC by patient
p2 <- VlnPlot(
  bc_bm,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  group.by = "sampleDemuxIDstrict",
  ncol = 1,
  pt.size = 0
)

ggsave(
  "figures/quality_control/patient_QC_distribution.png",
  p2,
  width = 16,
  height = 12,
  dpi = 300
)

# 3. QC by existing cell annotation
p3 <- VlnPlot(
  bc_bm,
  features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
  group.by = "common_cell_clusters",
  ncol = 3,
  pt.size = 0
)

ggsave(
  "figures/quality_control/cluster_QC_distribution.png",
  p3,
  width = 16,
  height = 7,
  dpi = 300
)

# 4. QC scatter plots
p4 <- FeatureScatter(
  bc_bm,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

p5 <- FeatureScatter(
  bc_bm,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)

ggsave(
  "figures/quality_control/UMI_vs_genes.png",
  p4,
  width = 7,
  height = 6,
  dpi = 300,
  bg = "white"
)

ggsave(
  "figures/quality_control/UMI_vs_mitochondrial.png",
  p5,
  width = 7,
  height = 6,
  dpi = 300,
  bg = "white"
)

# 5. Check mitochondrial values
cat("\n--- Mitochondrial percentage summary ---\n")
print(summary(meta$percent.mt))

cat("\n--- Maximum mitochondrial percentage ---\n")
print(max(meta$percent.mt, na.rm = TRUE))

cat("\n--- Mitochondrial values above 1 ---\n")
print(sum(meta$percent.mt > 1, na.rm = TRUE))

cat("\n--- Mitochondrial values above 0.2 ---\n")
print(sum(meta$percent.mt > 0.2, na.rm = TRUE))

cat("\nQC visualization completed.\n")
