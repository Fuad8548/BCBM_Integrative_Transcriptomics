# ============================================
# BCBM PROJECT
# 01_data_acquisition.R
# ============================================

# Install GEOquery if necessary
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("GEOquery", quietly = TRUE)) {
  BiocManager::install("GEOquery")
}

library(GEOquery)

# Project directories
project_dir <- "BCBM-Transcriptomic-Analysis"

data_dir <- file.path(project_dir, "data")
raw_dir <- file.path(data_dir, "raw")
processed_dir <- file.path(data_dir, "processed")

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

# Download GSE184869 supplementary files
getGEOSuppFiles(
  GEO = "GSE184869",
  baseDir = raw_dir,
  makeDirectory = TRUE
)

# Download GEO series metadata
gse <- getGEO(
  "GSE184869",
  GSEMatrix = TRUE,
  getGPL = FALSE
)

# Save GEO metadata
saveRDS(
  gse,
  file = file.path(processed_dir, "GSE184869_GEO_metadata.rds")
)

cat("GSE184869 acquisition completed.\n")
