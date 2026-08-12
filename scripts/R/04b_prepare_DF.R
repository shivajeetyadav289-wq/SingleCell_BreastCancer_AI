# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 04b_prepare_DF.R
# Purpose: Normalize and perform PCA before DoubletFinder
# ============================================================

library(Seurat)

project_dir <- "~/Projects/SingleCell_BreastCancer_AI"

processed_dir <- file.path(
  project_dir,
  "data",
  "processed"
)

figure_dir <- file.path(
  project_dir,
  "figures",
  "doublet_detection"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Load joined object
# ------------------------------------------------------------

combined <- readRDS(
  file.path(
    processed_dir,
    "GSE228499_joined.rds"
  )
)

cat("\nLoaded object:\n")
print(dim(combined))

cat("\nRNA layers:\n")
print(Layers(combined[["RNA"]]))


# ------------------------------------------------------------
# 2. Normalize
# ------------------------------------------------------------

cat("\nRunning normalization...\n")

combined <- NormalizeData(
  combined,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


# ------------------------------------------------------------
# 3. Find variable genes
# ------------------------------------------------------------

cat("\nFinding variable features...\n")

combined <- FindVariableFeatures(
  combined,
  selection.method = "vst",
  nfeatures = 2000
)

cat(
  "\nVariable features:",
  length(VariableFeatures(combined)),
  "\n"
)


# ------------------------------------------------------------
# 4. Scale
# ------------------------------------------------------------

cat("\nScaling data...\n")

combined <- ScaleData(
  combined,
  features = VariableFeatures(combined)
)


# ------------------------------------------------------------
# 5. PCA
# ------------------------------------------------------------

cat("\nRunning PCA...\n")

combined <- RunPCA(
  combined,
  features = VariableFeatures(combined),
  npcs = 30,
  verbose = FALSE
)


# ------------------------------------------------------------
# 6. Save elbow plot
# ------------------------------------------------------------

pdf(
  file.path(
    figure_dir,
    "PCA_elbow_plot.pdf"
  ),
  width = 7,
  height = 5
)

print(
  ElbowPlot(
    combined,
    ndims = 30
  )
)

dev.off()


# ------------------------------------------------------------
# 7. Save DoubletFinder-ready checkpoint
# ------------------------------------------------------------

output_file <- file.path(
  processed_dir,
  "GSE228499_DF_ready.rds"
)

saveRDS(
  combined,
  output_file
)


# ------------------------------------------------------------
# 8. Final confirmation
# ------------------------------------------------------------

cat("\n========================================\n")
cat("DoubletFinder-ready object saved\n")
cat("========================================\n")

cat("\nDimensions:\n")
print(dim(combined))

cat("\nRNA layers:\n")
print(Layers(combined[["RNA"]]))

cat("\nReductions:\n")
print(Reductions(combined))

cat("\nFile:\n")
cat(output_file, "\n")