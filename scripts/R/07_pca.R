# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 07_pca.R
#
# Purpose:
#   Scale the normalized scRNA-seq data and perform PCA.
# ============================================================

library(Seurat)
library(ggplot2)


# ------------------------------------------------------------
# 1. Project directories
# ------------------------------------------------------------

project_dir <- "~/Projects/SingleCell_BreastCancer_AI"

input_file <- file.path(
  project_dir,
  "data",
  "processed",
  "GSE228499_normalized.rds"
)

output_file <- file.path(
  project_dir,
  "data",
  "processed",
  "GSE228499_pca.rds"
)

figure_dir <- file.path(
  project_dir,
  "figures",
  "pca"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Load normalized object
# ------------------------------------------------------------

cat("\nLoading normalized dataset...\n")

singlets <- readRDS(input_file)

cat(
  "Dimensions:",
  nrow(singlets),
  "genes x",
  ncol(singlets),
  "cells\n"
)


# ------------------------------------------------------------
# 3. Scale data
# ------------------------------------------------------------

cat("\nScaling data using highly variable genes...\n")

singlets <- ScaleData(
  singlets,
  features = VariableFeatures(singlets),
  verbose = TRUE
)


# ------------------------------------------------------------
# 4. Run PCA
# ------------------------------------------------------------

cat("\nRunning PCA...\n")

singlets <- RunPCA(
  singlets,
  features = VariableFeatures(singlets),
  npcs = 50,
  verbose = TRUE
)


# ------------------------------------------------------------
# 5. PCA summary
# ------------------------------------------------------------

cat("\nPCA dimensions:\n")

print(
  dim(
    Embeddings(
      singlets,
      reduction = "pca"
    )
  )
)


# ------------------------------------------------------------
# 6. Print PCA loadings
# ------------------------------------------------------------

cat("\nTop genes contributing to first PCs:\n")

print(
  print(
    singlets[["pca"]],
    dims = 1:5,
    nfeatures = 10
  )
)


# ------------------------------------------------------------
# 7. Elbow plot
# ------------------------------------------------------------

pdf(
  file.path(
    figure_dir,
    "PCA_elbow_plot.pdf"
  ),
  width = 8,
  height = 6
)

print(
  ElbowPlot(
    singlets,
    ndims = 50
  )
)

dev.off()


# ------------------------------------------------------------
# 8. PCA plot by sample
# ------------------------------------------------------------

pdf(
  file.path(
    figure_dir,
    "PCA_by_sample.pdf"
  ),
  width = 9,
  height = 7
)

print(
  DimPlot(
    singlets,
    reduction = "pca",
    group.by = "sample_id"
  )
)

dev.off()


# ------------------------------------------------------------
# 9. PCA plot by patient
# ------------------------------------------------------------

pdf(
  file.path(
    figure_dir,
    "PCA_by_patient.pdf"
  ),
  width = 9,
  height = 7
)

print(
  DimPlot(
    singlets,
    reduction = "pca",
    group.by = "patient_id"
  )
)

dev.off()


# ------------------------------------------------------------
# 10. Save PCA object
# ------------------------------------------------------------

saveRDS(
  singlets,
  output_file
)


# ------------------------------------------------------------
# 11. Final summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("PCA COMPLETED\n")
cat("====================================================\n")

cat(
  "Cells:",
  ncol(singlets),
  "\n"
)

cat(
  "Genes:",
  nrow(singlets),
  "\n"
)

cat(
  "PCA dimensions:",
  ncol(
    Embeddings(
      singlets,
      reduction = "pca"
    )
  ),
  "\n"
)

cat(
  "Output:",
  output_file,
  "\n"
)