# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 04a_prepare_doubletfinder.R
# Purpose: Prepare Seurat v5 object for DoubletFinder
# ============================================================

library(Seurat)

project_dir <- "~/Projects/SingleCell_BreastCancer_AI"

processed_dir <- file.path(
  project_dir,
  "data",
  "processed"
)

# ------------------------------------------------------------
# 1. Load QC-filtered object
# ------------------------------------------------------------

combined <- readRDS(
  file.path(
    processed_dir,
    "GSE228499_filtered.rds"
  )
)

cat("\nOriginal dimensions:\n")
print(dim(combined))

cat("\nOriginal RNA layers:\n")
print(Layers(combined[["RNA"]]))


# ------------------------------------------------------------
# 2. Join Seurat v5 layers
# ------------------------------------------------------------

cat("\nJoining RNA layers...\n")

combined <- JoinLayers(
  object = combined,
  assay = "RNA"
)


# ------------------------------------------------------------
# 3. Verify joined object
# ------------------------------------------------------------

cat("\nRNA layers after JoinLayers:\n")
print(Layers(combined[["RNA"]]))

cat("\nDimensions after JoinLayers:\n")
print(dim(combined))


# ------------------------------------------------------------
# 4. Save joined object
# ------------------------------------------------------------

joined_file <- file.path(
  processed_dir,
  "GSE228499_joined.rds"
)

saveRDS(
  combined,
  joined_file
)

cat("\n========================================\n")
cat("Joined object saved successfully\n")
cat("========================================\n")

cat("\nFile:\n")
cat(joined_file, "\n")