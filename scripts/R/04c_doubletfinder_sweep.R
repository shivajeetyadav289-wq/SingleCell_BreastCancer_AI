# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 04c_doubletfinder_sweep.R
# Purpose: Optimize DoubletFinder pK
# ============================================================

library(Seurat)
library(DoubletFinder)

project_dir <- "~/Projects/SingleCell_BreastCancer_AI"

processed_dir <- file.path(
  project_dir,
  "data",
  "processed"
)

doublet_dir <- file.path(
  project_dir,
  "results",
  "doublet_detection"
)

dir.create(
  doublet_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 1. Load checkpoint
# ------------------------------------------------------------

combined <- readRDS(
  file.path(
    processed_dir,
    "GSE228499_DF_ready.rds"
  )
)


# ------------------------------------------------------------
# 2. Verify object
# ------------------------------------------------------------

cat("\nDimensions:\n")
print(dim(combined))

cat("\nRNA layers:\n")
print(Layers(combined[["RNA"]]))

cat("\nReductions:\n")
print(Reductions(combined))


# ------------------------------------------------------------
# 3. Run parameter sweep
# ------------------------------------------------------------

cat("\n========================================\n")
cat("Running DoubletFinder parameter sweep\n")
cat("========================================\n")

sweep_res <- paramSweep(
  combined,
  PCs = 1:20,
  sct = FALSE
)


# ------------------------------------------------------------
# 4. Summarize sweep
# ------------------------------------------------------------

sweep_stats <- summarizeSweep(
  sweep_res,
  GT = FALSE
)


# ------------------------------------------------------------
# 5. Find optimal pK
# ------------------------------------------------------------

bcmvn <- find.pK(
  sweep_stats
)

print(bcmvn)


# ------------------------------------------------------------
# 6. Save pK results
# ------------------------------------------------------------

write.csv(
  bcmvn,
  file.path(
    doublet_dir,
    "DoubletFinder_pK_results.csv"
  ),
  row.names = FALSE
)

cat("\n========================================\n")
cat("Parameter sweep completed\n")
cat("========================================\n")