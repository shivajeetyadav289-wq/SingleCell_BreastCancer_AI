# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 03_qc_filtering.R
# Purpose: Filter low-quality cells after QC assessment
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(Seurat)


# ------------------------------------------------------------
# 2. Define directories
# ------------------------------------------------------------

project_dir <- "~/Projects/SingleCell_BreastCancer_AI"

processed_dir <- file.path(
  project_dir,
  "data",
  "processed"
)

qc_dir <- file.path(
  project_dir,
  "results",
  "qc"
)


# ------------------------------------------------------------
# 3. Load raw merged object
# ------------------------------------------------------------

combined <- readRDS(
  file.path(
    processed_dir,
    "GSE228499_merged_raw.rds"
  )
)


# ------------------------------------------------------------
# 4. Calculate mitochondrial percentage
# ------------------------------------------------------------

combined[["percent.mt"]] <- PercentageFeatureSet(
  combined,
  pattern = "^MT-"
)


# ------------------------------------------------------------
# 5. Define QC thresholds
# ------------------------------------------------------------

min_features <- 200

max_mito <- 20


# ------------------------------------------------------------
# 6. Record cells before filtering
# ------------------------------------------------------------

cells_before <- ncol(combined)


# ------------------------------------------------------------
# 7. Filter cells
# ------------------------------------------------------------

filtered <- subset(
  combined,
  subset =
    nFeature_RNA >= min_features &
    percent.mt < max_mito
)


# ------------------------------------------------------------
# 8. Record cells after filtering
# ------------------------------------------------------------

cells_after <- ncol(filtered)

cells_removed <- cells_before - cells_after

retention_percent <- (
  cells_after / cells_before
) * 100


# ------------------------------------------------------------
# 9. Print filtering summary
# ------------------------------------------------------------

cat("\n========================================\n")
cat("QC Filtering Summary\n")
cat("========================================\n")

cat(
  "\nMinimum genes per cell:",
  min_features,
  "\n"
)

cat(
  "Maximum mitochondrial percentage:",
  max_mito,
  "%\n"
)

cat(
  "\nCells before filtering:",
  cells_before,
  "\n"
)

cat(
  "Cells after filtering:",
  cells_after,
  "\n"
)

cat(
  "Cells removed:",
  cells_removed,
  "\n"
)

cat(
  "Retention percentage:",
  round(retention_percent, 2),
  "%\n"
)


# ------------------------------------------------------------
# 10. Cell counts by sample after filtering
# ------------------------------------------------------------

cells_after_sample <- as.data.frame(
  table(filtered$sample_id)
)

colnames(cells_after_sample) <- c(
  "sample_id",
  "cells_after_qc"
)

print(cells_after_sample)


# ------------------------------------------------------------
# 11. Save cell counts after QC
# ------------------------------------------------------------

write.csv(
  cells_after_sample,
  file.path(
    qc_dir,
    "cells_per_sample_after_qc.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 12. Save filtered Seurat object
# ------------------------------------------------------------

saveRDS(
  filtered,
  file.path(
    processed_dir,
    "GSE228499_filtered.rds"
  )
)


# ------------------------------------------------------------
# 13. Final confirmation
# ------------------------------------------------------------

cat("\n========================================\n")
cat("Filtered Seurat object saved successfully\n")
cat("========================================\n")

cat(
  "\nFile:",
  file.path(
    processed_dir,
    "GSE228499_filtered.rds"
  ),
  "\n"
)

cat(
  "\nFinal dimensions:",
  nrow(filtered),
  "genes x",
  ncol(filtered),
  "cells\n"
)