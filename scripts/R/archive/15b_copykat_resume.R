# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 15b_copykat_resume.R
#
# Purpose:
# Resume CopyKAT CNV inference from the saved raw count matrix.
# ============================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(copykat)
})

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

project_dir <- path.expand(
  "~/Projects/SingleCell_BreastCancer_AI"
)

result_dir <- file.path(
  project_dir,
  "results",
  "cnv"
)

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

count_file <- file.path(
  result_dir,
  "GSE228499_tumor_candidate_counts.rds"
)

# ------------------------------------------------------------
# 2. Load checkpoint
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("COPYKAT RESUME\n")
cat("====================================================\n")

cat("Loading saved count matrix...\n")

counts_matrix <- readRDS(
  count_file
)

cat(
  "Genes:",
  nrow(counts_matrix),
  "\n"
)

cat(
  "Cells:",
  ncol(counts_matrix),
  "\n"
)

cat(
  "Matrix class:",
  class(counts_matrix),
  "\n"
)

# ------------------------------------------------------------
# 3. Validate checkpoint
# ------------------------------------------------------------

if (
  nrow(counts_matrix) != 36626 ||
  ncol(counts_matrix) != 17574
) {

  stop(
    "Count matrix dimensions are not as expected."
  )
}

if (
  !inherits(
    counts_matrix,
    "dgCMatrix"
  )
) {

  counts_matrix <- as(
    counts_matrix,
    "dgCMatrix"
  )
}

cat("\nCheckpoint validation: PASSED\n")

# ------------------------------------------------------------
# 4. Run CopyKAT
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("RUNNING COPYKAT\n")
cat("====================================================\n")

cat(
  "Input:",
  nrow(counts_matrix),
  "genes x",
  ncol(counts_matrix),
  "cells\n"
)

cat("CPU cores: 4\n")
cat("This may take a long time.\n")
cat("Do not start another R analysis simultaneously.\n")

copykat_result <- copykat(
  rawmat = counts_matrix,
  id.type = "S",
  cell.group = "",
  cell.line = "no",
  ngene.chr = 5,
  min.gene.per.cell = 200,
  LOW.DR = 0.05,
  UP.DR = 0.1,
  win.size = 25,
  norm.cell.names = "",
  KS.cut = 0.1,
  sam.name = "GSE228499",
  distance = "euclidean",
  test.emd = "FALSE",
  output.seg = "FALSE",
  plot.genes = "TRUE",
  genome = "hg20",
  n.cores = 2
)

# ------------------------------------------------------------
# 5. Save immediately
# ------------------------------------------------------------

cat("\n")
cat("CopyKAT finished.\n")
cat("Saving result immediately...\n")

copykat_file <- file.path(
  result_dir,
  "GSE228499_CopyKAT_result.rds"
)

saveRDS(
  copykat_result,
  copykat_file
)

cat(
  "Saved:\n",
  copykat_file,
  "\n"
)

# ------------------------------------------------------------
# 6. Prediction table
# ------------------------------------------------------------

if (
  "prediction" %in%
  names(copykat_result)
) {

  prediction <- copykat_result$prediction

  prediction_file <- file.path(
    result_dir,
    "CopyKAT_prediction.csv"
  )

  write.csv(
    prediction,
    prediction_file,
    row.names = FALSE
  )

  cat(
    "\nPrediction table saved:\n",
    prediction_file,
    "\n"
  )

  cat("\nCopyKAT predictions:\n")

  print(
    table(
      prediction$copykat.pred,
      useNA = "ifany"
    )
  )
}

# ------------------------------------------------------------
# 7. Save CNV matrix
# ------------------------------------------------------------

if (
  "CNAmat" %in%
  names(copykat_result)
) {

  cnv_file <- file.path(
    result_dir,
    "CopyKAT_CNAmat.rds"
  )

  saveRDS(
    copykat_result$CNAmat,
    cnv_file
  )

  cat(
    "\nCNV matrix saved:\n",
    cnv_file,
    "\n"
  )
}

# ------------------------------------------------------------
# 8. Final summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("COPYKAT COMPLETED\n")
cat("====================================================\n")

cat(
  "Input genes:",
  nrow(counts_matrix),
  "\n"
)

cat(
  "Input cells:",
  ncol(counts_matrix),
  "\n"
)

cat(
  "Result:",
  copykat_file,
  "\n"
)

cat("\nNEXT STEP:\n")
cat(
  "CNV-based malignant-cell validation.\n"
)

cat("====================================================\n")