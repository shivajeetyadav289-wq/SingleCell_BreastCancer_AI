# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 05b_clean_singlet_metadata.R
#
# Purpose:
#   Remove unnecessary DoubletFinder intermediate metadata
#   while preserving important QC and biological metadata.
# ============================================================

library(Seurat)

project_dir <- "~/Projects/SingleCell_BreastCancer_AI"

input_file <- file.path(
  project_dir,
  "data",
  "processed",
  "GSE228499_singlets.rds"
)

output_file <- file.path(
  project_dir,
  "data",
  "processed",
  "GSE228499_singlets_clean.rds"
)


# ------------------------------------------------------------
# Load object
# ------------------------------------------------------------

singlets <- readRDS(input_file)

cat("\nOriginal dimensions:\n")
print(dim(singlets))


# ------------------------------------------------------------
# Identify DoubletFinder intermediate columns
# ------------------------------------------------------------

df_columns <- grep(
  "^(pANN_|DF\\.classifications_)",
  colnames(singlets@meta.data),
  value = TRUE
)

cat("\nDoubletFinder intermediate columns to remove:\n")
print(df_columns)


# ------------------------------------------------------------
# Remove unnecessary columns
# ------------------------------------------------------------

if (length(df_columns) > 0) {

  singlets@meta.data <- singlets@meta.data[
    ,
    !colnames(singlets@meta.data) %in% df_columns,
    drop = FALSE
  ]
}


# ------------------------------------------------------------
# Keep important metadata
# ------------------------------------------------------------

cat("\nRemaining metadata:\n")
print(colnames(singlets@meta.data))


# ------------------------------------------------------------
# Verify doublet status
# ------------------------------------------------------------

cat("\nDoublet status:\n")
print(table(singlets$doublet_status))


# ------------------------------------------------------------
# Verify samples
# ------------------------------------------------------------

cat("\nCells per sample:\n")
print(table(singlets$sample_id))


# ------------------------------------------------------------
# Verify biological metadata
# ------------------------------------------------------------

cat("\nCancer subtype:\n")
print(table(singlets$cancer_subtype))

cat("\nTreatment:\n")
print(table(singlets$treatment))


# ------------------------------------------------------------
# Save clean object
# ------------------------------------------------------------

saveRDS(
  singlets,
  output_file
)

cat("\nClean singlet object saved:\n")
cat(output_file, "\n")


# ------------------------------------------------------------
# Final dimensions
# ------------------------------------------------------------

cat("\nFinal dimensions:\n")
print(dim(singlets))

cat("\nMetadata cleanup completed.\n")