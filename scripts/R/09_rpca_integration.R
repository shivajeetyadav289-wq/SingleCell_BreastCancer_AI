# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 09_rpca_integration.R
#
# Purpose:
#   Integrate the 9 breast cancer patient samples using
#   Seurat v5 Reciprocal PCA (RPCA) integration.
#
# Input:
#   GSE228499_integration_prepared.rds
#
# Output:
#   GSE228499_integrated_rpca.rds
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(Seurat)


# ------------------------------------------------------------
# 2. Project paths
# ------------------------------------------------------------

project_dir <- "~/Projects/SingleCell_BreastCancer_AI"

input_file <- file.path(
  project_dir,
  "data",
  "processed",
  "GSE228499_integration_prepared.rds"
)

output_file <- file.path(
  project_dir,
  "data",
  "processed",
  "GSE228499_integrated_rpca.rds"
)


# ------------------------------------------------------------
# 3. Check input file
# ------------------------------------------------------------

if (!file.exists(input_file)) {

  stop(
    paste(
      "Input file not found:",
      input_file
    )
  )

}


# ------------------------------------------------------------
# 4. Load integration-ready object
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("LOADING INTEGRATION-READY DATASET\n")
cat("====================================================\n")

obj <- readRDS(
  input_file
)


cat(
  "Genes:",
  nrow(obj),
  "\n"
)

cat(
  "Cells:",
  ncol(obj),
  "\n"
)


# ------------------------------------------------------------
# 5. Check samples
# ------------------------------------------------------------

cat("\nCells per sample:\n")

print(
  table(obj$sample_id)
)


# ------------------------------------------------------------
# 6. Check RNA layers
# ------------------------------------------------------------

cat("\nRNA layers:\n")

print(
  Layers(obj[["RNA"]])
)


# ------------------------------------------------------------
# 7. Check PCA
# ------------------------------------------------------------

cat("\nAvailable dimensional reductions:\n")

print(
  Reductions(obj)
)


if (!"pca" %in% Reductions(obj)) {

  stop(
    "PCA reduction was not found."
  )

}


cat("\nPCA dimensions:\n")

print(
  dim(
    Embeddings(
      obj,
      reduction = "pca"
    )
  )
)


# ------------------------------------------------------------
# 8. Define integration dimensions
# ------------------------------------------------------------

integration_dims <- 1:30

cat("\n")
cat(
  "Integration dimensions:",
  paste(
    integration_dims,
    collapse = ", "
  ),
  "\n"
)


# ------------------------------------------------------------
# 9. Run RPCA integration
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("RUNNING RPCA INTEGRATION\n")
cat("====================================================\n")

obj <- IntegrateLayers(
  object = obj,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  dims = integration_dims,
  verbose = TRUE
)


# ------------------------------------------------------------
# 10. Check integrated reduction
# ------------------------------------------------------------

cat("\n")
cat("Available dimensional reductions after integration:\n")

print(
  Reductions(obj)
)


if (
  !"integrated.rpca" %in% Reductions(obj)
) {

  stop(
    "integrated.rpca reduction was not created."
  )

}


cat("\nIntegrated RPCA dimensions:\n")

print(
  dim(
    Embeddings(
      obj,
      reduction = "integrated.rpca"
    )
  )
)


# ------------------------------------------------------------
# 11. Save integrated object
# ------------------------------------------------------------

cat("\nSaving integrated object...\n")

saveRDS(
  obj,
  output_file
)


# ------------------------------------------------------------
# 12. Final summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("RPCA INTEGRATION COMPLETED\n")
cat("====================================================\n")

cat(
  "Genes:",
  nrow(obj),
  "\n"
)

cat(
  "Cells:",
  ncol(obj),
  "\n"
)

cat(
  "Samples:",
  length(
    unique(obj$sample_id)
  ),
  "\n"
)

cat(
  "Integration method: RPCA\n"
)

cat(
  "Dimensions used:",
  length(integration_dims),
  "\n"
)

cat("\nCells per sample:\n")

print(
  table(obj$sample_id)
)

cat("\nDimensional reductions:\n")

print(
  Reductions(obj)
)

cat(
  "\nSaved to:\n",
  output_file,
  "\n"
)

cat("\n")
cat("====================================================\n")
cat("READY FOR NEIGHBORS + UMAP + CLUSTERING\n")
cat("====================================================\n")