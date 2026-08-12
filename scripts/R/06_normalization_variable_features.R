# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 08_prepare_integration.R
#
# Purpose:
#   Prepare the clean singlet dataset for Seurat v5
#   sample-aware RPCA integration.
#
# Strategy:
#   1. Load clean singlet object
#   2. Join old Seurat v5 layers
#   3. Extract ONLY the raw counts layer
#   4. Preserve metadata
#   5. Create a completely fresh Seurat object
#   6. Split RNA counts by sample
#   7. Normalize
#   8. Find highly variable genes
#   9. Scale
#  10. Run PCA
#  11. Save integration-ready object
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
  "GSE228499_singlets_clean.rds"
)

output_file <- file.path(
  project_dir,
  "data",
  "processed",
  "GSE228499_integration_prepared.rds"
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
# 4. Load clean singlet object
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("LOADING CLEAN SINGLET DATASET\n")
cat("====================================================\n")

old_obj <- readRDS(
  input_file
)

cat(
  "Genes:",
  nrow(old_obj),
  "\n"
)

cat(
  "Cells:",
  ncol(old_obj),
  "\n"
)


# ------------------------------------------------------------
# 5. Check sample metadata
# ------------------------------------------------------------

if (!"sample_id" %in% colnames(old_obj@meta.data)) {
  stop(
    "sample_id column is missing from metadata."
  )
}

cat("\nCells per sample:\n")

print(
  table(old_obj$sample_id)
)


# ------------------------------------------------------------
# 6. Inspect existing RNA layers
# ------------------------------------------------------------

cat("\n")
cat("RNA layers BEFORE JoinLayers:\n")

print(
  Layers(old_obj[["RNA"]])
)


# ------------------------------------------------------------
# 7. Join existing layers
# ------------------------------------------------------------

cat("\n")
cat("Joining existing RNA layers...\n")

old_obj <- JoinLayers(
  object = old_obj,
  assay = "RNA"
)


cat("\n")
cat("RNA layers AFTER JoinLayers:\n")

print(
  Layers(old_obj[["RNA"]])
)


# ------------------------------------------------------------
# IMPORTANT:
#
# After JoinLayers(), it is NORMAL to have:
#
# counts
# data
# scale.data
#
# We only need the counts layer to build a clean object.
# ------------------------------------------------------------


# ------------------------------------------------------------
# 8. Extract raw counts
# ------------------------------------------------------------

cat("\n")
cat("Extracting counts layer...\n")

counts_matrix <- LayerData(
  object = old_obj,
  assay = "RNA",
  layer = "counts"
)


cat(
  "Counts dimensions:",
  nrow(counts_matrix),
  "genes x",
  ncol(counts_matrix),
  "cells\n"
)


# ------------------------------------------------------------
# 9. Extract metadata
# ------------------------------------------------------------

cat("\n")
cat("Extracting cell metadata...\n")

metadata <- old_obj@meta.data


cat(
  "Metadata dimensions:",
  nrow(metadata),
  "cells x",
  ncol(metadata),
  "columns\n"
)


# ------------------------------------------------------------
# 10. Verify cell names
# ------------------------------------------------------------

if (!all(colnames(counts_matrix) %in% rownames(metadata))) {

  stop(
    "Some count-matrix cell names are missing from metadata."
  )

}


# ------------------------------------------------------------
# 11. Reorder metadata to match counts
# ------------------------------------------------------------

metadata <- metadata[
  colnames(counts_matrix),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 12. Verify exact cell order
# ------------------------------------------------------------

if (
  !identical(
    colnames(counts_matrix),
    rownames(metadata)
  )
) {

  stop(
    "Counts columns and metadata rows are not in identical order."
  )

}


cat(
  "\nCell names successfully matched.\n"
)


# ------------------------------------------------------------
# 13. Create fresh Seurat object
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("CREATING FRESH SEURAT OBJECT\n")
cat("====================================================\n")

obj <- CreateSeuratObject(
  counts = counts_matrix,
  meta.data = metadata,
  project = "GSE228499"
)


# ------------------------------------------------------------
# 14. Check fresh object
# ------------------------------------------------------------

cat("\nFresh object dimensions:\n")

print(
  dim(obj)
)


cat("\nFresh RNA layers:\n")

print(
  Layers(obj[["RNA"]])
)


# ------------------------------------------------------------
# 15. Check sample information again
# ------------------------------------------------------------

cat("\nCells per sample in fresh object:\n")

print(
  table(obj$sample_id)
)


# ------------------------------------------------------------
# 16. Verify number of samples
# ------------------------------------------------------------

number_of_samples <- length(
  unique(obj$sample_id)
)

if (number_of_samples != 9) {

  stop(
    paste(
      "Expected 9 samples but found",
      number_of_samples
    )
  )

}


# ------------------------------------------------------------
# 17. Split counts layer by sample
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("SPLITTING COUNTS BY SAMPLE\n")
cat("====================================================\n")

obj[["RNA"]] <- split(
  obj[["RNA"]],
  f = obj$sample_id
)


# ------------------------------------------------------------
# 18. Display split layers
# ------------------------------------------------------------

cat("\nRNA layers AFTER SPLITTING:\n")

print(
  Layers(obj[["RNA"]])
)


# ------------------------------------------------------------
# 19. Check that split occurred
# ------------------------------------------------------------

rna_layers <- Layers(
  obj[["RNA"]]
)

count_layers <- grep(
  "^counts",
  rna_layers,
  value = TRUE
)

if (length(count_layers) != 9) {

  stop(
    paste(
      "Expected 9 count layers after splitting, found",
      length(count_layers)
    )
  )

}


cat(
  "\nSuccessfully created",
  length(count_layers),
  "sample-specific count layers.\n"
)


# ------------------------------------------------------------
# 20. Normalize data
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("NORMALIZATION\n")
cat("====================================================\n")

obj <- NormalizeData(
  object = obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)


# ------------------------------------------------------------
# 21. Find highly variable genes
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("HIGHLY VARIABLE GENE SELECTION\n")
cat("====================================================\n")

obj <- FindVariableFeatures(
  object = obj,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = TRUE
)


cat(
  "\nNumber of variable genes:",
  length(
    VariableFeatures(obj)
  ),
  "\n"
)


# ------------------------------------------------------------
# 22. Scale data
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("SCALING DATA\n")
cat("====================================================\n")

obj <- ScaleData(
  object = obj,
  features = VariableFeatures(obj),
  verbose = TRUE
)


# ------------------------------------------------------------
# 23. Run PCA
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("PCA\n")
cat("====================================================\n")

obj <- RunPCA(
  object = obj,
  features = VariableFeatures(obj),
  npcs = 30,
  verbose = TRUE
)


# ------------------------------------------------------------
# 24. Verify PCA
# ------------------------------------------------------------

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
# 25. Save integration-ready object
# ------------------------------------------------------------

cat("\n")
cat("Saving integration-ready object...\n")

saveRDS(
  obj,
  output_file
)


# ------------------------------------------------------------
# 26. Final summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("INTEGRATION PREPARATION COMPLETED\n")
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
  "Highly variable genes:",
  length(
    VariableFeatures(obj)
  ),
  "\n"
)

cat(
  "PCA dimensions:",
  ncol(
    Embeddings(
      obj,
      reduction = "pca"
    )
  ),
  "\n"
)

cat("\nCells per sample:\n")

print(
  table(obj$sample_id)
)

cat("\nFinal RNA layers:\n")

print(
  Layers(obj[["RNA"]])
)

cat(
  "\nSaved to:\n",
  output_file,
  "\n"
)

cat("\n")
cat("====================================================\n")
cat("READY FOR RPCA INTEGRATION\n")
cat("====================================================\n")