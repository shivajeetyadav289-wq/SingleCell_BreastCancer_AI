# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 15_cnv_inference.R
#
# Purpose:
# CNV inference using CopyKAT to identify putative malignant
# epithelial cells.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(copykat)
})

# ------------------------------------------------------------
# 1. Paths
# ------------------------------------------------------------

project_dir <- path.expand(
  "~/Projects/SingleCell_BreastCancer_AI"
)

input_file <- file.path(
  project_dir,
  "data",
  "processed",
  "GSE228499_tumor_candidate_epithelial.rds"
)

result_dir <- file.path(
  project_dir,
  "results",
  "cnv"
)

figure_dir <- file.path(
  project_dir,
  "figures",
  "cnv"
)

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. Load candidate epithelial cells
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("COPYKAT CNV ANALYSIS\n")
cat("====================================================\n")

cat("Loading tumor-candidate epithelial object...\n")

obj <- readRDS(input_file)

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

cat("\nCells per sample:\n")
print(table(obj$sample_id))

# ------------------------------------------------------------
# 3. Identify raw count layers
# ------------------------------------------------------------

rna_layers <- Layers(obj[["RNA"]])

count_layers <- rna_layers[
  grepl("^counts\\.", rna_layers)
]

cat("\nRaw count layers:\n")
print(count_layers)

if (length(count_layers) == 0) {
  stop("No sample-specific count layers found.")
}

# ------------------------------------------------------------
# 4. Extract raw counts
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("EXTRACTING RAW COUNTS\n")
cat("====================================================\n")

count_list <- lapply(
  count_layers,
  function(layer_name) {

    cat(
      "Reading:",
      layer_name,
      "\n"
    )

    LayerData(
      obj,
      assay = "RNA",
      layer = layer_name
    )
  }
)

# ------------------------------------------------------------
# 5. Find common genes
# ------------------------------------------------------------

common_genes <- Reduce(
  intersect,
  lapply(
    count_list,
    rownames
  )
)

cat(
  "\nCommon genes:",
  length(common_genes),
  "\n"
)

count_list <- lapply(
  count_list,
  function(x) {
    x[
      common_genes,
      ,
      drop = FALSE
    ]
  }
)

# ------------------------------------------------------------
# 6. Combine count matrices
# ------------------------------------------------------------

cat("\nCombining count matrices...\n")

counts_matrix <- do.call(
  cbind,
  count_list
)

counts_matrix <- as(
  counts_matrix,
  "dgCMatrix"
)

cat(
  "Combined matrix:",
  nrow(counts_matrix),
  "genes x",
  ncol(counts_matrix),
  "cells\n"
)

# ------------------------------------------------------------
# 7. Verify cell names
# ------------------------------------------------------------

if (!all(
  colnames(counts_matrix) %in%
    colnames(obj)
)) {

  stop(
    "Count-matrix cell names do not match Seurat object."
  )
}

# Reorder metadata to match count matrix

metadata <- obj@meta.data[
  colnames(counts_matrix),
  ,
  drop = FALSE
]

if (!all(
  rownames(metadata) ==
    colnames(counts_matrix)
)) {

  stop(
    "Metadata and count matrix are not aligned."
  )
}

# ------------------------------------------------------------
# 8. Save count matrix checkpoint
# ------------------------------------------------------------

count_file <- file.path(
  result_dir,
  "GSE228499_tumor_candidate_counts.rds"
)

saveRDS(
  counts_matrix,
  count_file
)

cat(
  "\nCount matrix checkpoint saved:\n",
  count_file,
  "\n"
)

# ------------------------------------------------------------
# 9. Save metadata
# ------------------------------------------------------------

metadata_file <- file.path(
  result_dir,
  "GSE228499_tumor_candidate_metadata.csv"
)

write.csv(
  metadata,
  metadata_file,
  row.names = TRUE
)

cat(
  "Metadata saved:\n",
  metadata_file,
  "\n"
)

# ------------------------------------------------------------
# 10. Run CopyKAT
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("RUNNING COPYKAT\n")
cat("====================================================\n")

cat(
  "Cells entering CopyKAT:",
  ncol(counts_matrix),
  "\n"
)

cat(
  "Genes entering CopyKAT:",
  nrow(counts_matrix),
  "\n"
)

cat(
  "Using 4 CPU cores.\n"
)

cat(
  "This step may take a substantial amount of time.\n"
)

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
  n.cores = 4
)

# ------------------------------------------------------------
# 11. Save CopyKAT object
# ------------------------------------------------------------

copykat_file <- file.path(
  result_dir,
  "GSE228499_CopyKAT_result.rds"
)

saveRDS(
  copykat_result,
  copykat_file
)

cat(
  "\nCopyKAT result saved:\n",
  copykat_file,
  "\n"
)

# ------------------------------------------------------------
# 12. Inspect CopyKAT output
# ------------------------------------------------------------

cat("\n")
cat("CopyKAT result components:\n")

print(
  names(copykat_result)
)

# ------------------------------------------------------------
# 13. Save prediction table
# ------------------------------------------------------------

if ("prediction" %in% names(copykat_result)) {

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

  cat("\nCopyKAT prediction summary:\n")

  print(
    table(
      prediction$copykat.pred,
      useNA = "ifany"
    )
  )
}

# ------------------------------------------------------------
# 14. Save CNV matrix if available
# ------------------------------------------------------------

if ("CNAmat" %in% names(copykat_result)) {

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
# 15. Final summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("COPYKAT ANALYSIS COMPLETED\n")
cat("====================================================\n")

cat(
  "Input cells:",
  ncol(counts_matrix),
  "\n"
)

cat(
  "Input genes:",
  nrow(counts_matrix),
  "\n"
)

cat(
  "Results directory:\n",
  result_dir,
  "\n"
)

cat("\n")
cat("NEXT STEP:\n")
cat(
  "Validate CopyKAT diploid/aneuploid predictions and\n"
)

cat(
  "identify CNV-supported malignant epithelial cells.\n"
)

cat("====================================================\n")