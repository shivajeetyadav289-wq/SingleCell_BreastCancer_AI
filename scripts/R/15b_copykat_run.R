# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 15b_copykat_run.R
#
# Purpose:
# Run CopyKAT CNV inference on the validated 4,000-cell
# representative subset of tumor-candidate epithelial cells.
#
# Full tumor candidates: 17,574
# CopyKAT subset:         4,000
#
# Hardware:
# Physical RAM: ~7.5 GB
# CopyKAT cores: 2
#
# Important:
# The original 17,574-cell matrix is NEVER converted to dense.
# Only the already filtered 4,000-cell matrix is converted.
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

suppressPackageStartupMessages({

  library(copykat)
  library(Matrix)

})


# ------------------------------------------------------------
# 2. Paths
# ------------------------------------------------------------

project_dir <- path.expand(
  "~/Projects/SingleCell_BreastCancer_AI"
)

subset_dir <- file.path(
  project_dir,
  "results",
  "cnv",
  "subset"
)

output_dir <- file.path(
  project_dir,
  "results",
  "cnv",
  "copykat"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. Input files
# ------------------------------------------------------------

count_file <- file.path(
  subset_dir,
  "GSE228499_CopyKAT_subset_counts.rds"
)

metadata_file <- file.path(
  subset_dir,
  "GSE228499_CopyKAT_subset_metadata.csv"
)


# ------------------------------------------------------------
# 4. Settings
# ------------------------------------------------------------

N_CORES <- 2

RANDOM_SEED <- 228499

SAMPLE_NAME <- "GSE228499_CopyKAT_4000"


# ------------------------------------------------------------
# 5. Header
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("COPYKAT CNV INFERENCE\n")
cat("====================================================\n")

cat(
  "CopyKAT version:",
  as.character(
    packageVersion("copykat")
  ),
  "\n"
)

cat(
  "Matrix checkpoint:",
  count_file,
  "\n"
)

cat(
  "CPU cores:",
  N_CORES,
  "\n"
)


# ------------------------------------------------------------
# 6. Check files
# ------------------------------------------------------------

if (!file.exists(count_file)) {

  stop(
    "Count matrix not found:\n",
    count_file
  )

}

if (!file.exists(metadata_file)) {

  stop(
    "Metadata not found:\n",
    metadata_file
  )

}


# ------------------------------------------------------------
# 7. Load subset
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("LOADING 4,000-CELL CHECKPOINT\n")
cat("====================================================\n")

counts <- readRDS(
  count_file
)

cat(
  "Genes:",
  nrow(counts),
  "\n"
)

cat(
  "Cells:",
  ncol(counts),
  "\n"
)

cat(
  "Class:",
  class(counts)[1],
  "\n"
)


# ------------------------------------------------------------
# 8. Validate sparse matrix
# ------------------------------------------------------------

if (!inherits(counts, "dgCMatrix")) {

  stop(
    "Expected dgCMatrix."
  )

}


if (
  ncol(counts) != 4000
) {

  stop(
    "Expected 4,000 cells but found ",
    ncol(counts)
  )

}


# ------------------------------------------------------------
# 9. Load metadata
# ------------------------------------------------------------

metadata <- read.csv(
  metadata_file,
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

metadata$cell_id <- rownames(
  metadata
)


# ------------------------------------------------------------
# 10. Validate metadata
# ------------------------------------------------------------

if (
  nrow(metadata) !=
    ncol(counts)
) {

  stop(
    "Metadata and count matrix dimensions do not match."
  )

}


if (
  !identical(
    metadata$cell_id,
    colnames(counts)
  )
) {

  stop(
    "Metadata cell IDs do not match count matrix."
  )

}


cat(
  "Metadata validation: PASSED\n"
)


# ------------------------------------------------------------
# 11. Sample information
# ------------------------------------------------------------

cat("\n")
cat("Cells per sample:\n")

print(
  table(
    metadata$sample_id
  )
)


# ------------------------------------------------------------
# 12. LOW.DR filtering
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("COPYKAT GENE FILTERING\n")
cat("====================================================\n")

detection_rate <- Matrix::rowSums(
  counts > 0
) /
  ncol(counts)

filtered_counts <- counts[
  detection_rate > 0.05,
  ,
  drop = FALSE
]

cat(
  "Genes before filtering:",
  nrow(counts),
  "\n"
)

cat(
  "Genes after LOW.DR filtering:",
  nrow(filtered_counts),
  "\n"
)

cat(
  "Cells:",
  ncol(filtered_counts),
  "\n"
)

cat(
  "Class:",
  class(filtered_counts)[1],
  "\n"
)


# ------------------------------------------------------------
# 13. Save filtered sparse checkpoint
# ------------------------------------------------------------

filtered_file <- file.path(
  output_dir,
  "GSE228499_CopyKAT_filtered_matrix.rds"
)

saveRDS(
  filtered_counts,
  filtered_file
)

cat(
  "\nFiltered sparse checkpoint saved:\n",
  filtered_file,
  "\n"
)


# ------------------------------------------------------------
# 14. Memory cleanup
# ------------------------------------------------------------

gc()


# ------------------------------------------------------------
# 15. Convert ONLY the 4,000-cell filtered matrix
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("DENSE MATRIX CONVERSION\n")
cat("====================================================\n")

cat(
  "IMPORTANT: Only the filtered 4,000-cell matrix is being\n",
  "converted to dense format.\n",
  sep = ""
)

cat(
  "The original 17,574-cell matrix is NOT converted.\n"
)

cat(
  "Sparse dimensions:",
  nrow(filtered_counts),
  "x",
  ncol(filtered_counts),
  "\n"
)


# Estimate memory

estimated_mb <-
  nrow(filtered_counts) *
  ncol(filtered_counts) *
  8 /
  1024^2

cat(
  "Estimated dense matrix size:",
  round(
    estimated_mb,
    1
  ),
  "MB\n"
)


# ------------------------------------------------------------
# 16. Dense conversion
# ------------------------------------------------------------

dense_counts <- as.matrix(
  filtered_counts
)

cat(
  "Dense matrix conversion: PASSED\n"
)

cat(
  "Dense dimensions:",
  nrow(dense_counts),
  "x",
  ncol(dense_counts),
  "\n"
)

cat(
  "Dense class:",
  class(dense_counts),
  "\n"
)


# ------------------------------------------------------------
# 17. Remove sparse matrix
# ------------------------------------------------------------

rm(
  filtered_counts,
  counts,
  detection_rate
)

gc()


# ------------------------------------------------------------
# 18. Verify dense matrix
# ------------------------------------------------------------

if (
  !is.matrix(dense_counts)
) {

  stop(
    "Dense conversion failed."
  )

}


if (
  anyNA(dense_counts)
) {

  stop(
    "NA values detected in dense matrix."
  )

}


cat(
  "Dense matrix validation: PASSED\n"
)


# ------------------------------------------------------------
# 19. Save dense checkpoint
# ------------------------------------------------------------

dense_checkpoint <- file.path(
  output_dir,
  "GSE228499_CopyKAT_dense_4000.rds"
)

saveRDS(
  dense_counts,
  dense_checkpoint
)

cat(
  "Dense checkpoint saved:\n",
  dense_checkpoint,
  "\n"
)


# ------------------------------------------------------------
# 20. Test CopyKAT annotation
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("TESTING COPYKAT ANNOTATION\n")
cat("====================================================\n")

annotation_test <- tryCatch(

  {

    copykat:::annotateGenes.hg20(
      mat = dense_counts,
      ID.type = "S"
    )

  },

  error = function(e) {

    stop(
      "CopyKAT annotation failed:\n",
      conditionMessage(e)
    )

  }

)


cat(
  "CopyKAT annotation: PASSED\n"
)

cat(
  "Annotated genes:",
  nrow(annotation_test),
  "\n"
)


# ------------------------------------------------------------
# 21. Save annotation checkpoint
# ------------------------------------------------------------

annotation_file <- file.path(
  output_dir,
  "GSE228499_CopyKAT_annotation.rds"
)

saveRDS(
  annotation_test,
  annotation_file
)

cat(
  "Annotation checkpoint saved.\n"
)


# ------------------------------------------------------------
# 22. Remove annotation object
# ------------------------------------------------------------

rm(
  annotation_test
)

gc()


# ------------------------------------------------------------
# 23. Run CopyKAT
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("RUNNING COPYKAT\n")
cat("====================================================\n")

cat(
  "Input:",
  nrow(dense_counts),
  "genes x",
  ncol(dense_counts),
  "cells\n"
)

cat(
  "Input class:",
  class(dense_counts),
  "\n"
)

cat(
  "CPU cores:",
  N_CORES,
  "\n"
)

cat(
  "Genome: hg20\n"
)

cat(
  "KS.cut: 0.1\n"
)

cat(
  "This may take a long time.\n"
)

cat(
  "Do NOT start another R analysis.\n"
)


set.seed(
  RANDOM_SEED
)


# ------------------------------------------------------------
# 24. CopyKAT execution
# ------------------------------------------------------------

copykat_result <- tryCatch(

  {

    copykat(

      rawmat = dense_counts,

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

      sam.name = SAMPLE_NAME,

      distance = "euclidean",

      test.emd = "FALSE",

      output.seg = "FALSE",

      plot.genes = "TRUE",

      genome = "hg20",

      n.cores = N_CORES

    )

  },

  error = function(e) {

    cat("\n")
    cat("====================================================\n")
    cat("COPYKAT ERROR\n")
    cat("====================================================\n")

    cat(
      conditionMessage(e),
      "\n"
    )

    cat("\n")
    cat(
      "Dense checkpoint preserved:\n"
    )

    cat(
      dense_checkpoint,
      "\n"
    )

    stop(
      e
    )

  }

)


# ------------------------------------------------------------
# 25. Save complete result
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("SAVING COPYKAT RESULT\n")
cat("====================================================\n")

result_file <- file.path(
  output_dir,
  "GSE228499_CopyKAT_result.rds"
)

saveRDS(
  copykat_result,
  result_file
)

cat(
  "Saved:",
  result_file,
  "\n"
)


# ------------------------------------------------------------
# 26. Inspect result
# ------------------------------------------------------------

cat("\n")
cat("CopyKAT result components:\n")

print(
  names(
    copykat_result
  )
)


# ------------------------------------------------------------
# 27. Save result components
# ------------------------------------------------------------

if (
  is.list(copykat_result)
) {

  for (
    nm in names(copykat_result)
  ) {

    component <- copykat_result[[nm]]

    if (
      !is.null(component)
    ) {

      safe_name <- gsub(
        "[^A-Za-z0-9_.-]",
        "_",
        nm
      )

      component_file <- file.path(
        output_dir,
        paste0(
          "GSE228499_CopyKAT_",
          safe_name,
          ".rds"
        )
      )

      saveRDS(
        component,
        component_file
      )

      cat(
        "Saved component:",
        component_file,
        "\n"
      )

    }

  }

}


# ------------------------------------------------------------
# 28. Final summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("COPYKAT CNV ANALYSIS COMPLETED\n")
cat("====================================================\n")

cat(
  "Input cells:",
  ncol(dense_counts),
  "\n"
)

cat(
  "Filtered genes:",
  nrow(dense_counts),
  "\n"
)

cat(
  "Samples:",
  length(
    unique(
      metadata$sample_id
    )
  ),
  "\n"
)

cat(
  "Result:",
  result_file,
  "\n"
)

cat("\n")
cat("Next step:\n")
cat(
  "Inspect CopyKAT CNV predictions and determine\n",
  "diploid/aneuploid populations.\n",
  sep = ""
)

cat("\n")
cat("====================================================\n")
cat("DO NOT YET CALL CELLS MALIGNANT\n")
cat("====================================================\n")


# ------------------------------------------------------------
# 29. Cleanup
# ------------------------------------------------------------

rm(
  dense_counts,
  copykat_result
)

gc()