# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 15_copykat_subset.R
#
# Purpose:
# Prepare a representative, sample-stratified subset of
# tumor-candidate epithelial cells for CopyKAT CNV analysis.
#
# Full dataset:        17,574 cells
# Target subset:        4,000 cells
#
# Hardware:
# Physical RAM: ~7.5 GB
# CPU for later CopyKAT: 2 cores
#
# IMPORTANT:
# This script ONLY prepares the subset.
# CopyKAT is NOT run here.
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(Matrix)
})


# ------------------------------------------------------------
# 2. Project paths
# ------------------------------------------------------------

project_dir <- path.expand(
  "~/Projects/SingleCell_BreastCancer_AI"
)

cnv_dir <- file.path(
  project_dir,
  "results",
  "cnv"
)

subset_dir <- file.path(
  cnv_dir,
  "subset"
)

dir.create(
  subset_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. Input files
# ------------------------------------------------------------

count_file <- file.path(
  cnv_dir,
  "GSE228499_tumor_candidate_counts.rds"
)

metadata_file <- file.path(
  cnv_dir,
  "GSE228499_tumor_candidate_metadata.csv"
)


# ------------------------------------------------------------
# 4. Settings
# ------------------------------------------------------------

TARGET_CELLS <- 4000

RANDOM_SEED <- 228499


# ------------------------------------------------------------
# 5. Header
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("COPYKAT STRATIFIED SUBSET PREPARATION\n")
cat("====================================================\n")

cat(
  "Target cells:",
  TARGET_CELLS,
  "\n"
)

cat(
  "Random seed:",
  RANDOM_SEED,
  "\n"
)


# ------------------------------------------------------------
# 6. Check input files
# ------------------------------------------------------------

if (!file.exists(count_file)) {

  stop(
    "Count matrix not found:\n",
    count_file
  )

}

if (!file.exists(metadata_file)) {

  stop(
    "Metadata file not found:\n",
    metadata_file
  )

}


# ------------------------------------------------------------
# 7. Load count matrix
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("LOADING COUNT MATRIX\n")
cat("====================================================\n")

counts <- readRDS(
  count_file
)

cat(
  "Dimensions: ",
  nrow(counts),
  " genes x ",
  ncol(counts),
  " cells\n",
  sep = ""
)

cat(
  "Matrix class: ",
  class(counts)[1],
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 8. Ensure sparse matrix
# ------------------------------------------------------------

if (!inherits(counts, "dgCMatrix")) {

  cat(
    "Converting count matrix to dgCMatrix...\n"
  )

  counts <- as(
    counts,
    "dgCMatrix"
  )

}


# ------------------------------------------------------------
# 9. Load metadata
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("LOADING METADATA\n")
cat("====================================================\n")

# IMPORTANT:
# The CSV stores cell barcodes as row names.
# Therefore row.names = 1 is required.

metadata <- read.csv(
  metadata_file,
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Store row names explicitly as cell IDs
metadata$cell_id <- rownames(
  metadata
)

cat(
  "Metadata rows:",
  nrow(metadata),
  "\n"
)

cat(
  "Metadata columns:",
  ncol(metadata),
  "\n"
)

cat(
  "Cell ID source: metadata row names\n"
)


# ------------------------------------------------------------
# 10. Show metadata columns
# ------------------------------------------------------------

cat("\nMetadata columns:\n")

print(
  colnames(metadata)
)


# ------------------------------------------------------------
# 11. Validate sample_id
# ------------------------------------------------------------

if (!"sample_id" %in% colnames(metadata)) {

  stop(
    "sample_id column is missing from metadata."
  )

}


# ------------------------------------------------------------
# 12. Validate cell IDs
# ------------------------------------------------------------

metadata_cells <- metadata$cell_id

cat("\n")
cat("====================================================\n")
cat("VALIDATING CELL IDENTIFIERS\n")
cat("====================================================\n")

cat(
  "Count matrix cells:",
  ncol(counts),
  "\n"
)

cat(
  "Metadata cells:",
  length(metadata_cells),
  "\n"
)


# Check duplicates

if (
  anyDuplicated(
    metadata_cells
  )
) {

  stop(
    "Duplicate cell IDs found in metadata."
  )

}


# Check number of cells

if (
  length(metadata_cells) !=
    ncol(counts)
) {

  stop(
    "Number of metadata cells does not match count matrix."
  )

}


# Check that all count-matrix cells exist in metadata

missing_metadata <- setdiff(
  colnames(counts),
  metadata_cells
)

if (
  length(missing_metadata) > 0
) {

  cat(
    "Missing metadata cells:",
    length(missing_metadata),
    "\n"
  )

  print(
    head(
      missing_metadata,
      10
    )
  )

  stop(
    "Some count-matrix cell IDs are missing from metadata."
  )

}


# Check that all metadata cells exist in count matrix

missing_counts <- setdiff(
  metadata_cells,
  colnames(counts)
)

if (
  length(missing_counts) > 0
) {

  cat(
    "Missing count-matrix cells:",
    length(missing_counts),
    "\n"
  )

  print(
    head(
      missing_counts,
      10
    )
  )

  stop(
    "Some metadata cell IDs are missing from count matrix."
  )

}


# ------------------------------------------------------------
# 13. Reorder metadata to count matrix order
# ------------------------------------------------------------

metadata <- metadata[
  match(
    colnames(counts),
    metadata$cell_id
  ),
  ,
  drop = FALSE
]


# Final exact check

if (
  !identical(
    metadata$cell_id,
    colnames(counts)
  )
) {

  stop(
    "Metadata and count-matrix cell order do not match."
  )

}


cat(
  "Cell ID validation: PASSED\n"
)


# ------------------------------------------------------------
# 14. Check sample distribution
# ------------------------------------------------------------

sample_ids <- as.character(
  metadata$sample_id
)

sample_counts <- table(
  sample_ids
)

cat("\n")
cat("====================================================\n")
cat("CELLS PER SAMPLE\n")
cat("====================================================\n")

print(
  sample_counts
)

samples <- sort(
  unique(
    sample_ids
  )
)

cat(
  "\nNumber of samples:",
  length(samples),
  "\n"
)


# ------------------------------------------------------------
# 15. Verify all expected samples
# ------------------------------------------------------------

expected_samples <- c(
  "BC03",
  "BC05",
  "BC06",
  "BC08",
  "BC11",
  "BC12",
  "BC14",
  "BC15",
  "BC17"
)

missing_samples <- setdiff(
  expected_samples,
  samples
)

if (
  length(missing_samples) > 0
) {

  stop(
    "Expected samples are missing:\n",
    paste(
      missing_samples,
      collapse = ", "
    )
  )

}


# ------------------------------------------------------------
# 16. Calculate proportional allocation
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("CALCULATING STRATIFIED SAMPLING\n")
cat("====================================================\n")

set.seed(
  RANDOM_SEED
)

available <- as.numeric(
  sample_counts
)

names(available) <- names(
  sample_counts
)


# Proportional allocation

exact_allocation <- (
  available /
    sum(available)
) *
  TARGET_CELLS


# Start with floor values

allocation <- floor(
  exact_allocation
)

names(allocation) <- names(
  available
)


# Guarantee at least one cell per sample

allocation[
  allocation < 1
] <- 1


# ------------------------------------------------------------
# 17. Distribute remaining cells
# ------------------------------------------------------------

remaining <- TARGET_CELLS -
  sum(allocation)

if (
  remaining > 0
) {

  fractional_part <- exact_allocation -
    floor(exact_allocation)

  names(fractional_part) <-
    names(exact_allocation)

  order_fraction <- order(
    fractional_part,
    decreasing = TRUE
  )

  for (
    i in seq_len(remaining)
  ) {

    target_sample <- names(
      fractional_part
    )[
      order_fraction[
        ((i - 1) %%
          length(order_fraction)) + 1
      ]
    ]

    allocation[
      target_sample
    ] <-
      allocation[
        target_sample
      ] + 1

  }

}


# ------------------------------------------------------------
# 18. Check allocation against availability
# ------------------------------------------------------------

if (
  any(
    allocation >
      available
  )
) {

  stop(
    "Sampling allocation exceeds available cells."
  )

}


# ------------------------------------------------------------
# 19. Sampling plan
# ------------------------------------------------------------

sampling_plan <- data.frame(
  sample_id = names(
    available
  ),

  cells_available = as.integer(
    available
  ),

  cells_selected = as.integer(
    allocation[
      names(available)
    ]
  ),

  stringsAsFactors = FALSE
)


sampling_plan$fraction_selected <-
  sampling_plan$cells_selected /
  sampling_plan$cells_available


cat("\n")
cat("Sampling plan:\n")

print(
  sampling_plan
)


cat(
  "\nTotal selected:",
  sum(
    sampling_plan$cells_selected
  ),
  "\n"
)


# ------------------------------------------------------------
# 20. Select cells
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("SELECTING REPRESENTATIVE CELLS\n")
cat("====================================================\n")

selected_cells <- character(
  0
)

for (
  sample_name in samples
) {

  available_cells <- metadata$cell_id[
    metadata$sample_id ==
      sample_name
  ]

  n_select <- allocation[
    sample_name
  ]

  selected <- sample(
    available_cells,
    size = n_select,
    replace = FALSE
  )

  selected_cells <- c(
    selected_cells,
    selected
  )

}


# Preserve count-matrix ordering

selected_cells <- colnames(
  counts
)[
  colnames(counts) %in%
    selected_cells
]


# ------------------------------------------------------------
# 21. Validate selected cells
# ------------------------------------------------------------

if (
  length(selected_cells) !=
    TARGET_CELLS
) {

  stop(
    "Selected cell count is not equal to TARGET_CELLS."
  )

}


cat(
  "Selected cells:",
  length(selected_cells),
  "\n"
)


# ------------------------------------------------------------
# 22. Create subset count matrix
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("CREATING COPYKAT SUBSET\n")
cat("====================================================\n")

subset_counts <- counts[
  ,
  selected_cells,
  drop = FALSE
]

cat(
  "Subset dimensions: ",
  nrow(subset_counts),
  " genes x ",
  ncol(subset_counts),
  " cells\n",
  sep = ""
)

cat(
  "Subset class: ",
  class(subset_counts)[1],
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 23. Create subset metadata
# ------------------------------------------------------------

subset_metadata <- metadata[
  match(
    selected_cells,
    metadata$cell_id
  ),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 24. Validate subset metadata
# ------------------------------------------------------------

if (
  !identical(
    subset_metadata$cell_id,
    colnames(subset_counts)
  )
) {

  stop(
    "Subset metadata and count matrix do not match."
  )

}


cat(
  "Subset metadata validation: PASSED\n"
)


# ------------------------------------------------------------
# 25. Check sample representation
# ------------------------------------------------------------

cat("\n")
cat("Selected cells per sample:\n")

print(
  table(
    subset_metadata$sample_id
  )
)


# ------------------------------------------------------------
# 26. Save selected cell IDs
# ------------------------------------------------------------

selected_cell_file <- file.path(
  subset_dir,
  "GSE228499_CopyKAT_selected_cells.csv"
)

write.csv(
  data.frame(
    cell_id = selected_cells,
    sample_id =
      subset_metadata$sample_id
  ),
  selected_cell_file,
  row.names = FALSE
)


# ------------------------------------------------------------
# 27. Save sampling plan
# ------------------------------------------------------------

sampling_plan_file <- file.path(
  subset_dir,
  "GSE228499_CopyKAT_sampling_plan.csv"
)

write.csv(
  sampling_plan,
  sampling_plan_file,
  row.names = FALSE
)


# ------------------------------------------------------------
# 28. Save subset metadata
# ------------------------------------------------------------

subset_metadata_file <- file.path(
  subset_dir,
  "GSE228499_CopyKAT_subset_metadata.csv"
)

write.csv(
  subset_metadata,
  subset_metadata_file,
  row.names = TRUE
)


# ------------------------------------------------------------
# 29. Save subset count matrix
# ------------------------------------------------------------

subset_count_file <- file.path(
  subset_dir,
  "GSE228499_CopyKAT_subset_counts.rds"
)

saveRDS(
  subset_counts,
  subset_count_file
)


# ------------------------------------------------------------
# 30. Memory information
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("SUBSET MEMORY CHECK\n")
cat("====================================================\n")

cat(
  "Non-zero entries:",
  Matrix::nnzero(
    subset_counts
  ),
  "\n"
)

cat(
  "Object size:",
  format(
    object.size(
      subset_counts
    ),
    units = "MB"
  ),
  "\n"
)


# ------------------------------------------------------------
# 31. Save summary
# ------------------------------------------------------------

summary_file <- file.path(
  subset_dir,
  "GSE228499_CopyKAT_subset_summary.txt"
)

summary_lines <- c(

  "COPYKAT SUBSET SUMMARY",

  "======================",

  paste(
    "Original genes:",
    nrow(counts)
  ),

  paste(
    "Original cells:",
    ncol(counts)
  ),

  paste(
    "Subset genes:",
    nrow(subset_counts)
  ),

  paste(
    "Subset cells:",
    ncol(subset_counts)
  ),

  paste(
    "Samples:",
    length(samples)
  ),

  paste(
    "Random seed:",
    RANDOM_SEED
  ),

  "",

  "Cells per sample:",

  capture.output(
    print(
      table(
        subset_metadata$sample_id
      )
    )
  )
)


writeLines(
  summary_lines,
  summary_file
)


# ------------------------------------------------------------
# 32. Final output
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("SUBSET PREPARATION COMPLETED\n")
cat("====================================================\n")

cat(
  "Original dataset:",
  nrow(counts),
  "genes x",
  ncol(counts),
  "cells\n"
)

cat(
  "CopyKAT subset:",
  nrow(subset_counts),
  "genes x",
  ncol(subset_counts),
  "cells\n"
)

cat(
  "Samples represented:",
  length(
    unique(
      subset_metadata$sample_id
    )
  ),
  "\n"
)

cat("\nSaved files:\n")

cat(
  subset_count_file,
  "\n"
)

cat(
  subset_metadata_file,
  "\n"
)

cat(
  selected_cell_file,
  "\n"
)

cat(
  sampling_plan_file,
  "\n"
)

cat(
  summary_file,
  "\n"
)

cat("\n")
cat("CopyKAT has NOT been run yet.\n")
cat("The subset is ready for the next step.\n")

cat("\n")
cat("====================================================\n")
cat("READY FOR COPYKAT\n")
cat("====================================================\n")


# ------------------------------------------------------------
# 33. Cleanup
# ------------------------------------------------------------

rm(
  counts,
  metadata,
  subset_counts,
  subset_metadata
)

gc()