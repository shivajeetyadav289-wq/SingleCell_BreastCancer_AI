# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 04d_doublet_detection_all_samples.R
#
# Purpose:
#   Automated sample-wise DoubletFinder analysis
#
# Strategy:
#   1. Process each patient independently
#   2. Join Seurat v5 RNA layers
#   3. Normalize and identify variable genes
#   4. Scale and run PCA
#   5. Determine optimal pK automatically
#   6. Run DoubletFinder
#   7. Save doublet results
#   8. Save singlet cells
#   9. Skip samples already completed
#  10. Save summary after each sample
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(Seurat)
library(DoubletFinder)


# ------------------------------------------------------------
# 2. Project directories
# ------------------------------------------------------------

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
# 3. Load QC-filtered dataset
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("Loading QC-filtered dataset\n")
cat("====================================================\n")

combined <- readRDS(
  file.path(
    processed_dir,
    "GSE228499_filtered.rds"
  )
)

cat(
  "Total cells:",
  ncol(combined),
  "\n"
)

cat(
  "Total genes:",
  nrow(combined),
  "\n"
)


# ------------------------------------------------------------
# 4. Split by patient/sample
# ------------------------------------------------------------

sample_objects <- SplitObject(
  combined,
  split.by = "sample_id"
)

samples <- names(sample_objects)

cat("\nSamples to process:\n")
print(samples)


# ------------------------------------------------------------
# 5. Process each sample
# ------------------------------------------------------------

for (sample_name in samples) {

  cat("\n\n")
  cat("====================================================\n")
  cat("Processing sample:", sample_name, "\n")
  cat("====================================================\n")


  # ----------------------------------------------------------
  # Define output files
  # ----------------------------------------------------------

  singlet_file <- file.path(
    doublet_dir,
    paste0(
      sample_name,
      "_singlets.rds"
    )
  )

  result_file <- file.path(
    doublet_dir,
    paste0(
      sample_name,
      "_doublet_results.rds"
    )
  )

  pK_file <- file.path(
    doublet_dir,
    paste0(
      sample_name,
      "_pK_results.csv"
    )
  )

  summary_file <- file.path(
    doublet_dir,
    paste0(
      sample_name,
      "_summary.csv"
    )
  )


  # ----------------------------------------------------------
  # Skip completed samples
  # ----------------------------------------------------------

  if (
    file.exists(singlet_file) &&
    file.exists(result_file) &&
    file.exists(summary_file)
  ) {

    cat(
      "Sample already completed:",
      sample_name,
      "\n"
    )

    cat("Skipping...\n")

    next
  }


  # ----------------------------------------------------------
  # Extract sample
  # ----------------------------------------------------------

  obj <- sample_objects[[sample_name]]

  cat(
    "Cells before DoubletFinder:",
    ncol(obj),
    "\n"
  )


  # ----------------------------------------------------------
  # Join Seurat v5 RNA layers
  # ----------------------------------------------------------

  rna_layers <- Layers(
    obj[["RNA"]]
  )

  cat(
    "RNA layers before JoinLayers:",
    paste(
      rna_layers,
      collapse = ", "
    ),
    "\n"
  )


  if (
    length(rna_layers) > 1
  ) {

    cat(
      "Joining RNA layers...\n"
    )

    obj <- JoinLayers(
      object = obj,
      assay = "RNA"
    )
  }


  cat(
    "RNA layers after JoinLayers:",
    paste(
      Layers(obj[["RNA"]]),
      collapse = ", "
    ),
    "\n"
  )


  # ----------------------------------------------------------
  # Normalize
  # ----------------------------------------------------------

  cat(
    "Running NormalizeData...\n"
  )

  obj <- NormalizeData(
    object = obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )


  # ----------------------------------------------------------
  # Find variable genes
  # ----------------------------------------------------------

  cat(
    "Finding variable genes...\n"
  )

  obj <- FindVariableFeatures(
    object = obj,
    selection.method = "vst",
    nfeatures = 2000
  )


  # ----------------------------------------------------------
  # Scale
  # ----------------------------------------------------------

  cat(
    "Running ScaleData...\n"
  )

  obj <- ScaleData(
    object = obj,
    features = VariableFeatures(obj)
  )


  # ----------------------------------------------------------
  # PCA
  # ----------------------------------------------------------

  cat(
    "Running PCA...\n"
  )

  obj <- RunPCA(
    object = obj,
    features = VariableFeatures(obj),
    npcs = 20,
    verbose = FALSE
  )


  # ----------------------------------------------------------
  # Verify Seurat command history
  # ----------------------------------------------------------

  required_commands <- c(
    "NormalizeData.RNA",
    "FindVariableFeatures.RNA",
    "ScaleData.RNA",
    "RunPCA.RNA"
  )

  missing_commands <- setdiff(
    required_commands,
    names(obj@commands)
  )

  if (
    length(missing_commands) > 0
  ) {

    stop(
      paste(
        "Required Seurat commands missing for",
        sample_name,
        ":",
        paste(
          missing_commands,
          collapse = ", "
        )
      )
    )
  }


  cat(
    "Seurat preprocessing commands verified.\n"
  )


  # ----------------------------------------------------------
  # DoubletFinder parameter sweep
  # ----------------------------------------------------------

  cat("\n")
  cat(
    "Running DoubletFinder parameter sweep...\n"
  )

  sweep_res <- paramSweep(
    obj,
    PCs = 1:20,
    sct = FALSE
  )


  sweep_stats <- summarizeSweep(
    sweep_res,
    GT = FALSE
  )


  bcmvn <- find.pK(
    sweep_stats
  )


  # ----------------------------------------------------------
  # Select optimal pK
  # ----------------------------------------------------------

  best_index <- which.max(
    bcmvn$BCmetric
  )

  best_row <- bcmvn[
    best_index,
    ,
    drop = FALSE
  ]

  best_pK <- as.numeric(
    as.character(
      best_row$pK
    )
  )

  best_BCmetric <- as.numeric(
    best_row$BCmetric
  )


  cat(
    "\nOptimal pK:",
    best_pK,
    "\n"
  )

  cat(
    "Best BCmetric:",
    best_BCmetric,
    "\n"
  )


  # ----------------------------------------------------------
  # Save complete pK table
  # ----------------------------------------------------------

  write.csv(
    bcmvn,
    pK_file,
    row.names = FALSE
  )


  # ----------------------------------------------------------
  # Expected doublet rate
  #
  # Initial project estimate = 7.5%
  # ----------------------------------------------------------

  expected_rate <- 0.075

  nExp <- round(
    expected_rate * ncol(obj)
  )


  cat(
    "Expected doublets:",
    nExp,
    "\n"
  )


  # ----------------------------------------------------------
  # Run DoubletFinder
  # ----------------------------------------------------------

  cat("\n")
  cat(
    "Running DoubletFinder...\n"
  )

  obj <- doubletFinder(
    obj,
    PCs = 1:20,
    pN = 0.25,
    pK = best_pK,
    nExp = nExp,
    reuse.pANN = NULL,
    sct = FALSE
  )


  # ----------------------------------------------------------
  # Identify DoubletFinder classification column
  # ----------------------------------------------------------

  df_columns <- grep(
    "^DF.classifications",
    colnames(obj@meta.data),
    value = TRUE
  )

  if (
    length(df_columns) == 0
  ) {

    stop(
      paste(
        "DoubletFinder classification column not found for",
        sample_name
      )
    )
  }


  df_column <- tail(
    df_columns,
    1
  )


  cat(
    "DoubletFinder classification column:",
    df_column,
    "\n"
  )


  # ----------------------------------------------------------
  # Store simplified classification
  # ----------------------------------------------------------

  obj$doublet_status <- obj@meta.data[[df_column]]


  # ----------------------------------------------------------
  # Classification summary
  # ----------------------------------------------------------

  classification_table <- table(
    obj$doublet_status
  )

  cat("\n")
  cat("Classification:\n")
  print(
    classification_table
  )


  classification_percent <- prop.table(
    classification_table
  ) * 100

  cat("\n")
  cat("Classification percentage:\n")
  print(
    round(
      classification_percent,
      2
    )
  )


  # ----------------------------------------------------------
  # Save complete DoubletFinder object
  # ----------------------------------------------------------

  saveRDS(
    obj,
    result_file
  )


  # ----------------------------------------------------------
  # Create singlet object
  # ----------------------------------------------------------

  singlets <- subset(
    x = obj,
    subset = doublet_status == "Singlet"
  )


  # ----------------------------------------------------------
  # Save singlets
  # ----------------------------------------------------------

  saveRDS(
    singlets,
    singlet_file
  )


  # ----------------------------------------------------------
  # Create sample summary
  # ----------------------------------------------------------

  summary_row <- data.frame(
    sample_id = sample_name,

    cells_before = ncol(obj),

    doublets = sum(
      obj$doublet_status == "Doublet"
    ),

    singlets = sum(
      obj$doublet_status == "Singlet"
    ),

    doublet_percent =
      mean(
        obj$doublet_status == "Doublet"
      ) * 100,

    pK = best_pK,

    BCmetric = best_BCmetric,

    expected_doublet_rate =
      expected_rate,

    stringsAsFactors = FALSE
  )


  # ----------------------------------------------------------
  # Save sample summary
  # ----------------------------------------------------------

  write.csv(
    summary_row,
    summary_file,
    row.names = FALSE
  )


  # ----------------------------------------------------------
  # Confirm checkpoint files
  # ----------------------------------------------------------

  cat("\n")
  cat(
    "Checkpoint files saved for:",
    sample_name,
    "\n"
  )

  cat(
    "  Doublet result:",
    result_file,
    "\n"
  )

  cat(
    "  Singlets:",
    singlet_file,
    "\n"
  )

  cat(
    "  pK results:",
    pK_file,
    "\n"
  )

  cat(
    "  Summary:",
    summary_file,
    "\n"
  )


  # ----------------------------------------------------------
  # Free memory
  # ----------------------------------------------------------

  rm(
    obj,
    singlets,
    sweep_res,
    sweep_stats,
    bcmvn
  )

  gc()


  cat("\n")
  cat(
    "Completed:",
    sample_name,
    "\n"
  )

  cat(
    "Moving to next sample...\n"
  )
}


# ------------------------------------------------------------
# 6. Combine all available sample summaries
# ------------------------------------------------------------

summary_files <- list.files(
  doublet_dir,
  pattern = "_summary\\.csv$",
  full.names = TRUE
)

if (
  length(summary_files) > 0
) {

  summary_list <- lapply(
    summary_files,
    read.csv
  )

  summary_table <- do.call(
    rbind,
    summary_list
  )

  write.csv(
    summary_table,
    file.path(
      doublet_dir,
      "DoubletFinder_all_samples_summary.csv"
    ),
    row.names = FALSE
  )


  cat("\n\n")
  cat("====================================================\n")
  cat("DOUBLETFINDER SUMMARY\n")
  cat("====================================================\n")

  print(
    summary_table
  )
}


# ------------------------------------------------------------
# 7. Final message
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("DOUBLETFINDER PIPELINE COMPLETED\n")
cat("====================================================\n")