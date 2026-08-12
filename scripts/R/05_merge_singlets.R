# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 05_merge_singlets.R
#
# Purpose:
#   Merge the 9 sample-level singlet Seurat objects
# ============================================================


# ------------------------------------------------------------
# 1. Load Seurat
# ------------------------------------------------------------

library(Seurat)


# ------------------------------------------------------------
# 2. Project directories
# ------------------------------------------------------------

project_dir <- "~/Projects/SingleCell_BreastCancer_AI"

doublet_dir <- file.path(
  project_dir,
  "results",
  "doublet_detection"
)

processed_dir <- file.path(
  project_dir,
  "data",
  "processed"
)

results_dir <- file.path(
  project_dir,
  "results"
)

dir.create(
  processed_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. Define samples
# ------------------------------------------------------------

samples <- c(
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


# ------------------------------------------------------------
# 4. Load singlet objects
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("Loading singlet objects\n")
cat("====================================================\n")


singlet_objects <- lapply(
  samples,
  function(sample_name) {

    file_path <- file.path(
      doublet_dir,
      paste0(
        sample_name,
        "_singlets.rds"
      )
    )

    cat(
      "Loading:",
      sample_name,
      "\n"
    )

    if (!file.exists(file_path)) {

      stop(
        paste(
          "Missing file:",
          file_path
        )
      )
    }

    readRDS(file_path)
  }
)


names(singlet_objects) <- samples


# ------------------------------------------------------------
# 5. Check cells in each sample
# ------------------------------------------------------------

cat("\n")
cat("Cells per sample:\n")

cells_per_sample <- sapply(
  singlet_objects,
  ncol
)

print(
  cells_per_sample
)


# ------------------------------------------------------------
# 6. Check total cells
# ------------------------------------------------------------

total_singlets <- sum(
  cells_per_sample
)

cat(
  "\nTotal singlet cells:",
  total_singlets,
  "\n"
)


# ------------------------------------------------------------
# 7. Merge all samples
# ------------------------------------------------------------

cat("\n")
cat("Merging singlet objects...\n")


merged_singlets <- merge(
  x = singlet_objects[[1]],
  y = singlet_objects[2:length(singlet_objects)],
  add.cell.ids = samples,
  project = "GSE228499"
)


# ------------------------------------------------------------
# 8. Check dimensions
# ------------------------------------------------------------

cat("\n")
cat("Merged object dimensions:\n")

print(
  dim(merged_singlets)
)


# ------------------------------------------------------------
# 9. Check sample metadata
# ------------------------------------------------------------

cat("\n")
cat("Cells by sample after merging:\n")

merged_sample_counts <- table(
  merged_singlets$sample_id
)

print(
  merged_sample_counts
)


# ------------------------------------------------------------
# 10. Check patient metadata
# ------------------------------------------------------------

cat("\n")
cat("Patients:\n")

print(
  unique(
    merged_singlets$patient_id
  )
)


# ------------------------------------------------------------
# 11. Check RNA layers
# ------------------------------------------------------------

cat("\n")
cat("RNA layers:\n")

print(
  Layers(
    merged_singlets[["RNA"]]
  )
)


# ------------------------------------------------------------
# 12. Save merged singlet object
# ------------------------------------------------------------

output_file <- file.path(
  processed_dir,
  "GSE228499_singlets.rds"
)

saveRDS(
  merged_singlets,
  output_file
)


cat("\n")
cat("Saved:\n")
cat(
  output_file,
  "\n"
)


# ------------------------------------------------------------
# 13. Save sample cell counts
# ------------------------------------------------------------

sample_count_df <- data.frame(
  sample_id = names(
    merged_sample_counts
  ),
  cells = as.numeric(
    merged_sample_counts
  )
)

write.csv(
  sample_count_df,
  file.path(
    results_dir,
    "singlet_cells_by_sample.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. Final summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("MERGING COMPLETED\n")
cat("====================================================\n")

cat(
  "Genes:",
  nrow(merged_singlets),
  "\n"
)

cat(
  "Cells:",
  ncol(merged_singlets),
  "\n"
)

cat(
  "Samples:",
  length(
    unique(
      merged_singlets$sample_id
    )
  ),
  "\n"
)

cat(
  "Output:",
  output_file,
  "\n"
)