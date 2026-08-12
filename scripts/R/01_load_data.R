# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 01_load_data.R
# Purpose: Load all 9 scRNA-seq samples into Seurat
# ============================================================


# ------------------------------------------------------------
# 1. Load required package
# ------------------------------------------------------------

library(Seurat)


# ------------------------------------------------------------
# 2. Define project directories
# ------------------------------------------------------------

project_dir <- "~/Projects/SingleCell_BreastCancer_AI"

data_dir <- file.path(
  project_dir,
  "data",
  "raw",
  "GSE228499"
)

metadata_file <- file.path(
  project_dir,
  "data",
  "metadata",
  "sample_metadata.csv"
)


# ------------------------------------------------------------
# 3. Load sample metadata
# ------------------------------------------------------------

sample_metadata <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE
)

print(sample_metadata)


# ------------------------------------------------------------
# 4. Define sample IDs
# ------------------------------------------------------------

sample_ids <- sample_metadata$sample_id

print(sample_ids)


# ------------------------------------------------------------
# 5. Create empty list for Seurat objects
# ------------------------------------------------------------

seurat_list <- list()


# ------------------------------------------------------------
# 6. Load each sample
# ------------------------------------------------------------

for (sample_id in sample_ids) {

  cat("\n========================================\n")
  cat("Loading sample:", sample_id, "\n")
  cat("========================================\n")

  # Path to sample directory
  sample_dir <- file.path(
    data_dir,
    sample_id
  )

  # Read 10x expression matrix
  counts <- Read10X(
    data.dir = sample_dir
  )

  # Display matrix dimensions
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

  # Create Seurat object
  seu <- CreateSeuratObject(
    counts = counts,
    project = paste0("GSE228499_", sample_id)
  )

  # Add sample information
  seu$sample_id <- sample_id

  seu$patient_id <- sample_id

  # Add cancer subtype
  seu$cancer_subtype <-
    sample_metadata$cancer_subtype[
      sample_metadata$sample_id == sample_id
    ]

  # Add treatment information
  seu$treatment <-
    sample_metadata$treatment[
      sample_metadata$sample_id == sample_id
    ]

  # Add GSM accession
  seu$GSM <-
    sample_metadata$GSM[
      sample_metadata$sample_id == sample_id
    ]

  # Store Seurat object in list
  seurat_list[[sample_id]] <- seu

  # Print object summary
  print(seu)
}


# ------------------------------------------------------------
# 7. Check that all samples were loaded
# ------------------------------------------------------------

cat("\n========================================\n")
cat("All samples loaded successfully\n")
cat("========================================\n")

print(names(seurat_list))


# ------------------------------------------------------------
# 8. Display cells per sample
# ------------------------------------------------------------

cell_counts <- sapply(
  seurat_list,
  ncol
)

print(cell_counts)


# ------------------------------------------------------------
# 9. Calculate total number of cells
# ------------------------------------------------------------

total_cells <- sum(cell_counts)

cat(
  "\nTotal cells across all samples:",
  total_cells,
  "\n"
)


# ------------------------------------------------------------
# 10. Merge all samples
# ------------------------------------------------------------

combined <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "GSE228499"
)


# ------------------------------------------------------------
# 11. Display combined object
# ------------------------------------------------------------

cat("\n========================================\n")
cat("Combined Seurat object\n")
cat("========================================\n")

print(combined)


# ------------------------------------------------------------
# 12. Check combined metadata
# ------------------------------------------------------------

head(combined@meta.data)

# ------------------------------------------------------------
# 13. Check number of cells by sample
# ------------------------------------------------------------

table(combined$sample_id)
saveRDS(
  combined,
  file = file.path(
    project_dir,
    "data",
    "processed",
    "GSE228499_merged_raw.rds"
  )
)