# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 04_doublet_detection.R
# Purpose: Detect potential doublets using DoubletFinder
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(Seurat)
library(DoubletFinder)


# ------------------------------------------------------------
# 2. Define project directories
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

figure_dir <- file.path(
  project_dir,
  "figures",
  "doublet_detection"
)


# ------------------------------------------------------------
# 3. Create output directories
# ------------------------------------------------------------

dir.create(
  doublet_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 4. Load QC-filtered Seurat object
# ------------------------------------------------------------

combined <- readRDS(
  file.path(
    processed_dir,
    "GSE228499_filtered.rds"
  )
)


# ------------------------------------------------------------
# 5. Confirm input object
# ------------------------------------------------------------

cat("\n========================================\n")
cat("Input Seurat object\n")
cat("========================================\n")

print(combined)

cat(
  "\nCells:",
  ncol(combined),
  "\n"
)

cat(
  "Genes:",
  nrow(combined),
  "\n"
)


# ------------------------------------------------------------
# 6. Normalize RNA data
# ------------------------------------------------------------

cat("\n========================================\n")
cat("Normalization\n")
cat("========================================\n")

combined <- NormalizeData(
  combined,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


# ------------------------------------------------------------
# 7. Identify highly variable genes
# ------------------------------------------------------------

cat("\n========================================\n")
cat("Finding variable features\n")
cat("========================================\n")

combined <- FindVariableFeatures(
  combined,
  selection.method = "vst",
  nfeatures = 2000
)


cat(
  "\nNumber of variable features:",
  length(VariableFeatures(combined)),
  "\n"
)


# ------------------------------------------------------------
# 8. Scale the data
# ------------------------------------------------------------

cat("\n========================================\n")
cat("Scaling data\n")
cat("========================================\n")

combined <- ScaleData(
  combined,
  features = VariableFeatures(combined)
)


# ------------------------------------------------------------
# 9. Principal Component Analysis
# ------------------------------------------------------------

cat("\n========================================\n")
cat("Running PCA\n")
cat("========================================\n")

combined <- RunPCA(
  combined,
  features = VariableFeatures(combined),
  npcs = 30,
  verbose = FALSE
)


# ------------------------------------------------------------
# 10. Save PCA elbow plot
# ------------------------------------------------------------

pdf(
  file.path(
    figure_dir,
    "PCA_elbow_plot.pdf"
  ),
  width = 7,
  height = 5
)

print(
  ElbowPlot(
    combined,
    ndims = 30
  )
)

dev.off()


# ------------------------------------------------------------
# 11. Parameter sweep for pK
# ------------------------------------------------------------

cat("\n========================================\n")
cat("Running DoubletFinder parameter sweep\n")
cat("========================================\n")

sweep_res <- paramSweep(
  combined,
  PCs = 1:20,
  sct = FALSE
)


# ------------------------------------------------------------
# 12. Summarize parameter sweep
# ------------------------------------------------------------

sweep_stats <- summarizeSweep(
  sweep_res,
  GT = FALSE
)


# ------------------------------------------------------------
# 13. Identify optimal pK
# ------------------------------------------------------------

bcmvn <- find.pK(
  sweep_stats
)


# ------------------------------------------------------------
# 14. Save pK results
# ------------------------------------------------------------

write.csv(
  bcmvn,
  file.path(
    doublet_dir,
    "DoubletFinder_pK_results.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 15. Save pK diagnostic plot
# ------------------------------------------------------------

pdf(
  file.path(
    figure_dir,
    "DoubletFinder_pK_BCmvn.pdf"
  ),
  width = 8,
  height = 6
)

plot(
  bcmvn$BCmetric,
  type = "b",
  xlab = "pK index",
  ylab = "BCmetric",
  main = "DoubletFinder pK Optimization"
)

dev.off()


# ------------------------------------------------------------
# 16. Print pK results
# ------------------------------------------------------------

cat("\n========================================\n")
cat("pK optimization results\n")
cat("========================================\n")

print(bcmvn)


# ------------------------------------------------------------
# 17. Select pK with maximum BCmetric
# ------------------------------------------------------------

best_pK_index <- which.max(
  bcmvn$BCmetric
)

best_pK <- as.numeric(
  as.character(
    bcmvn$pK[best_pK_index]
  )
)

cat(
  "\nSelected pK:",
  best_pK,
  "\n"
)


# ------------------------------------------------------------
# 18. Estimate expected doublets
# ------------------------------------------------------------

# Initial expected doublet rate.
# This will be refined using the experimental loading information.

expected_doublet_rate <- 0.075

nExp_poi <- round(
  expected_doublet_rate * ncol(combined)
)

cat(
  "\nInitial expected doublets:",
  nExp_poi,
  "\n"
)


# ------------------------------------------------------------
# 19. Run DoubletFinder
# ------------------------------------------------------------

cat("\n========================================\n")
cat("Running DoubletFinder\n")
cat("========================================\n")

combined <- doubletFinder(
  combined,
  PCs = 1:20,
  pN = 0.25,
  pK = best_pK,
  nExp = nExp_poi,
  reuse.pANN = NULL,
  sct = FALSE
)


# ------------------------------------------------------------
# 20. Identify DoubletFinder classification column
# ------------------------------------------------------------

df_columns <- grep(
  "^DF.classifications",
  colnames(combined@meta.data),
  value = TRUE
)

df_column <- tail(
  df_columns,
  1
)

cat(
  "\nDoubletFinder classification column:",
  df_column,
  "\n"
)


# ------------------------------------------------------------
# 21. Display classification results
# ------------------------------------------------------------

cat("\n========================================\n")
cat("DoubletFinder classification\n")
cat("========================================\n")

print(
  table(
    combined@meta.data[[df_column]]
  )
)


# ------------------------------------------------------------
# 22. Save classification summary
# ------------------------------------------------------------

classification_table <- as.data.frame(
  table(
    combined@meta.data[[df_column]]
  )
)

colnames(classification_table) <- c(
  "classification",
  "cell_count"
)

write.csv(
  classification_table,
  file.path(
    doublet_dir,
    "doublet_classification_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 23. Add a simple standardized doublet column
# ------------------------------------------------------------

combined$doublet_status <- combined@meta.data[[df_column]]


# ------------------------------------------------------------
# 24. Save object containing DoubletFinder results
# ------------------------------------------------------------

saveRDS(
  combined,
  file.path(
    processed_dir,
    "GSE228499_doublet_results.rds"
  )
)


# ------------------------------------------------------------
# 25. Final message
# ------------------------------------------------------------

cat("\n========================================\n")
cat("DoubletFinder analysis completed\n")
cat("========================================\n")

cat(
  "\nResult saved to:\n",
  file.path(
    processed_dir,
    "GSE228499_doublet_results.rds"
  ),
  "\n"
)