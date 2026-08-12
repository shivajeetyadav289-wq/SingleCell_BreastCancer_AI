# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 02_quality_control.R
# Purpose: Calculate and visualize scRNA-seq QC metrics
# ============================================================


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

library(Seurat)
library(ggplot2)


# ------------------------------------------------------------
# 2. Define project directories
# ------------------------------------------------------------

project_dir <- "~/Projects/SingleCell_BreastCancer_AI"

processed_dir <- file.path(
  project_dir,
  "data",
  "processed"
)

qc_dir <- file.path(
  project_dir,
  "results",
  "qc"
)

figure_dir <- file.path(
  project_dir,
  "figures",
  "qc"
)


# ------------------------------------------------------------
# 3. Create output directories if they do not exist
# ------------------------------------------------------------

dir.create(
  qc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 4. Load merged raw Seurat object
# ------------------------------------------------------------

combined <- readRDS(
  file.path(
    processed_dir,
    "GSE228499_merged_raw.rds"
  )
)


# ------------------------------------------------------------
# 5. Confirm the loaded object
# ------------------------------------------------------------

cat("\n========================================\n")
cat("Loaded Seurat object\n")
cat("========================================\n")

print(combined)

cat("\nNumber of genes:",
    nrow(combined),
    "\n")

cat("Number of cells:",
    ncol(combined),
    "\n")


# ------------------------------------------------------------
# 6. Calculate mitochondrial percentage
# ------------------------------------------------------------

combined[["percent.mt"]] <- PercentageFeatureSet(
  combined,
  pattern = "^MT-"
)


# ------------------------------------------------------------
# 7. Display overall QC summaries
# ------------------------------------------------------------

cat("\n========================================\n")
cat("Overall QC Summary\n")
cat("========================================\n")

cat("\nnFeature_RNA:\n")
print(summary(combined$nFeature_RNA))

cat("\nnCount_RNA:\n")
print(summary(combined$nCount_RNA))

cat("\npercent.mt:\n")
print(summary(combined$percent.mt))


# ------------------------------------------------------------
# 8. Save overall QC summary
# ------------------------------------------------------------

qc_summary <- data.frame(
  metric = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  minimum = c(
    min(combined$nFeature_RNA),
    min(combined$nCount_RNA),
    min(combined$percent.mt)
  ),
  Q1 = c(
    quantile(combined$nFeature_RNA, 0.25),
    quantile(combined$nCount_RNA, 0.25),
    quantile(combined$percent.mt, 0.25)
  ),
  median = c(
    median(combined$nFeature_RNA),
    median(combined$nCount_RNA),
    median(combined$percent.mt)
  ),
  mean = c(
    mean(combined$nFeature_RNA),
    mean(combined$nCount_RNA),
    mean(combined$percent.mt)
  ),
  Q3 = c(
    quantile(combined$nFeature_RNA, 0.75),
    quantile(combined$nCount_RNA, 0.75),
    quantile(combined$percent.mt, 0.75)
  ),
  maximum = c(
    max(combined$nFeature_RNA),
    max(combined$nCount_RNA),
    max(combined$percent.mt)
  )
)

write.csv(
  qc_summary,
  file.path(
    qc_dir,
    "overall_qc_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. QC summary by sample
# ------------------------------------------------------------

qc_by_sample <- aggregate(
  cbind(
    nFeature_RNA,
    nCount_RNA,
    percent.mt
  ) ~ sample_id,
  data = combined@meta.data,
  FUN = median
)

cat("\n========================================\n")
cat("Median QC Metrics by Sample\n")
cat("========================================\n")

print(qc_by_sample)


# ------------------------------------------------------------
# 10. Save QC summary by sample
# ------------------------------------------------------------

write.csv(
  qc_by_sample,
  file.path(
    qc_dir,
    "qc_summary_by_sample.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 11. Number of cells per sample before QC
# ------------------------------------------------------------

cells_per_sample <- as.data.frame(
  table(combined$sample_id)
)

colnames(cells_per_sample) <- c(
  "sample_id",
  "cell_count"
)

cat("\n========================================\n")
cat("Cells per Sample Before QC\n")
cat("========================================\n")

print(cells_per_sample)


# ------------------------------------------------------------
# 12. Save cell counts
# ------------------------------------------------------------

write.csv(
  cells_per_sample,
  file.path(
    qc_dir,
    "cells_per_sample_before_qc.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. Violin plots — all cells
# ------------------------------------------------------------

pdf(
  file.path(
    figure_dir,
    "QC_violin_before_filtering.pdf"
  ),
  width = 12,
  height = 5
)

print(
  VlnPlot(
    combined,
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt"
    ),
    ncol = 3,
    pt.size = 0
  )
)

dev.off()


# ------------------------------------------------------------
# 14. QC by sample
# ------------------------------------------------------------

pdf(
  file.path(
    figure_dir,
    "QC_by_sample.pdf"
  ),
  width = 12,
  height = 10
)

print(
  VlnPlot(
    combined,
    features = c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt"
    ),
    group.by = "sample_id",
    ncol = 1,
    pt.size = 0
  )
)

dev.off()


# ------------------------------------------------------------
# 15. Scatter plot:
#     nCount_RNA vs nFeature_RNA
# ------------------------------------------------------------

pdf(
  file.path(
    figure_dir,
    "QC_scatter_counts_vs_features.pdf"
  ),
  width = 7,
  height = 6
)

print(
  FeatureScatter(
    combined,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
  )
)

dev.off()


# ------------------------------------------------------------
# 16. Scatter plot:
#     percent.mt vs nFeature_RNA
# ------------------------------------------------------------

pdf(
  file.path(
    figure_dir,
    "QC_scatter_mt_vs_features.pdf"
  ),
  width = 7,
  height = 6
)

print(
  FeatureScatter(
    combined,
    feature1 = "percent.mt",
    feature2 = "nFeature_RNA"
  )
)

dev.off()


# ------------------------------------------------------------
# 17. Final QC message
# ------------------------------------------------------------

cat("\n========================================\n")
cat("QC assessment completed successfully\n")
cat("========================================\n")

cat("\nNo cells were filtered in this script.\n")

cat("\nQC results saved to:\n")
cat(qc_dir, "\n")

cat("\nQC figures saved to:\n")
cat(figure_dir, "\n")

# ------------------------------------------------------------
# 18. Evaluate proposed QC thresholds
# ------------------------------------------------------------

min_features <- 200
max_mito <- 20

# Determine which cells pass QC
pass_qc <- (
  combined$nFeature_RNA >= min_features &
  combined$percent.mt < max_mito
)


# ------------------------------------------------------------
# 19. Overall QC retention
# ------------------------------------------------------------

cat("\n========================================\n")
cat("Proposed QC Thresholds\n")
cat("========================================\n")

cat(
  "\nMinimum genes per cell:",
  min_features,
  "\n"
)

cat(
  "Maximum mitochondrial percentage:",
  max_mito,
  "%\n"
)

cat(
  "\nCells before QC:",
  ncol(combined),
  "\n"
)

cat(
  "Cells passing QC:",
  sum(pass_qc),
  "\n"
)

cat(
  "Cells removed:",
  sum(!pass_qc),
  "\n"
)

cat(
  "Percentage retained:",
  round(mean(pass_qc) * 100, 2),
  "%\n"
)


# ------------------------------------------------------------
# 20. QC retention by sample
# ------------------------------------------------------------

retention_data <- data.frame(
  sample_id = combined$sample_id,
  pass_qc = pass_qc
)

qc_retention_summary <- aggregate(
  pass_qc ~ sample_id,
  data = retention_data,
  FUN = function(x) {
    c(
      total = length(x),
      retained = sum(x),
      removed = sum(!x),
      retention_percent = mean(x) * 100
    )
  }
)

# Convert aggregate matrix into readable columns
qc_retention_summary <- data.frame(
  sample_id = qc_retention_summary$sample_id,
  total_cells = qc_retention_summary$pass_qc[, "total"],
  retained_cells = qc_retention_summary$pass_qc[, "retained"],
  removed_cells = qc_retention_summary$pass_qc[, "removed"],
  retention_percent = qc_retention_summary$pass_qc[, "retention_percent"]
)

print(qc_retention_summary)


# ------------------------------------------------------------
# 21. Save QC retention summary
# ------------------------------------------------------------

write.csv(
  qc_retention_summary,
  file.path(
    qc_dir,
    "qc_retention_by_sample.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 22. Determine reason for QC failure
# ------------------------------------------------------------

low_features <- combined$nFeature_RNA < min_features

high_mito <- combined$percent.mt >= max_mito

qc_failure <- data.frame(
  sample_id = combined$sample_id,
  low_features = low_features,
  high_mito = high_mito,
  fail_qc = !pass_qc
)


# ------------------------------------------------------------
# 23. Summarize QC failure reasons by sample
# ------------------------------------------------------------

qc_failure_summary <- aggregate(
  cbind(
    low_features,
    high_mito,
    fail_qc
  ) ~ sample_id,
  data = qc_failure,
  FUN = sum
)

print(qc_failure_summary)


# ------------------------------------------------------------
# 24. Save QC failure summary
# ------------------------------------------------------------

write.csv(
  qc_failure_summary,
  file.path(
    qc_dir,
    "qc_failure_reasons_by_sample.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 25. QC metric quantiles by sample
# ------------------------------------------------------------

qc_quantiles <- do.call(
  rbind,
  lapply(
    split(
      combined@meta.data,
      combined$sample_id
    ),
    function(df) {
      data.frame(
        nFeature_Q01 = quantile(df$nFeature_RNA, 0.01),
        nFeature_Q05 = quantile(df$nFeature_RNA, 0.05),
        nFeature_Q10 = quantile(df$nFeature_RNA, 0.10),
        nFeature_median = median(df$nFeature_RNA),
        nFeature_Q90 = quantile(df$nFeature_RNA, 0.90),
        mt_Q90 = quantile(df$percent.mt, 0.90),
        mt_Q95 = quantile(df$percent.mt, 0.95),
        mt_Q99 = quantile(df$percent.mt, 0.99)
      )
    }
  )
)

qc_quantiles$sample_id <- rownames(qc_quantiles)

rownames(qc_quantiles) <- NULL

print(qc_quantiles)

write.csv(
  qc_quantiles,
  file.path(
    qc_dir,
    "qc_quantiles_by_sample.csv"
  ),
  row.names = FALSE
)