# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 11_find_cluster_markers.R
#
# Purpose:
#   Identify marker genes for every scRNA-seq cluster.
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)


# ------------------------------------------------------------
# 1. Project paths
# ------------------------------------------------------------

project_dir <- "~/Projects/SingleCell_BreastCancer_AI"

input_file <- file.path(
  project_dir,
  "data",
  "processed",
  "GSE228499_clustered.rds"
)

result_dir <- file.path(
  project_dir,
  "results",
  "markers"
)

figure_dir <- file.path(
  project_dir,
  "figures",
  "markers"
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
# 2. Load clustered object
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("LOADING CLUSTERED DATASET\n")
cat("====================================================\n")

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


# ------------------------------------------------------------
# 3. Cluster information
# ------------------------------------------------------------

cat("\nCluster sizes:\n")

print(
  table(
    obj$seurat_clusters
  )
)


# ------------------------------------------------------------
# 4. Check RNA layers
# ------------------------------------------------------------

cat("\nRNA layers:\n")

print(
  Layers(
    obj[["RNA"]]
  )
)


# ------------------------------------------------------------
# 5. Join RNA layers for differential expression
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("JOINING RNA LAYERS FOR MARKER ANALYSIS\n")
cat("====================================================\n")

obj <- JoinLayers(
  obj,
  assay = "RNA"
)

cat("\nRNA layers after JoinLayers:\n")

print(
  Layers(
    obj[["RNA"]]
  )
)


# ------------------------------------------------------------
# 6. Find markers for every cluster
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("FINDING CLUSTER MARKERS\n")
cat("====================================================\n")

markers <- FindAllMarkers(
  object = obj,
  assay = "RNA",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = TRUE
)


# ------------------------------------------------------------
# 7. Save complete marker table
# ------------------------------------------------------------

write.csv(
  markers,
  file.path(
    result_dir,
    "all_cluster_markers.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 8. Top 10 markers per cluster
# ------------------------------------------------------------

top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10
  ) %>%
  ungroup()


write.csv(
  top10,
  file.path(
    result_dir,
    "top10_markers_per_cluster.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 9. Print top markers
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("TOP MARKERS PER CLUSTER\n")
cat("====================================================\n")

print(
  top10 %>%
    select(
      cluster,
      gene,
      avg_log2FC,
      pct.1,
      pct.2,
      p_val_adj
    )
)


# ------------------------------------------------------------
# 10. Heatmap of top markers
# ------------------------------------------------------------

top_heatmap_genes <- markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 5
  ) %>%
  pull(gene) %>%
  unique()


pdf(
  file.path(
    figure_dir,
    "cluster_marker_heatmap_top5.pdf"
  ),
  width = 12,
  height = 10
)

print(
  DoHeatmap(
    obj,
    features = top_heatmap_genes,
    group.by = "seurat_clusters"
  )
)

dev.off()


# ------------------------------------------------------------
# 11. Dot plot of top markers
# ------------------------------------------------------------

pdf(
  file.path(
    figure_dir,
    "cluster_marker_dotplot.pdf"
  ),
  width = 14,
  height = 10
)

print(
  DotPlot(
    obj,
    features = top_heatmap_genes,
    group.by = "seurat_clusters"
  ) +
    RotatedAxis()
)

dev.off()


# ------------------------------------------------------------
# 12. Final summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("MARKER ANALYSIS COMPLETED\n")
cat("====================================================\n")

cat(
  "Total marker genes:",
  nrow(markers),
  "\n"
)

cat(
  "Clusters analyzed:",
  length(
    unique(
      markers$cluster
    )
  ),
  "\n"
)

cat("\nResults:\n")

cat(
  file.path(
    result_dir,
    "all_cluster_markers.csv"
  ),
  "\n"
)

cat(
  file.path(
    result_dir,
    "top10_markers_per_cluster.csv"
  ),
  "\n"
)

cat("\nFigures:\n")

cat(
  file.path(
    figure_dir,
    "cluster_marker_heatmap_top5.pdf"
  ),
  "\n"
)

cat(
  file.path(
    figure_dir,
    "cluster_marker_dotplot.pdf"
  ),
  "\n"
)

cat("\n====================================================\n")
cat("READY FOR CELL-TYPE ANNOTATION\n")
cat("====================================================\n")