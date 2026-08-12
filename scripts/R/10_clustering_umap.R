# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 10_clustering_umap.R
#
# Purpose:
#   Perform graph-based clustering and UMAP visualization
#   using the integrated RPCA representation.
# ============================================================

library(Seurat)
library(ggplot2)


# ------------------------------------------------------------
# 1. Project paths
# ------------------------------------------------------------

project_dir <- "~/Projects/SingleCell_BreastCancer_AI"

input_file <- file.path(
  project_dir,
  "data",
  "processed",
  "GSE228499_integrated_rpca.rds"
)

output_file <- file.path(
  project_dir,
  "data",
  "processed",
  "GSE228499_clustered.rds"
)

figure_dir <- file.path(
  project_dir,
  "figures",
  "clustering"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 2. Load integrated object
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("LOADING INTEGRATED DATASET\n")
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
# 3. Check integrated reduction
# ------------------------------------------------------------

cat("\nAvailable reductions:\n")

print(
  Reductions(obj)
)

if (
  !"integrated.rpca" %in% Reductions(obj)
) {
  stop(
    "integrated.rpca reduction not found."
  )
}


# ------------------------------------------------------------
# 4. Define dimensions
# ------------------------------------------------------------

dims_use <- 1:20

cat(
  "\nUsing",
  length(dims_use),
  "integrated PCs for clustering.\n"
)


# ------------------------------------------------------------
# 5. Find nearest neighbors
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("FINDING NEIGHBORS\n")
cat("====================================================\n")

obj <- FindNeighbors(
  object = obj,
  reduction = "integrated.rpca",
  dims = dims_use,
  verbose = TRUE
)


# ------------------------------------------------------------
# 6. Cluster cells
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("CLUSTERING CELLS\n")
cat("====================================================\n")

obj <- FindClusters(
  object = obj,
  resolution = 0.5,
  verbose = TRUE
)


# ------------------------------------------------------------
# 7. Display cluster sizes
# ------------------------------------------------------------

cat("\nCluster sizes:\n")

print(
  table(
    obj$seurat_clusters
  )
)


# ------------------------------------------------------------
# 8. Run UMAP
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("RUNNING UMAP\n")
cat("====================================================\n")

obj <- RunUMAP(
  object = obj,
  reduction = "integrated.rpca",
  dims = dims_use,
  reduction.name = "umap.rpca",
  reduction.key = "UMAPR",
  verbose = TRUE
)


# ------------------------------------------------------------
# 9. UMAP by cluster
# ------------------------------------------------------------

pdf(
  file.path(
    figure_dir,
    "UMAP_by_cluster.pdf"
  ),
  width = 9,
  height = 7
)

print(
  DimPlot(
    obj,
    reduction = "umap.rpca",
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE
  )
)

dev.off()


# ------------------------------------------------------------
# 10. UMAP by sample
# ------------------------------------------------------------

pdf(
  file.path(
    figure_dir,
    "UMAP_by_sample.pdf"
  ),
  width = 9,
  height = 7
)

print(
  DimPlot(
    obj,
    reduction = "umap.rpca",
    group.by = "sample_id"
  )
)

dev.off()


# ------------------------------------------------------------
# 11. UMAP by patient
# ------------------------------------------------------------

pdf(
  file.path(
    figure_dir,
    "UMAP_by_patient.pdf"
  ),
  width = 9,
  height = 7
)

print(
  DimPlot(
    obj,
    reduction = "umap.rpca",
    group.by = "patient_id"
  )
)

dev.off()


# ------------------------------------------------------------
# 12. Cluster composition by sample
# ------------------------------------------------------------

cluster_sample_table <- table(
  obj$seurat_clusters,
  obj$sample_id
)

write.csv(
  as.data.frame(cluster_sample_table),
  file.path(
    project_dir,
    "results",
    "clustering",
    "cluster_by_sample_counts.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. Save clustered object
# ------------------------------------------------------------

saveRDS(
  obj,
  output_file
)


# ------------------------------------------------------------
# 14. Final summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("CLUSTERING + UMAP COMPLETED\n")
cat("====================================================\n")

cat(
  "Cells:",
  ncol(obj),
  "\n"
)

cat(
  "Clusters:",
  length(
    unique(
      obj$seurat_clusters
    )
  ),
  "\n"
)

cat(
  "Dimensions used:",
  length(dims_use),
  "\n"
)

cat(
  "Resolution:",
  0.5,
  "\n"
)

cat("\nCluster sizes:\n")

print(
  table(
    obj$seurat_clusters
  )
)

cat(
  "\nSaved to:\n",
  output_file,
  "\n"
)