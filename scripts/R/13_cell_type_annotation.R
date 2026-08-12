# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 13_cell_type_annotation.R
#
# Purpose:
#   Annotate scRNA-seq clusters using canonical markers and
#   cluster-specific marker evidence.
# ============================================================


# ============================================================
# 1. Libraries
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)


# ============================================================
# 2. Project paths
# ============================================================

project_dir <- "~/Projects/SingleCell_BreastCancer_AI"

input_file <- file.path(
  project_dir,
  "data",
  "processed",
  "GSE228499_clustered.rds"
)

output_file <- file.path(
  project_dir,
  "data",
  "processed",
  "GSE228499_annotated.rds"
)

result_dir <- file.path(
  project_dir,
  "results",
  "annotation"
)

figure_dir <- file.path(
  project_dir,
  "figures",
  "annotation"
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


# ============================================================
# 3. Load clustered object
# ============================================================

cat("\n")
cat("====================================================\n")
cat("CELL-TYPE ANNOTATION\n")
cat("====================================================\n")

cat("\nLoading clustered dataset...\n")

obj <- readRDS(input_file)

cat(
  "Dimensions:",
  nrow(obj),
  "genes x",
  ncol(obj),
  "cells\n"
)


# ============================================================
# 4. Verify clusters
# ============================================================

if (
  !"seurat_clusters" %in%
  colnames(obj@meta.data)
) {

  stop(
    "ERROR: seurat_clusters is not present in the object."
  )

}

Idents(obj) <- "seurat_clusters"

cat("\nCluster sizes:\n")

print(
  table(obj$seurat_clusters)
)


# ============================================================
# 5. Check dimensional reductions
# ============================================================

cat("\n")
cat("Available dimensional reductions:\n")

print(
  names(obj@reductions)
)


# ============================================================
# 6. Create UMAP if missing
# ============================================================

if (
  !"umap" %in%
  names(obj@reductions)
) {

  cat("\n")
  cat("====================================================\n")
  cat("UMAP NOT FOUND\n")
  cat("Creating UMAP from integrated RPCA reduction...\n")
  cat("====================================================\n")

  if (
    !"integrated.rpca" %in%
    names(obj@reductions)
  ) {

    stop(
      paste(
        "ERROR: Neither UMAP nor integrated.rpca",
        "reduction is available."
      )
    )

  }


  obj <- RunUMAP(
    object = obj,
    reduction = "integrated.rpca",
    dims = 1:30,
    reduction.name = "umap",
    reduction.key = "UMAP_"
  )

  cat("\nUMAP successfully created.\n")

} else {

  cat("\n")
  cat("Existing UMAP found. Using existing UMAP.\n")

}


# ============================================================
# 7. Verify UMAP
# ============================================================

cat("\nDimensional reductions after UMAP check:\n")

print(
  names(obj@reductions)
)


# ============================================================
# 8. Cluster annotations
# ============================================================

cluster_annotation <- c(

  "0"  = "Luminal_Epithelial",

  "1"  = "Luminal_Epithelial",

  "2"  = "Basal_Like_Epithelial",

  "3"  = "Luminal_Epithelial",

  "4"  = "T_Cells",

  "5"  = "Cytotoxic_NK_T",

  "6"  = "Cycling_Proliferating",

  "7"  = "Fibroblast_CAF",

  "8"  = "Specialized_Epithelial",

  "9"  = "Macrophages",

  "10" = "Endothelial",

  "11" = "B_Cells",

  "12" = "Basal_Epithelial",

  "13" = "Luminal_Epithelial"

)


# ============================================================
# 9. Verify all clusters have annotations
# ============================================================

clusters_present <- sort(
  unique(
    as.character(
      obj$seurat_clusters
    )
  )
)

missing_clusters <- setdiff(
  clusters_present,
  names(cluster_annotation)
)

if (
  length(missing_clusters) > 0
) {

  stop(
    paste(
      "Missing annotations for clusters:",
      paste(
        missing_clusters,
        collapse = ", "
      )
    )
  )

}


# ============================================================
# 10. Add cell type
# ============================================================

obj$cell_type <- unname(
  cluster_annotation[
    as.character(
      obj$seurat_clusters
    )
  ]
)


# ============================================================
# 11. Annotation confidence
# ============================================================

annotation_confidence <- c(

  "0"  = "Medium",
  "1"  = "Medium",
  "2"  = "High",
  "3"  = "High",
  "4"  = "High",
  "5"  = "Very_High",
  "6"  = "Very_High",
  "7"  = "Very_High",
  "8"  = "Medium",
  "9"  = "Very_High",
  "10" = "Very_High",
  "11" = "Very_High",
  "12" = "Very_High",
  "13" = "High"

)

obj$annotation_confidence <- unname(
  annotation_confidence[
    as.character(
      obj$seurat_clusters
    )
  ]
)


# ============================================================
# 12. Annotation evidence table
# ============================================================

annotation_evidence <- data.frame(

  cluster = 0:13,

  cell_type = c(

    "Luminal_Epithelial",
    "Luminal_Epithelial",
    "Basal_Like_Epithelial",
    "Luminal_Epithelial",
    "T_Cells",
    "Cytotoxic_NK_T",
    "Cycling_Proliferating",
    "Fibroblast_CAF",
    "Specialized_Epithelial",
    "Macrophages",
    "Endothelial",
    "B_Cells",
    "Basal_Epithelial",
    "Luminal_Epithelial"

  ),

  confidence = c(

    "Medium",
    "Medium",
    "High",
    "High",
    "High",
    "Very_High",
    "Very_High",
    "Very_High",
    "Medium",
    "Very_High",
    "Very_High",
    "Very_High",
    "Very_High",
    "High"

  ),

  key_markers = c(

    "EPCAM, KRT8, KRT18, KRT19, GATA3",

    "EPCAM, KRT8, KRT18, KRT19, NQO1, FBP1",

    "SPNS2, CA9, KRT80, S100A2, LAMC2",

    "ANKRD30A, VMP1, ELF3, CDK5RAP3",

    "TRAC, CD3D, CD3E, IL7R, LTB, CD2",

    "NKG7, GNLY, GZMB, GZMH, CD8A, CD8B, KLRK1",

    "MKI67, TOP2A, UBE2C, BIRC5, PBK, FOXM1",

    "COL14A1, COL3A1, COL1A2, SFRP2, SFRP4, MFAP4",

    "LACRT, CXCL14, CXCL13, CA2, MUC5B",

    "FOLR2, C1QA, C1QB, C1QC, TREM2, MSR1, CSF1R",

    "EMCN, CLEC14A, CDH5, PLVAP, ESM1, CLDN5",

    "MS4A1, CD79A, CD37, HLA-DRA, HLA-DRB1, IGKC",

    "KLK5, KRT5, KRT14, KRT6B, GABRP, MMP7",

    "APOD, CALML5, MUCL1, SPINK8, IGFBP5, PIP"

  ),

  stringsAsFactors = FALSE
)


# ============================================================
# 13. Save annotation evidence
# ============================================================

write.csv(
  annotation_evidence,
  file.path(
    result_dir,
    "cluster_annotation_evidence.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 14. Print annotation table
# ============================================================

cat("\n")
cat("====================================================\n")
cat("CLUSTER ANNOTATIONS\n")
cat("====================================================\n")

print(
  annotation_evidence
)


# ============================================================
# 15. Cell-type counts
# ============================================================

cell_type_counts <- as.data.frame(
  table(
    obj$cell_type
  )
)

colnames(
  cell_type_counts
) <- c(
  "cell_type",
  "cell_count"
)

cell_type_counts$percentage <-
  cell_type_counts$cell_count /
  ncol(obj) *
  100

cell_type_counts <- cell_type_counts %>%
  arrange(
    desc(cell_count)
  )


write.csv(
  cell_type_counts,
  file.path(
    result_dir,
    "cell_type_counts.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. Print cell-type composition
# ============================================================

cat("\n")
cat("====================================================\n")
cat("CELL-TYPE COMPOSITION\n")
cat("====================================================\n")

print(
  cell_type_counts
)


# ============================================================
# 17. Cell-type composition plot
# ============================================================

composition_plot <- ggplot(
  cell_type_counts,
  aes(
    x = reorder(
      cell_type,
      cell_count
    ),
    y = cell_count
  )
) +

  geom_col() +

  coord_flip() +

  labs(
    title = "GSE228499 Cell-Type Composition",
    x = "Cell Type",
    y = "Number of Cells"
  ) +

  theme_classic() +

  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )


ggsave(
  file.path(
    figure_dir,
    "cell_type_composition.pdf"
  ),
  composition_plot,
  width = 10,
  height = 7
)

ggsave(
  file.path(
    figure_dir,
    "cell_type_composition.png"
  ),
  composition_plot,
  width = 10,
  height = 7,
  dpi = 300
)


# ============================================================
# 18. Annotated UMAP
# ============================================================

cat("\n")
cat("Creating annotated UMAP...\n")

umap_celltype <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "cell_type",
  label = TRUE,
  repel = TRUE,
  label.size = 4
) +

  ggtitle(
    "GSE228499 — Cell-Type Annotation"
  ) +

  theme_classic() +

  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    legend.title = element_blank()
  )


print(
  umap_celltype
)


# ============================================================
# 19. Save annotated UMAP
# ============================================================

ggsave(
  file.path(
    figure_dir,
    "UMAP_cell_type_annotation.pdf"
  ),
  umap_celltype,
  width = 12,
  height = 8
)

ggsave(
  file.path(
    figure_dir,
    "UMAP_cell_type_annotation.png"
  ),
  umap_celltype,
  width = 12,
  height = 8,
  dpi = 300
)


# ============================================================
# 20. Original cluster UMAP
# ============================================================

umap_cluster <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE
) +

  ggtitle(
    "GSE228499 — Seurat Clusters"
  ) +

  theme_classic()


ggsave(
  file.path(
    figure_dir,
    "UMAP_original_clusters.png"
  ),
  umap_cluster,
  width = 10,
  height = 8,
  dpi = 300
)


# ============================================================
# 21. Cell type by sample
# ============================================================

celltype_by_sample <- as.data.frame(
  table(
    obj$sample_id,
    obj$cell_type
  )
)

colnames(
  celltype_by_sample
) <- c(
  "sample_id",
  "cell_type",
  "cell_count"
)


write.csv(
  celltype_by_sample,
  file.path(
    result_dir,
    "cell_type_by_sample.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 22. Cell type percentage by sample
# ============================================================

celltype_percentage <- celltype_by_sample %>%

  group_by(
    sample_id
  ) %>%

  mutate(
    percentage =
      cell_count /
      sum(cell_count) *
      100
  ) %>%

  ungroup()


write.csv(
  celltype_percentage,
  file.path(
    result_dir,
    "cell_type_percentage_by_sample.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 23. Save annotated object
# ============================================================

saveRDS(
  obj,
  output_file
)


# ============================================================
# 24. Final verification
# ============================================================

cat("\n")
cat("====================================================\n")
cat("ANNOTATION COMPLETED\n")
cat("====================================================\n")

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
  "Cell types:",
  length(
    unique(
      obj$cell_type
    )
  ),
  "\n"
)

cat("\nCell types:\n")

print(
  table(
    obj$cell_type
  )
)

cat("\nDimensional reductions:\n")

print(
  names(obj@reductions)
)

cat("\nSaved object:\n")

cat(
  output_file,
  "\n"
)

cat("\n====================================================\n")
cat("READY FOR MALIGNANT CELL IDENTIFICATION\n")
cat("====================================================\n")