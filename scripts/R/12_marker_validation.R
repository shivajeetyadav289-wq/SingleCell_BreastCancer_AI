# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 12_marker_validation.R
#
# Purpose:
#   Validate computational clusters using canonical
#   cell-type marker genes.
#
# Input:
#   GSE228499_clustered.rds
#
# Outputs:
#   - Canonical marker DotPlot
#   - Canonical marker Heatmap
#   - Marker availability table
#   - Canonical marker panel
# ============================================================


# ============================================================
# 1. Load packages
# ============================================================

library(Seurat)
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

result_dir <- file.path(
  project_dir,
  "results",
  "marker_validation"
)

figure_dir <- file.path(
  project_dir,
  "figures",
  "marker_validation"
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
cat("CANONICAL MARKER VALIDATION\n")
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
# 4. Verify cluster information
# ============================================================

cat("\nMetadata columns:\n")

print(
  colnames(obj@meta.data)
)


if (!"seurat_clusters" %in% colnames(obj@meta.data)) {

  stop(
    "ERROR: seurat_clusters is not present in GSE228499_clustered.rds"
  )

}


# Use Seurat clusters as active identities

Idents(obj) <- "seurat_clusters"


cat("\nCluster sizes:\n")

print(
  table(
    obj$seurat_clusters
  )
)


# ============================================================
# 5. Set RNA assay
# ============================================================

DefaultAssay(obj) <- "RNA"

cat("\nDefault assay:\n")

cat(
  DefaultAssay(obj),
  "\n"
)


# ============================================================
# 6. Inspect RNA layers
# ============================================================

cat("\nRNA layers before preparation:\n")

print(
  Layers(
    obj[["RNA"]]
  )
)


# ============================================================
# 7. Join RNA layers ONLY if required
# ============================================================

rna_layers <- Layers(
  obj[["RNA"]]
)


if (
  any(grepl("^counts\\.", rna_layers)) ||
  any(grepl("^data\\.", rna_layers))
) {

  cat("\nJoining RNA layers...\n")

  obj <- JoinLayers(
    obj,
    assay = "RNA"
  )

} else {

  cat("\nRNA layers are already joined.\n")

}


cat("\nRNA layers after preparation:\n")

print(
  Layers(
    obj[["RNA"]]
  )
)


# ============================================================
# 8. Check normalized expression layer
# ============================================================

rna_layers <- Layers(
  obj[["RNA"]]
)


if (!"data" %in% rna_layers) {

  cat("\nNormalized RNA data layer not found.\n")

  cat("Running NormalizeData...\n")

  obj <- NormalizeData(
    obj,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

} else {

  cat("\nNormalized RNA data layer detected.\n")

}


# ============================================================
# 9. Define canonical marker panel
# ============================================================

marker_list <- list(

  Epithelial = c(
    "EPCAM",
    "KRT8",
    "KRT18",
    "KRT19",
    "KRT7"
  ),

  Luminal_Epithelial = c(
    "ESR1",
    "PGR",
    "GATA3",
    "FOXA1"
  ),

  Basal_Epithelial = c(
    "KRT5",
    "KRT14",
    "KRT17",
    "KRT6A",
    "KRT6B"
  ),

  T_cells = c(
    "CD3D",
    "CD3E",
    "TRAC",
    "IL7R",
    "LTB"
  ),

  Cytotoxic_NK = c(
    "NKG7",
    "GNLY",
    "GZMB",
    "GZMH",
    "CD8A",
    "CD8B"
  ),

  B_cells = c(
    "MS4A1",
    "CD79A",
    "CD37",
    "CD74",
    "HLA-DRA"
  ),

  Macrophages = c(
    "LYZ",
    "LST1",
    "C1QA",
    "C1QB",
    "C1QC",
    "TREM2",
    "FOLR2"
  ),

  Fibroblast_CAF = c(
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "DCN",
    "LUM",
    "COL6A1"
  ),

  Endothelial = c(
    "PECAM1",
    "VWF",
    "KDR",
    "EMCN",
    "CLDN5",
    "PLVAP"
  ),

  Cycling = c(
    "MKI67",
    "TOP2A",
    "UBE2C",
    "BIRC5",
    "FOXM1"
  ),

  Plasma_cells = c(
    "MZB1",
    "JCHAIN",
    "SDC1",
    "XBP1"
  ),

  Myeloid = c(
    "TYROBP",
    "FCER1G",
    "LILRB1",
    "AIF1"
  ),

  Pericytes = c(
    "RGS5",
    "CSPG4",
    "MCAM",
    "PDGFRB"
  ),

  Dendritic_cells = c(
    "FCER1A",
    "CD1C",
    "CLEC10A"
  )

)


# ============================================================
# 10. Check marker availability
# ============================================================

all_markers <- unique(
  unlist(marker_list)
)


genes_present <- intersect(
  all_markers,
  rownames(obj)
)


genes_missing <- setdiff(
  all_markers,
  rownames(obj)
)


cat("\n")
cat("====================================================\n")
cat("MARKER AVAILABILITY\n")
cat("====================================================\n")


cat(
  "Total unique markers:",
  length(all_markers),
  "\n"
)


cat(
  "Markers present:",
  length(genes_present),
  "\n"
)


cat(
  "Markers missing:",
  length(genes_missing),
  "\n"
)


# ============================================================
# 11. Save marker availability
# ============================================================

marker_availability <- data.frame(
  gene = all_markers,
  present = all_markers %in% rownames(obj),
  stringsAsFactors = FALSE
)


write.csv(
  marker_availability,
  file.path(
    result_dir,
    "marker_gene_availability.csv"
  ),
  row.names = FALSE
)


write.csv(
  data.frame(
    missing_gene = genes_missing
  ),
  file.path(
    result_dir,
    "missing_marker_genes.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 12. Keep only markers present in dataset
# ============================================================

marker_list_present <- lapply(
  marker_list,
  function(x) {

    x[
      x %in% rownames(obj)
    ]

  }
)


# Remove empty groups

marker_list_present <- marker_list_present[
  lengths(marker_list_present) > 0
]


# ============================================================
# 13. Create UNIQUE marker vector
# ============================================================

marker_genes <- unique(
  unlist(
    marker_list_present
  )
)


cat("\nMarkers used for validation:\n")

print(
  marker_genes
)


# ============================================================
# 14. Save canonical marker panel
# ============================================================

marker_panel <- do.call(
  rbind,
  lapply(
    names(marker_list_present),
    function(cell_type) {

      data.frame(
        cell_type = cell_type,
        gene = marker_list_present[[cell_type]],
        stringsAsFactors = FALSE
      )

    }
  )
)


write.csv(
  marker_panel,
  file.path(
    result_dir,
    "canonical_marker_panel.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 15. Canonical marker DotPlot
# ============================================================

cat("\n")
cat("Creating canonical marker DotPlot...\n")


dotplot <- DotPlot(
  obj,
  assay = "RNA",
  features = marker_genes,
  group.by = "seurat_clusters",
  dot.scale = 6
) +
  RotatedAxis() +
  ggtitle(
    "Canonical Cell-Type Marker Validation"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )


print(dotplot)


# ============================================================
# 16. Save DotPlot
# ============================================================

ggsave(
  filename = file.path(
    figure_dir,
    "canonical_marker_DotPlot.pdf"
  ),
  plot = dotplot,
  width = 18,
  height = 10
)


ggsave(
  filename = file.path(
    figure_dir,
    "canonical_marker_DotPlot.png"
  ),
  plot = dotplot,
  width = 18,
  height = 10,
  dpi = 300
)


# ============================================================
# 17. Scale canonical markers for heatmap
# ============================================================

cat("\n")
cat("Scaling canonical marker genes...\n")


obj <- ScaleData(
  obj,
  assay = "RNA",
  features = marker_genes,
  verbose = FALSE
)


# ============================================================
# 18. Canonical marker Heatmap
# ============================================================

cat("\nCreating canonical marker Heatmap...\n")


heatmap <- DoHeatmap(
  obj,
  assay = "RNA",
  features = marker_genes,
  group.by = "seurat_clusters",
  size = 3
) +
  ggtitle(
    "Canonical Cell-Type Marker Expression"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    )
  )


print(heatmap)


# ============================================================
# 19. Save Heatmap
# ============================================================

ggsave(
  filename = file.path(
    figure_dir,
    "canonical_marker_Heatmap.pdf"
  ),
  plot = heatmap,
  width = 18,
  height = 12
)


ggsave(
  filename = file.path(
    figure_dir,
    "canonical_marker_Heatmap.png"
  ),
  plot = heatmap,
  width = 18,
  height = 12,
  dpi = 300
)


# ============================================================
# 20. Save validation object
# ============================================================

saveRDS(
  obj,
  file.path(
    result_dir,
    "GSE228499_marker_validation.rds"
  )
)


# ============================================================
# 21. Final summary
# ============================================================

cat("\n")
cat("====================================================\n")
cat("MARKER VALIDATION COMPLETED\n")
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
  "Canonical markers used:",
  length(marker_genes),
  "\n"
)


cat("\nCluster sizes:\n")

print(
  table(
    obj$seurat_clusters
  )
)


cat("\nFigures saved:\n")

cat(
  file.path(
    figure_dir,
    "canonical_marker_DotPlot.png"
  ),
  "\n"
)

cat(
  file.path(
    figure_dir,
    "canonical_marker_Heatmap.png"
  ),
  "\n"
)


cat("\n====================================================\n")
cat("READY FOR BIOLOGICAL CLUSTER ANNOTATION\n")
cat("====================================================\n")