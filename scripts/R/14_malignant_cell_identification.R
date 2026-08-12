# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 14_malignant_cell_identification.R
#
# Purpose:
#   Identify epithelial and tumor-candidate populations
#   before CNV-based malignant-cell confirmation.
#
# IMPORTANT:
#   Marker expression alone is NOT used to definitively call
#   cells malignant. CNV inference will be used subsequently
#   for stronger malignant-cell identification.
# ============================================================


# ============================================================
# 1. Load packages
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
  "GSE228499_annotated.rds"
)

output_file <- file.path(
  project_dir,
  "data",
  "processed",
  "GSE228499_malignant_candidates.rds"
)

result_dir <- file.path(
  project_dir,
  "results",
  "malignant_identification"
)

figure_dir <- file.path(
  project_dir,
  "figures",
  "malignant_identification"
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
# 3. Load annotated object
# ============================================================

cat("\n")
cat("====================================================\n")
cat("MALIGNANT CELL IDENTIFICATION\n")
cat("====================================================\n")

cat("\nLoading annotated dataset...\n")

obj <- readRDS(input_file)

cat(
  "Dimensions:",
  nrow(obj),
  "genes x",
  ncol(obj),
  "cells\n"
)


# ============================================================
# 4. Verify required metadata
# ============================================================

required_metadata <- c(
  "seurat_clusters",
  "cell_type",
  "sample_id"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(obj@meta.data)
)

if (
  length(missing_metadata) > 0
) {

  stop(
    paste(
      "Missing required metadata:",
      paste(
        missing_metadata,
        collapse = ", "
      )
    )
  )

}


# ============================================================
# 5. Set identities
# ============================================================

Idents(obj) <- "seurat_clusters"


# ============================================================
# 6. Print cell-type composition
# ============================================================

cat("\n")
cat("====================================================\n")
cat("CURRENT CELL-TYPE COMPOSITION\n")
cat("====================================================\n")

print(
  table(
    obj$cell_type
  )
)


# ============================================================
# 7. Define canonical marker sets
# ============================================================

cat("\n")
cat("Defining canonical marker sets...\n")


# ------------------------------------------------------------
# Epithelial markers
# ------------------------------------------------------------

epithelial_markers <- c(

  "EPCAM",
  "KRT8",
  "KRT18",
  "KRT19",
  "KRT7",
  "KRT17",
  "KRT5",
  "KRT14",
  "KRT6A",
  "KRT6B",
  "KRT6C",
  "KRT4",
  "KRT13",
  "MUC1",
  "KRT23"

)


# ------------------------------------------------------------
# Luminal epithelial markers
# ------------------------------------------------------------

luminal_markers <- c(

  "EPCAM",
  "KRT8",
  "KRT18",
  "KRT19",
  "KRT7",
  "MUC1",
  "KRT23",
  "ESR1",
  "PGR",
  "GATA3",
  "FOXA1"

)


# ------------------------------------------------------------
# Basal epithelial markers
# ------------------------------------------------------------

basal_markers <- c(

  "KRT5",
  "KRT14",
  "KRT17",
  "KRT6A",
  "KRT6B",
  "KRT6C",
  "KRT23",
  "KRT16",
  "TP63",
  "SFN"

)


# ------------------------------------------------------------
# Immune markers
# ------------------------------------------------------------

immune_markers <- c(

  "PTPRC",
  "CD3D",
  "CD3E",
  "CD3G",
  "TRAC",
  "MS4A1",
  "CD79A",
  "CD74",
  "HLA-DRA",
  "NKG7",
  "GNLY",
  "LYZ",
  "FCER1G",
  "TYROBP"

)


# ------------------------------------------------------------
# Stromal markers
# ------------------------------------------------------------

stromal_markers <- c(

  "COL1A1",
  "COL1A2",
  "COL3A1",
  "COL6A1",
  "COL6A2",
  "COL6A3",
  "DCN",
  "LUM",
  "COL14A1",
  "SFRP2",
  "SFRP4",
  "PDGFRA"

)


# ------------------------------------------------------------
# Endothelial markers
# ------------------------------------------------------------

endothelial_markers <- c(

  "PECAM1",
  "VWF",
  "EMCN",
  "KDR",
  "ESM1",
  "CDH5",
  "CLDN5",
  "PLVAP",
  "ENG",
  "RAMP2"

)


# ============================================================
# 8. Keep genes actually present in dataset
# ============================================================

epithelial_markers <- intersect(
  epithelial_markers,
  rownames(obj)
)

luminal_markers <- intersect(
  luminal_markers,
  rownames(obj)
)

basal_markers <- intersect(
  basal_markers,
  rownames(obj)
)

immune_markers <- intersect(
  immune_markers,
  rownames(obj)
)

stromal_markers <- intersect(
  stromal_markers,
  rownames(obj)
)

endothelial_markers <- intersect(
  endothelial_markers,
  rownames(obj)
)


cat("\nMarkers available in dataset:\n")

cat(
  "Epithelial:",
  length(epithelial_markers),
  "\n"
)

cat(
  "Luminal:",
  length(luminal_markers),
  "\n"
)

cat(
  "Basal:",
  length(basal_markers),
  "\n"
)

cat(
  "Immune:",
  length(immune_markers),
  "\n"
)

cat(
  "Stromal:",
  length(stromal_markers),
  "\n"
)

cat(
  "Endothelial:",
  length(endothelial_markers),
  "\n"
)


# ============================================================
# 9. Calculate module scores
# ============================================================

cat("\n")
cat("====================================================\n")
cat("CALCULATING CELL-TYPE MODULE SCORES\n")
cat("====================================================\n")


obj <- AddModuleScore(
  object = obj,
  features = list(
    epithelial_markers,
    luminal_markers,
    basal_markers,
    immune_markers,
    stromal_markers,
    endothelial_markers
  ),
  name = "MarkerScore"
)


# ============================================================
# 10. Rename module-score columns
# ============================================================

score_columns <- c(
  "MarkerScore1",
  "MarkerScore2",
  "MarkerScore3",
  "MarkerScore4",
  "MarkerScore5",
  "MarkerScore6"
)

new_score_names <- c(
  "epithelial_score",
  "luminal_score",
  "basal_score",
  "immune_score",
  "stromal_score",
  "endothelial_score"
)

for (
  i in seq_along(score_columns)
) {

  if (
    score_columns[i] %in%
    colnames(obj@meta.data)
  ) {

    colnames(
      obj@meta.data
    )[
      colnames(
        obj@meta.data
      ) == score_columns[i]
    ] <- new_score_names[i]

  }

}


# ============================================================
# 11. Identify epithelial cells
# ============================================================

cat("\n")
cat("====================================================\n")
cat("IDENTIFYING EPITHELIAL CELLS\n")
cat("====================================================\n")


# Use annotated cell type as the primary biological definition
# and epithelial module score as supporting evidence.

epithelial_cell_types <- c(

  "Luminal_Epithelial",
  "Basal_Like_Epithelial",
  "Specialized_Epithelial",
  "Basal_Epithelial"

)


obj$epithelial_status <- ifelse(

  obj$cell_type %in%
    epithelial_cell_types,

  "Epithelial",

  "Non_Epithelial"

)


cat("\nEpithelial status:\n")

print(
  table(
    obj$epithelial_status
  )
)


# ============================================================
# 12. Identify epithelial clusters
# ============================================================

epithelial_clusters <- sort(
  unique(
    obj$seurat_clusters[
      obj$epithelial_status ==
        "Epithelial"
    ]
  )
)


cat("\nEpithelial clusters:\n")

print(
  epithelial_clusters
)


# ============================================================
# 13. Identify tumor candidate epithelial cells
# ============================================================

cat("\n")
cat("====================================================\n")
cat("IDENTIFYING TUMOR-CANDIDATE EPITHELIAL CELLS\n")
cat("====================================================\n")


# Tumor-candidate designation is deliberately conservative.
#
# These cells are epithelial populations that should be
# investigated further for malignant CNV patterns.
#
# This is NOT a final malignant-cell call.

tumor_candidate_clusters <- c(
  "0",
  "1",
  "2",
  "3",
  "8",
  "12",
  "13"
)


obj$tumor_candidate <- ifelse(

  obj$seurat_clusters %in%
    tumor_candidate_clusters,

  "Tumor_Candidate",

  "Other"

)


# ============================================================
# 14. Tumor candidate summary
# ============================================================

tumor_candidate_summary <- obj@meta.data %>%

  count(
    seurat_clusters,
    cell_type,
    tumor_candidate,
    name = "cell_count"
  ) %>%

  arrange(
    as.numeric(
      as.character(
        seurat_clusters
      )
    )
  )


write.csv(
  tumor_candidate_summary,
  file.path(
    result_dir,
    "tumor_candidate_summary.csv"
  ),
  row.names = FALSE
)


cat("\nTumor candidate summary:\n")

print(
  tumor_candidate_summary
)


# ============================================================
# 15. Epithelial marker DotPlot
# ============================================================

cat("\nCreating epithelial marker DotPlot...\n")


epithelial_dotplot_genes <- unique(
  c(
    epithelial_markers,
    luminal_markers,
    basal_markers
  )
)


pdf(
  file.path(
    figure_dir,
    "epithelial_marker_dotplot.pdf"
  ),
  width = 15,
  height = 9
)

print(

  DotPlot(
    obj,
    features = epithelial_dotplot_genes,
    group.by = "seurat_clusters"
  ) +

    RotatedAxis() +

    ggtitle(
      "Epithelial Marker Expression Across Clusters"
    ) +

    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      )
    )

)

dev.off()


# ============================================================
# 16. Marker score UMAP
# ============================================================

cat("\nCreating marker-score UMAPs...\n")


pdf(
  file.path(
    figure_dir,
    "epithelial_marker_score_umap.pdf"
  ),
  width = 12,
  height = 8
)

print(

  FeaturePlot(
    obj,
    features = c(
      "epithelial_score",
      "luminal_score",
      "basal_score"
    ),
    reduction = "umap",
    ncol = 3
  )

)

dev.off()


# ============================================================
# 17. Tumor candidate UMAP
# ============================================================

cat("\nCreating tumor-candidate UMAP...\n")


tumor_umap <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "tumor_candidate",
  label = FALSE
) +

  ggtitle(
    "Tumor-Candidate Epithelial Populations"
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
    "UMAP_tumor_candidate_cells.png"
  ),
  tumor_umap,
  width = 11,
  height = 8,
  dpi = 300
)

ggsave(
  file.path(
    figure_dir,
    "UMAP_tumor_candidate_cells.pdf"
  ),
  tumor_umap,
  width = 11,
  height = 8
)


# ============================================================
# 18. Epithelial-only object
# ============================================================

cat("\n")
cat("Creating epithelial-only object...\n")


epithelial_obj <- subset(
  obj,
  subset =
    epithelial_status ==
    "Epithelial"
)


cat(
  "Epithelial cells:",
  ncol(epithelial_obj),
  "\n"
)


cat("\nEpithelial cells by cluster:\n")

print(
  table(
    epithelial_obj$seurat_clusters
  )
)


# ============================================================
# 19. Save epithelial-only object
# ============================================================

epithelial_output <- file.path(

  project_dir,
  "data",
  "processed",
  "GSE228499_epithelial_cells.rds"

)


saveRDS(
  epithelial_obj,
  epithelial_output
)


# ============================================================
# 20. Tumor candidate object
# ============================================================

cat("\n")
cat("Creating tumor-candidate object...\n")


tumor_candidate_obj <- subset(
  obj,
  subset =
    tumor_candidate ==
    "Tumor_Candidate"
)


cat(
  "Tumor-candidate cells:",
  ncol(tumor_candidate_obj),
  "\n"
)


cat("\nTumor-candidate cells by cluster:\n")

print(
  table(
    tumor_candidate_obj$seurat_clusters
  )
)


# ============================================================
# 21. Save tumor candidate object
# ============================================================

tumor_candidate_output <- file.path(

  project_dir,
  "data",
  "processed",
  "GSE228499_tumor_candidate_epithelial.rds"

)


saveRDS(
  tumor_candidate_obj,
  tumor_candidate_output
)


# ============================================================
# 22. Export cell-level scoring table
# ============================================================

cell_score_table <- obj@meta.data %>%

  select(
    seurat_clusters,
    cell_type,
    sample_id,
    epithelial_status,
    tumor_candidate,
    epithelial_score,
    luminal_score,
    basal_score,
    immune_score,
    stromal_score,
    endothelial_score
  )


write.csv(
  cell_score_table,
  file.path(
    result_dir,
    "cell_type_marker_scores.csv"
  ),
  row.names = TRUE
)


# ============================================================
# 23. Cluster-level marker score summary
# ============================================================

cluster_score_summary <- obj@meta.data %>%

  group_by(
    seurat_clusters,
    cell_type
  ) %>%

  summarise(
    cells = n(),

    epithelial_score =
      median(
        epithelial_score,
        na.rm = TRUE
      ),

    luminal_score =
      median(
        luminal_score,
        na.rm = TRUE
      ),

    basal_score =
      median(
        basal_score,
        na.rm = TRUE
      ),

    immune_score =
      median(
        immune_score,
        na.rm = TRUE
      ),

    stromal_score =
      median(
        stromal_score,
        na.rm = TRUE
      ),

    endothelial_score =
      median(
        endothelial_score,
        na.rm = TRUE
      ),

    .groups = "drop"
  )


write.csv(
  cluster_score_summary,
  file.path(
    result_dir,
    "cluster_marker_score_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 24. Save annotated object with malignant candidates
# ============================================================

saveRDS(
  obj,
  output_file
)


# ============================================================
# 25. Final summary
# ============================================================

cat("\n")
cat("====================================================\n")
cat("MALIGNANT CELL PRE-SCREENING COMPLETED\n")
cat("====================================================\n")

cat(
  "Total cells:",
  ncol(obj),
  "\n"
)

cat(
  "Epithelial cells:",
  sum(
    obj$epithelial_status ==
      "Epithelial"
  ),
  "\n"
)

cat(
  "Tumor-candidate cells:",
  sum(
    obj$tumor_candidate ==
      "Tumor_Candidate"
  ),
  "\n"
)

cat(
  "Tumor-candidate percentage:",
  round(
    mean(
      obj$tumor_candidate ==
        "Tumor_Candidate"
    ) * 100,
    2
  ),
  "%\n"
)


cat("\nTumor-candidate clusters:\n")

print(
  tumor_candidate_clusters
)


cat("\nSaved annotated object:\n")

cat(
  output_file,
  "\n"
)

cat("\nSaved epithelial object:\n")

cat(
  epithelial_output,
  "\n"
)

cat("\nSaved tumor-candidate object:\n")

cat(
  tumor_candidate_output,
  "\n"
)


cat("\nFigures:\n")

cat(
  file.path(
    figure_dir,
    "epithelial_marker_dotplot.pdf"
  ),
  "\n"
)

cat(
  file.path(
    figure_dir,
    "epithelial_marker_score_umap.pdf"
  ),
  "\n"
)

cat(
  file.path(
    figure_dir,
    "UMAP_tumor_candidate_cells.pdf"
  ),
  "\n"
)


cat("\n")
cat("====================================================\n")
cat("NEXT STEP: CNV-BASED MALIGNANT CELL CONFIRMATION\n")
cat("====================================================\n")

cat(
  "Do NOT treat Tumor_Candidate as confirmed malignant cells yet.\n"
)

cat(
  "Use CNV inference to distinguish malignant epithelial cells\n"
)

cat(
  "from normal epithelial populations.\n"
)