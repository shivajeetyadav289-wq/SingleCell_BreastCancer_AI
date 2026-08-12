# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 19_final_malignant_cell_labeling.R
#
# Purpose:
# Define CNV-supported malignant-cell candidates and transfer
# the annotation back to the full Seurat object.
#
# Primary definition:
#
#   Tumor_Candidate == TRUE
#   AND epithelial_evidence == TRUE
#   AND CopyKAT == "aneuploid"
#   AND sample is supported by sufficient aneuploid cells
#
# Primary supported samples:
#   BC11
#   BC12
#   BC17
#
# IMPORTANT:
# These are research annotations and are NOT clinical
# diagnostic calls.
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

suppressPackageStartupMessages({

  library(Seurat)
  library(dplyr)
  library(ggplot2)

})


# ------------------------------------------------------------
# 2. Project paths
# ------------------------------------------------------------

project_dir <- path.expand(
  "~/Projects/SingleCell_BreastCancer_AI"
)

seurat_file <- file.path(
  project_dir,
  "data",
  "processed",
  "GSE228499_malignant_candidates.rds"
)

integration_file <- file.path(
  project_dir,
  "results",
  "cnv",
  "integration",
  "GSE228499_CNV_marker_integrated_results.csv"
)

prediction_file <- file.path(
  project_dir,
  "results",
  "cnv",
  "copykat",
  "GSE228499_CopyKAT_prediction.rds"
)

result_dir <- file.path(
  project_dir,
  "results",
  "malignant"
)

figure_dir <- file.path(
  project_dir,
  "figures",
  "malignant"
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
# 3. Header
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("FINAL MALIGNANT-CELL LABELING\n")
cat("====================================================\n")


# ------------------------------------------------------------
# 4. Check input files
# ------------------------------------------------------------

required_files <- c(
  seurat_file,
  integration_file,
  prediction_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (
  length(missing_files) > 0
) {

  stop(
    "Missing required files:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )

}


# ------------------------------------------------------------
# 5. Load Seurat object
# ------------------------------------------------------------

cat("\n")
cat("Loading annotated Seurat object...\n")

obj <- readRDS(
  seurat_file
)

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
  "Assays:",
  paste(
    Assays(obj),
    collapse = ", "
  ),
  "\n"
)


# ------------------------------------------------------------
# 6. Load integration results
# ------------------------------------------------------------

cat("\n")
cat("Loading CNV + marker integration...\n")

integration <- read.csv(
  integration_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "Integration rows:",
  nrow(integration),
  "\n"
)


# ------------------------------------------------------------
# 7. Load CopyKAT prediction
# ------------------------------------------------------------

prediction <- readRDS(
  prediction_file
)

cat(
  "CopyKAT prediction rows:",
  nrow(prediction),
  "\n"
)

cat("\nCopyKAT prediction summary:\n")

print(
  table(
    prediction$copykat.pred,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# 8. Validate required Seurat metadata
# ------------------------------------------------------------

required_metadata <- c(
  "sample_id",
  "cell_type",
  "tumor_candidate"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(
    obj@meta.data
  )
)

if (
  length(missing_metadata) > 0
) {

  stop(
    "Missing Seurat metadata columns:\n",
    paste(
      missing_metadata,
      collapse = "\n"
    )
  )

}


# ------------------------------------------------------------
# 9. Normalize cell IDs
# ------------------------------------------------------------

normalize_cell_id <- function(x) {

  x <- as.character(x)

  x <- gsub(
    "\\.([0-9]+)$",
    "-\\1",
    x
  )

  x

}


seurat_ids <- colnames(
  obj
)

seurat_ids_normalized <- normalize_cell_id(
  seurat_ids
)

integration_ids <- normalize_cell_id(
  integration$cell_id
)

prediction_ids <- normalize_cell_id(
  prediction$cell.names
)


# ------------------------------------------------------------
# 10. Create CopyKAT lookup
# ------------------------------------------------------------

copykat_lookup <- data.frame(

  cell_id =
    prediction_ids,

  copykat.pred =
    prediction$copykat.pred,

  stringsAsFactors = FALSE

)


# ------------------------------------------------------------
# 11. Match CopyKAT predictions
# ------------------------------------------------------------

prediction_match <- match(
  seurat_ids_normalized,
  copykat_lookup$cell_id
)

obj$copykat_prediction <- NA_character_

matched_prediction <- !is.na(
  prediction_match
)

obj$copykat_prediction[
  matched_prediction
] <-

  copykat_lookup$copykat.pred[
    prediction_match[
      matched_prediction
    ]
  ]


cat("\n")
cat(
  "Seurat cells with CopyKAT prediction:",
  sum(matched_prediction),
  "\n"
)


# ------------------------------------------------------------
# 12. Normalize tumor-candidate field
# ------------------------------------------------------------

tumor_value <- tolower(
  trimws(
    as.character(
      obj$tumor_candidate
    )
  )
)

tumor_value <- gsub(
  "_",
  "",
  tumor_value
)

tumor_value <- gsub(
  " ",
  "",
  tumor_value
)

obj$tumor_candidate_flag <-
  tumor_value %in% c(
    "true",
    "1",
    "yes",
    "tumor",
    "candidate",
    "tumorcandidate"
  )


# ------------------------------------------------------------
# 13. Define epithelial evidence
# ------------------------------------------------------------

epithelial_cell_types <- c(
  "Basal_Epithelial",
  "Basal_Like_Epithelial",
  "Luminal_Epithelial",
  "Specialized_Epithelial"
)

obj$epithelial_evidence <-
  obj$cell_type %in%
  epithelial_cell_types


# ------------------------------------------------------------
# 14. Primary supported samples
# ------------------------------------------------------------

supported_samples <- c(
  "BC11",
  "BC12",
  "BC17"
)


# ------------------------------------------------------------
# 15. Define malignant candidate
# ------------------------------------------------------------

obj$malignant_candidate <-
  obj$tumor_candidate_flag &

  obj$epithelial_evidence &

  obj$copykat_prediction ==
    "aneuploid" &

  obj$sample_id %in%
    supported_samples


# ------------------------------------------------------------
# 16. Define broader CNV-supported category
# ------------------------------------------------------------

obj$cnv_supported_status <- "Not_supported"

obj$cnv_supported_status[
  obj$copykat_prediction ==
    "aneuploid"
] <- "CopyKAT_aneuploid"

obj$cnv_supported_status[
  obj$copykat_prediction ==
    "diploid"
] <- "CopyKAT_diploid"

obj$cnv_supported_status[
  obj$copykat_prediction ==
    "not.defined"
] <- "CopyKAT_not_defined"


# ------------------------------------------------------------
# 17. More descriptive malignant label
# ------------------------------------------------------------

obj$malignant_status <- "Non_malignant_or_unresolved"

obj$malignant_status[
  obj$malignant_candidate
] <- "CNV_supported_malignant_candidate"


# ------------------------------------------------------------
# 18. Summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("FINAL MALIGNANT CANDIDATE SUMMARY\n")
cat("====================================================\n")

cat("\nTumor-candidate flag:\n")

print(
  table(
    obj$tumor_candidate_flag,
    useNA = "ifany"
  )
)

cat("\nEpithelial evidence:\n")

print(
  table(
    obj$epithelial_evidence,
    useNA = "ifany"
  )
)

cat("\nCopyKAT prediction:\n")

print(
  table(
    obj$copykat_prediction,
    useNA = "ifany"
  )
)

cat("\nFinal malignant status:\n")

print(
  table(
    obj$malignant_status,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# 19. Malignant candidates by sample
# ------------------------------------------------------------

sample_summary <- obj@meta.data %>%

  mutate(
    cell_id =
      rownames(
        obj@meta.data
      )
  ) %>%

  group_by(
    sample_id
  ) %>%

  summarise(

    total_cells =
      n(),

    tumor_candidates =
      sum(
        tumor_candidate_flag,
        na.rm = TRUE
      ),

    epithelial_cells =
      sum(
        epithelial_evidence,
        na.rm = TRUE
      ),

    aneuploid_cells =
      sum(
        copykat_prediction ==
          "aneuploid",
        na.rm = TRUE
      ),

    malignant_candidates =
      sum(
        malignant_candidate,
        na.rm = TRUE
      ),

    malignant_percent =
      100 *
      malignant_candidates /
      total_cells,

    .groups = "drop"

  )


cat("\n")
cat("Sample-specific malignant candidates:\n")

print(
  sample_summary
)


write.csv(
  sample_summary,
  file.path(
    result_dir,
    "19_malignant_candidates_by_sample.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 20. Cell-type composition
# ------------------------------------------------------------

celltype_summary <- obj@meta.data %>%

  filter(
    malignant_candidate
  ) %>%

  count(
    cell_type,
    name = "malignant_cells"
  ) %>%

  arrange(
    desc(
      malignant_cells
    )
  )


cat("\n")
cat("Malignant candidate cell types:\n")

print(
  celltype_summary
)


write.csv(
  celltype_summary,
  file.path(
    result_dir,
    "19_malignant_candidates_by_cell_type.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 21. Create malignant subset
# ------------------------------------------------------------

malignant_cells <- colnames(
  obj
)[
  obj$malignant_candidate
]

cat("\n")
cat(
  "Final CNV-supported malignant candidates:",
  length(malignant_cells),
  "\n"
)


if (
  length(malignant_cells) == 0
) {

  stop(
    "No malignant candidates were identified."
  )

}


malignant_obj <- subset(
  obj,
  cells = malignant_cells
)


# ------------------------------------------------------------
# 22. Save malignant object
# ------------------------------------------------------------

malignant_file <- file.path(
  result_dir,
  "GSE228499_CNV_supported_malignant_candidates.rds"
)

saveRDS(
  malignant_obj,
  malignant_file
)

cat(
  "Saved malignant object:\n",
  malignant_file,
  "\n"
)


# ------------------------------------------------------------
# 23. Save full annotated object
# ------------------------------------------------------------

annotated_file <- file.path(
  project_dir,
  "data",
  "processed",
  "GSE228499_final_malignant_annotated.rds"
)

saveRDS(
  obj,
  annotated_file
)

cat(
  "Saved full annotated object:\n",
  annotated_file,
  "\n"
)


# ------------------------------------------------------------
# 24. Save malignant cell IDs
# ------------------------------------------------------------

malignant_cell_table <- data.frame(

  cell_id =
    malignant_cells,

  sample_id =
    obj$sample_id[
      malignant_cells
    ],

  cell_type =
    obj$cell_type[
      malignant_cells
    ],

  copykat_prediction =
    obj$copykat_prediction[
      malignant_cells
    ],

  tumor_candidate =
    obj$tumor_candidate_flag[
      malignant_cells
    ],

  epithelial_evidence =
    obj$epithelial_evidence[
      malignant_cells
    ],

  malignant_status =
    obj$malignant_status[
      malignant_cells
    ],

  stringsAsFactors = FALSE

)


write.csv(
  malignant_cell_table,
  file.path(
    result_dir,
    "19_malignant_candidate_cells.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 25. UMAP
# ------------------------------------------------------------

if (
  "umap" %in%
  Reductions(obj)
) {

  cat("\n")
  cat("Creating malignant-candidate UMAP...\n")

  pdf(
    file.path(
      figure_dir,
      "19_malignant_candidate_UMAP.pdf"
    ),
    width = 10,
    height = 8
  )

  print(

    DimPlot(
      obj,
      reduction = "umap",
      group.by = "malignant_status",
      raster = TRUE
    ) +

      ggtitle(
        "CNV-Supported Malignant Cell Candidates"
      )

  )

  dev.off()


  # ----------------------------------------------------------
  # UMAP by sample
  # ----------------------------------------------------------

  pdf(
    file.path(
      figure_dir,
      "19_malignant_candidate_UMAP_by_sample.pdf"
    ),
    width = 10,
    height = 8
  )

  print(

    DimPlot(
      malignant_obj,
      reduction = "umap",
      group.by = "sample_id",
      raster = TRUE
    ) +

      ggtitle(
        "CNV-Supported Malignant Candidates by Sample"
      )

  )

  dev.off()

} else {

  cat(
    "\nWARNING: UMAP reduction not present.\n"
  )

}


# ------------------------------------------------------------
# 26. Marker/cell-type composition plot
# ------------------------------------------------------------

if (
  nrow(
    celltype_summary
  ) > 0
) {

  pdf(
    file.path(
      figure_dir,
      "19_malignant_candidate_celltype_composition.pdf"
    ),
    width = 10,
    height = 6
  )

  print(

    ggplot(
      celltype_summary,
      aes(
        x = reorder(
          cell_type,
          malignant_cells
        ),
        y = malignant_cells
      )
    ) +

      geom_col() +

      coord_flip() +

      labs(
        title =
          "Cell-Type Composition of CNV-Supported Malignant Candidates",
        x =
          "Cell type",
        y =
          "Cells"
      ) +

      theme_classic()

  )

  dev.off()

}


# ------------------------------------------------------------
# 27. Final report
# ------------------------------------------------------------

report_file <- file.path(
  result_dir,
  "19_final_malignant_cell_labeling_summary.txt"
)

report <- c(

  "FINAL MALIGNANT CELL LABELING SUMMARY",

  "====================================",

  "",

  paste(
    "Total Seurat cells:",
    ncol(obj)
  ),

  paste(
    "CopyKAT aneuploid cells:",
    sum(
      obj$copykat_prediction ==
        "aneuploid",
      na.rm = TRUE
    )
  ),

  paste(
    "Final CNV-supported malignant candidates:",
    sum(
      obj$malignant_candidate,
      na.rm = TRUE
    )
  ),

  "",

  "Supported samples:",

  paste(
    supported_samples,
    collapse = ", "
  ),

  "",

  "Definition:",

  "Tumor-candidate epithelial cells with",

  "CopyKAT aneuploid status in samples with",

  "substantial aneuploid representation.",

  "",

  "Important:",

  "These are research annotations.",

  "They are not clinical diagnostic classifications."

)

writeLines(
  report,
  report_file
)


# ------------------------------------------------------------
# 28. Final summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("STEP 19 COMPLETED\n")
cat("====================================================\n")

cat(
  "Final malignant candidates:",
  sum(
    obj$malignant_candidate,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Malignant object:",
  malignant_file,
  "\n"
)

cat(
  "Full annotated object:",
  annotated_file,
  "\n"
)

cat(
  "Cell table:",
  file.path(
    result_dir,
    "19_malignant_candidate_cells.csv"
  ),
  "\n"
)

cat(
  "Summary:",
  report_file,
  "\n"
)

cat("\n")
cat("====================================================\n")
cat("READY FOR MALIGNANT-CELL CHARACTERIZATION\n")
cat("====================================================\n")