# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 17_cnv_marker_integration.R
#
# Purpose:
# Integrate:
#   1. CopyKAT CNV prediction
#   2. Tumor-candidate status
#   3. Epithelial identity
#   4. Luminal/basal scores
#   5. Cell type annotation
#   6. Sample information
#
# Goal:
# Establish a conservative malignant-cell candidate population.
#
# IMPORTANT:
# Aneuploid != automatically malignant.
# Final classification requires CNV + epithelial/tumor evidence.
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

suppressPackageStartupMessages({

  library(dplyr)
  library(ggplot2)

})


# ------------------------------------------------------------
# 2. Project paths
# ------------------------------------------------------------

project_dir <- path.expand(
  "~/Projects/SingleCell_BreastCancer_AI"
)

validation_dir <- file.path(
  project_dir,
  "results",
  "cnv",
  "validation"
)

copykat_dir <- file.path(
  project_dir,
  "results",
  "cnv",
  "copykat"
)

result_dir <- file.path(
  project_dir,
  "results",
  "cnv",
  "integration"
)

figure_dir <- file.path(
  project_dir,
  "figures",
  "cnv",
  "integration"
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
# 3. Input files
# ------------------------------------------------------------

validation_file <- file.path(
  validation_dir,
  "GSE228499_CopyKAT_validation_table.csv"
)

prediction_file <- file.path(
  copykat_dir,
  "GSE228499_CopyKAT_prediction.rds"
)

cnv_file <- file.path(
  copykat_dir,
  "GSE228499_CopyKAT_CNAmat.rds"
)


# ------------------------------------------------------------
# 4. Header
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("CNV + MARKER INTEGRATION\n")
cat("====================================================\n")

cat(
  "Project: AI-Assisted Discovery of Breast Cancer Biomarkers\n"
)

cat(
  "Dataset: GSE228499\n"
)


# ------------------------------------------------------------
# 5. Check input files
# ------------------------------------------------------------

required_files <- c(
  validation_file,
  prediction_file,
  cnv_file
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
# 6. Load validation table
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("LOADING CNV VALIDATION TABLE\n")
cat("====================================================\n")

validation <- read.csv(
  validation_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "Rows:",
  nrow(validation),
  "\n"
)

cat(
  "Columns:",
  ncol(validation),
  "\n"
)

cat("\nAvailable columns:\n")

print(
  colnames(validation)
)


# ------------------------------------------------------------
# 7. Validate essential columns
# ------------------------------------------------------------

required_columns <- c(
  "cell.names",
  "cell_id",
  "copykat.pred",
  "sample_id"
)

missing_columns <- setdiff(
  required_columns,
  colnames(validation)
)

if (
  length(missing_columns) > 0
) {

  stop(
    "Required columns missing:\n",
    paste(
      missing_columns,
      collapse = ", "
    )
  )

}


# ------------------------------------------------------------
# 8. Normalize categorical fields
# ------------------------------------------------------------

validation$copykat.pred <-
  trimws(
    as.character(
      validation$copykat.pred
    )
  )

validation$sample_id <-
  trimws(
    as.character(
      validation$sample_id
    )
  )


# ------------------------------------------------------------
# 9. Tumor-candidate flag
# ------------------------------------------------------------
if (
  "tumor_candidate" %in%
    colnames(validation)
) {

  tumor_value <- tolower(
    trimws(
      as.character(
        validation$tumor_candidate
      )
    )
  )

  # Normalize separators so that:
  # Tumor_Candidate
  # tumor_candidate
  # Tumor Candidate
  # tumorcandidate
  # are treated consistently.

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

  validation$tumor_candidate_flag <-
    tumor_value %in% c(
      "true",
      "1",
      "yes",
      "tumor",
      "candidate",
      "tumorcandidate"
    )

} else {

  validation$tumor_candidate_flag <-
    FALSE

}
# ------------------------------------------------------------
# 10. Epithelial flag
# ------------------------------------------------------------

if (
  "epithelial_status" %in%
    colnames(validation)
) {

  epithelial_value <- tolower(
    trimws(
      as.character(
        validation$epithelial_status
      )
    )
  )

  validation$epithelial_flag <-
    epithelial_value %in% c(
      "true",
      "1",
      "yes",
      "epithelial",
      "positive"
    )

} else {

  validation$epithelial_flag <-
    FALSE

}


# ------------------------------------------------------------
# 11. Cell-type epithelial detection
# ------------------------------------------------------------

if (
  "cell_type" %in%
    colnames(validation)
) {

  cell_type_lower <- tolower(
    trimws(
      as.character(
        validation$cell_type
      )
    )
  )

  cell_type_epithelial <-
    grepl(
      "epithelial|luminal|basal|tumor|malignant",
      cell_type_lower
    )

} else {

  cell_type_epithelial <-
    FALSE

}


# ------------------------------------------------------------
# 12. Combine epithelial evidence
# ------------------------------------------------------------

validation$epithelial_evidence <-
  validation$epithelial_flag |
  cell_type_epithelial


# ------------------------------------------------------------
# 13. Marker score evidence
# ------------------------------------------------------------

validation$marker_score_evidence <-
  FALSE


if (
  "epithelial_score" %in%
    colnames(validation)
) {

  epithelial_score <- suppressWarnings(
    as.numeric(
      validation$epithelial_score
    )
  )

  # Do not impose an arbitrary biological threshold here.
  # We only use the existence of a positive epithelial score
  # as supporting evidence.

  validation$marker_score_evidence <-
    !is.na(epithelial_score) &
    epithelial_score > 0

}


# ------------------------------------------------------------
# 14. CNV evidence
# ------------------------------------------------------------

validation$cnv_aneuploid <-
  validation$copykat.pred ==
  "aneuploid"

validation$cnv_diploid <-
  validation$copykat.pred ==
  "diploid"

validation$cnv_not_defined <-
  validation$copykat.pred ==
  "not.defined"


# ------------------------------------------------------------
# 15. Evidence score
# ------------------------------------------------------------

# Evidence components:
#
# CNV aneuploid       = 2 points
# Tumor candidate     = 2 points
# Epithelial identity = 1 point
# Positive epithelial
# marker score        = 1 point
#
# Maximum = 6
#
# Classification:
#
# 6 = High-confidence malignant candidate
# 5 = High-confidence malignant candidate
# 4 = Malignant candidate
# 3 = Possible malignant candidate
# 0-2 = Not supported as malignant
#
# IMPORTANT:
# This is an evidence-integration framework, not a clinical
# diagnostic classifier.


validation$malignant_evidence_score <- 0


validation$malignant_evidence_score <-
  validation$malignant_evidence_score +
  ifelse(
    validation$cnv_aneuploid,
    2,
    0
  )


validation$malignant_evidence_score <-
  validation$malignant_evidence_score +
  ifelse(
    validation$tumor_candidate_flag,
    2,
    0
  )


validation$malignant_evidence_score <-
  validation$malignant_evidence_score +
  ifelse(
    validation$epithelial_evidence,
    1,
    0
  )


validation$malignant_evidence_score <-
  validation$malignant_evidence_score +
  ifelse(
    validation$marker_score_evidence,
    1,
    0
  )


# ------------------------------------------------------------
# 16. Final integrated classification
# ------------------------------------------------------------

validation$malignant_class <- "Not supported"


validation$malignant_class[
  validation$malignant_evidence_score >= 3
] <- "Possible malignant candidate"


validation$malignant_class[
  validation$malignant_evidence_score >= 4
] <- "Malignant candidate"


validation$malignant_class[
  validation$malignant_evidence_score >= 5
] <- "High-confidence malignant candidate"


# ------------------------------------------------------------
# 17. Protect against non-aneuploid malignant calls
# ------------------------------------------------------------

# A cell cannot receive a malignant classification solely
# from marker evidence if CopyKAT is diploid/not-defined.

validation$malignant_class[
  !validation$cnv_aneuploid
] <- "Not supported"


# ------------------------------------------------------------
# 18. Final malignant candidate flag
# ------------------------------------------------------------

validation$final_malignant_candidate <-
  validation$malignant_class %in% c(
    "Malignant candidate",
    "High-confidence malignant candidate"
  )


# ------------------------------------------------------------
# 19. Summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("INTEGRATED EVIDENCE SUMMARY\n")
cat("====================================================\n")

cat("\nCopyKAT:\n")

print(
  table(
    validation$copykat.pred,
    useNA = "ifany"
  )
)


cat("\nTumor candidate:\n")

print(
  table(
    validation$tumor_candidate_flag,
    useNA = "ifany"
  )
)


cat("\nEpithelial evidence:\n")

print(
  table(
    validation$epithelial_evidence,
    useNA = "ifany"
  )
)


cat("\nIntegrated malignant classification:\n")

print(
  table(
    validation$malignant_class,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# 20. Final malignant candidates
# ------------------------------------------------------------

malignant_candidates <- validation %>%

  filter(
    final_malignant_candidate
  )


cat("\n")
cat(
  "Final malignant candidates:",
  nrow(malignant_candidates),
  "\n"
)


# ------------------------------------------------------------
# 21. High-confidence malignant candidates
# ------------------------------------------------------------

high_confidence <- validation %>%

  filter(
    malignant_class ==
      "High-confidence malignant candidate"
  )


cat(
  "High-confidence malignant candidates:",
  nrow(high_confidence),
  "\n"
)


# ------------------------------------------------------------
# 22. Sample distribution
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("MALIGNANT CANDIDATES BY SAMPLE\n")
cat("====================================================\n")

malignant_by_sample <- validation %>%

  group_by(
    sample_id
  ) %>%

  summarise(

    total_cells = n(),

    aneuploid =
      sum(
        cnv_aneuploid,
        na.rm = TRUE
      ),

    malignant_candidates =
      sum(
        final_malignant_candidate,
        na.rm = TRUE
      ),

    high_confidence =
      sum(
        malignant_class ==
          "High-confidence malignant candidate",
        na.rm = TRUE
      ),

    malignant_percent =
      100 *
      malignant_candidates /
      total_cells,

    .groups = "drop"

  )


print(
  malignant_by_sample
)


# ------------------------------------------------------------
# 23. Save sample summary
# ------------------------------------------------------------

write.csv(
  malignant_by_sample,
  file.path(
    result_dir,
    "malignant_candidates_by_sample.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 24. Save complete integration table
# ------------------------------------------------------------

integration_file <- file.path(
  result_dir,
  "GSE228499_CNV_marker_integrated_results.csv"
)

write.csv(
  validation,
  integration_file,
  row.names = FALSE
)


# ------------------------------------------------------------
# 25. Save malignant candidates
# ------------------------------------------------------------

malignant_file <- file.path(
  result_dir,
  "GSE228499_final_malignant_candidates.csv"
)

write.csv(
  malignant_candidates,
  malignant_file,
  row.names = FALSE
)


# ------------------------------------------------------------
# 26. Save high-confidence candidates
# ------------------------------------------------------------

high_confidence_file <- file.path(
  result_dir,
  "GSE228499_high_confidence_malignant_candidates.csv"
)

write.csv(
  high_confidence,
  high_confidence_file,
  row.names = FALSE
)


# ------------------------------------------------------------
# 27. Evidence score distribution
# ------------------------------------------------------------

score_summary <- validation %>%

  count(
    malignant_evidence_score
  ) %>%

  arrange(
    malignant_evidence_score
  )


write.csv(
  score_summary,
  file.path(
    result_dir,
    "malignant_evidence_score_distribution.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 28. Plot: CopyKAT prediction
# ------------------------------------------------------------

prediction_plot <- validation %>%

  count(
    copykat.pred
  )


pdf(
  file.path(
    figure_dir,
    "CNV_prediction_summary.pdf"
  ),
  width = 8,
  height = 6
)

print(

  ggplot(
    prediction_plot,
    aes(
      x = copykat.pred,
      y = n
    )
  ) +

    geom_col() +

    geom_text(
      aes(
        label = n
      ),
      vjust = -0.4
    ) +

    labs(
      title = "CopyKAT CNV Prediction",
      x = "CopyKAT prediction",
      y = "Cells"
    ) +

    theme_classic()

)

dev.off()


# ------------------------------------------------------------
# 29. Plot: Malignant classification
# ------------------------------------------------------------

classification_plot <- validation %>%

  count(
    malignant_class
  )


pdf(
  file.path(
    figure_dir,
    "CNV_marker_malignant_classification.pdf"
  ),
  width = 10,
  height = 6
)

print(

  ggplot(
    classification_plot,
    aes(
      x = malignant_class,
      y = n
    )
  ) +

    geom_col() +

    geom_text(
      aes(
        label = n
      ),
      vjust = -0.4
    ) +

    labs(
      title =
        "Integrated CNV + Marker Malignant Classification",
      x = "Classification",
      y = "Cells"
    ) +

    theme_classic() +

    theme(
      axis.text.x =
        element_text(
          angle = 25,
          hjust = 1
        )
    )

)

dev.off()


# ------------------------------------------------------------
# 30. Plot: malignant candidates by sample
# ------------------------------------------------------------

pdf(
  file.path(
    figure_dir,
    "Malignant_candidates_by_sample.pdf"
  ),
  width = 10,
  height = 6
)

print(

  ggplot(
    malignant_by_sample,
    aes(
      x = sample_id,
      y = malignant_percent
    )
  ) +

    geom_col() +

    geom_text(
      aes(
        label =
          paste0(
            round(
              malignant_percent,
              1
            ),
            "%"
          )
      ),
      vjust = -0.4
    ) +

    labs(
      title =
        "Integrated Malignant Candidate Fraction by Sample",
      x = "Sample",
      y = "Malignant candidates (%)"
    ) +

    theme_classic()

)

dev.off()


# ------------------------------------------------------------
# 31. Plot: evidence score
# ------------------------------------------------------------

pdf(
  file.path(
    figure_dir,
    "Malignant_evidence_score.pdf"
  ),
  width = 8,
  height = 6
)

print(

  ggplot(
    score_summary,
    aes(
      x = factor(
        malignant_evidence_score
      ),
      y = n
    )
  ) +

    geom_col() +

    geom_text(
      aes(
        label = n
      ),
      vjust = -0.4
    ) +

    labs(
      title =
        "Integrated Malignant Evidence Score",
      x = "Evidence score",
      y = "Cells"
    ) +

    theme_classic()

)

dev.off()


# ------------------------------------------------------------
# 32. Create summary report
# ------------------------------------------------------------

report_file <- file.path(
  result_dir,
  "CNV_marker_integration_summary.txt"
)


report_lines <- c(

  "CNV + MARKER INTEGRATION SUMMARY",

  "================================",

  "",

  paste(
    "Total CopyKAT predictions:",
    nrow(validation)
  ),

  paste(
    "Aneuploid:",
    sum(
      validation$cnv_aneuploid,
      na.rm = TRUE
    )
  ),

  paste(
    "Diploid:",
    sum(
      validation$cnv_diploid,
      na.rm = TRUE
    )
  ),

  paste(
    "Not defined:",
    sum(
      validation$cnv_not_defined,
      na.rm = TRUE
    )
  ),

  "",

  paste(
    "Final malignant candidates:",
    nrow(malignant_candidates)
  ),

  paste(
    "High-confidence malignant candidates:",
    nrow(high_confidence)
  ),

  "",

  "Interpretation:",

  "Aneuploid status provides CNV evidence but does not",

  "by itself establish malignancy.",

  "Final candidates require integration with tumor-candidate",

  "and epithelial/marker evidence.",

  "These classifications are research annotations and",

  "not clinical diagnostic calls."

)


writeLines(
  report_lines,
  report_file
)


# ------------------------------------------------------------
# 33. Final output
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("CNV + MARKER INTEGRATION COMPLETED\n")
cat("====================================================\n")

cat(
  "Integrated table:\n",
  integration_file,
  "\n",
  sep = ""
)

cat(
  "Final malignant candidates:\n",
  malignant_file,
  "\n",
  sep = ""
)

cat(
  "High-confidence candidates:\n",
  high_confidence_file,
  "\n",
  sep = ""
)

cat(
  "Summary report:\n",
  report_file,
  "\n",
  sep = ""
)

cat("\nFigures:\n")

cat(
  figure_dir,
  "\n"
)

cat("\n")
cat("====================================================\n")
cat("MALIGNANT POPULATION VALIDATION COMPLETED\n")
cat("====================================================\n")

cat(
  "Next step: transfer validated malignant labels back\n",
  "to the Seurat object for downstream biomarker analysis.\n",
  sep = ""
)