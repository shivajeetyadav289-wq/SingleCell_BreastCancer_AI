# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 16_copykat_validation.R
#
# Purpose:
# Validate CopyKAT CNV predictions against epithelial/tumor
# candidate identity and sample information.
#
# CopyKAT input:
#   4,000 representative tumor-candidate cells
#
# CopyKAT predictions:
#   diploid
#   aneuploid
#   not.defined
#
# IMPORTANT:
# This script DOES NOT rerun CopyKAT.
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

suppressPackageStartupMessages({

  library(dplyr)
  library(ggplot2)
  library(readr)

})


# ------------------------------------------------------------
# 2. Project paths
# ------------------------------------------------------------

project_dir <- path.expand(
  "~/Projects/SingleCell_BreastCancer_AI"
)

cnv_dir <- file.path(
  project_dir,
  "results",
  "cnv",
  "copykat"
)

subset_dir <- file.path(
  project_dir,
  "results",
  "cnv",
  "subset"
)

result_dir <- file.path(
  project_dir,
  "results",
  "cnv",
  "validation"
)

figure_dir <- file.path(
  project_dir,
  "figures",
  "cnv",
  "validation"
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

prediction_file <- file.path(
  cnv_dir,
  "GSE228499_CopyKAT_prediction.rds"
)

cnv_file <- file.path(
  cnv_dir,
  "GSE228499_CopyKAT_CNAmat.rds"
)

metadata_file <- file.path(
  subset_dir,
  "GSE228499_CopyKAT_subset_metadata.csv"
)

selected_cells_file <- file.path(
  subset_dir,
  "GSE228499_CopyKAT_selected_cells.csv"
)


# ------------------------------------------------------------
# 4. Header
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("COPYKAT CNV VALIDATION\n")
cat("====================================================\n")


# ------------------------------------------------------------
# 5. Check files
# ------------------------------------------------------------

required_files <- c(
  prediction_file,
  cnv_file,
  metadata_file,
  selected_cells_file
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
# 6. Load CopyKAT predictions
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("LOADING COPYKAT PREDICTIONS\n")
cat("====================================================\n")

prediction <- readRDS(
  prediction_file
)

cat(
  "Prediction dimensions:",
  nrow(prediction),
  "rows x",
  ncol(prediction),
  "columns\n"
)

cat(
  "Prediction columns:\n"
)

print(
  colnames(prediction)
)


# ------------------------------------------------------------
# 7. Validate prediction structure
# ------------------------------------------------------------

required_prediction_columns <- c(
  "cell.names",
  "copykat.pred"
)

if (
  !all(
    required_prediction_columns %in%
      colnames(prediction)
  )
) {

  stop(
    "Expected CopyKAT prediction columns were not found."
  )

}


# ------------------------------------------------------------
# 8. Prediction summary
# ------------------------------------------------------------

cat("\n")
cat("CopyKAT prediction summary:\n")

print(
  table(
    prediction$copykat.pred,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# 9. Load metadata
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("LOADING CELL METADATA\n")
cat("====================================================\n")

metadata <- read.csv(
  metadata_file,
  row.names = 1,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

metadata$cell_id <- rownames(
  metadata
)

cat(
  "Metadata dimensions:",
  nrow(metadata),
  "rows x",
  ncol(metadata),
  "columns\n"
)


# ------------------------------------------------------------
# 10. Load selected cell list
# ------------------------------------------------------------

selected_cells <- read.csv(
  selected_cells_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "Selected cell records:",
  nrow(selected_cells),
  "\n"
)


# ------------------------------------------------------------
# 11. Validate cell IDs
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("MATCHING COPYKAT CELLS TO METADATA\n")
cat("====================================================\n")


# CopyKAT stores the original cell names.

prediction$cell_id <- prediction$cell.names


# Direct matching first

match_index <- match(
  prediction$cell_id,
  metadata$cell_id
)


# ------------------------------------------------------------
# 12. Handle possible cell-name formatting
# ------------------------------------------------------------

# CopyKAT/R may convert:
#
# BC03_xxxxx-1
#
# into:
#
# BC03_xxxxx.1
#
# for CNV matrix column names.
#
# Prediction cell.names should normally retain the original
# names, but we keep a robust matching procedure.


unmatched <- is.na(
  match_index
)


if (
  any(unmatched)
) {

  cat(
    "Direct metadata matches:",
    sum(!unmatched),
    "\n"
  )

  cat(
    "Unmatched CopyKAT cells:",
    sum(unmatched),
    "\n"
  )

  # Try removing possible R name conversion

  prediction_clean <- gsub(
    "\\.",
    "-",
    prediction$cell_id
  )

  metadata_clean <- gsub(
    "\\.",
    "-",
    metadata$cell_id
  )

  alternative_match <- match(
    prediction_clean,
    metadata_clean
  )

  match_index[unmatched] <-
    alternative_match[unmatched]

}


if (
  any(
    is.na(match_index)
  )
) {

  cat(
    "WARNING: Some CopyKAT cells could not be matched.\n"
  )

  cat(
    "Unmatched:",
    sum(
      is.na(match_index)
    ),
    "\n"
  )

}


# ------------------------------------------------------------
# 13. Build validation table
# ------------------------------------------------------------

validation <- prediction

validation$sample_id <- NA_character_

validation$patient_id <- NA_character_

validation$cancer_subtype <- NA_character_

validation$treatment <- NA_character_

validation$cell_type <- NA_character_

validation$annotation_confidence <- NA_character_

validation$epithelial_status <- NA_character_

validation$tumor_candidate <- NA_character_

validation$epithelial_score <- NA_real_

validation$luminal_score <- NA_real_

validation$basal_score <- NA_real_


matched <- !is.na(
  match_index
)


validation$sample_id[
  matched
] <-
  metadata$sample_id[
    match_index[matched]
  ]

validation$patient_id[
  matched
] <-
  metadata$patient_id[
    match_index[matched]
  ]

validation$cancer_subtype[
  matched
] <-
  metadata$cancer_subtype[
    match_index[matched]
  ]

validation$treatment[
  matched
] <-
  metadata$treatment[
    match_index[matched]
  ]


# ------------------------------------------------------------
# 14. Optional annotation columns
# ------------------------------------------------------------

if (
  "cell_type" %in%
    colnames(metadata)
) {

  validation$cell_type[
    matched
  ] <-
    metadata$cell_type[
      match_index[matched]
    ]

}


if (
  "annotation_confidence" %in%
    colnames(metadata)
) {

  validation$annotation_confidence[
    matched
  ] <-
    metadata$annotation_confidence[
      match_index[matched]
    ]

}


if (
  "epithelial_status" %in%
    colnames(metadata)
) {

  validation$epithelial_status[
    matched
  ] <-
    metadata$epithelial_status[
      match_index[matched]
    ]

}


if (
  "tumor_candidate" %in%
    colnames(metadata)
) {

  validation$tumor_candidate[
    matched
  ] <-
    metadata$tumor_candidate[
      match_index[matched]
    ]

}


if (
  "epithelial_score" %in%
    colnames(metadata)
) {

  validation$epithelial_score[
    matched
  ] <-
    metadata$epithelial_score[
      match_index[matched]
    ]

}


if (
  "luminal_score" %in%
    colnames(metadata)
) {

  validation$luminal_score[
    matched
  ] <-
    metadata$luminal_score[
      match_index[matched]
    ]

}


if (
  "basal_score" %in%
    colnames(metadata)
) {

  validation$basal_score[
    matched
  ] <-
    metadata$basal_score[
      match_index[matched]
    ]

}


# ------------------------------------------------------------
# 15. Matching summary
# ------------------------------------------------------------

cat("\n")
cat("Metadata matching summary:\n")

cat(
  "CopyKAT cells:",
  nrow(validation),
  "\n"
)

cat(
  "Matched:",
  sum(matched),
  "\n"
)

cat(
  "Unmatched:",
  sum(!matched),
  "\n"
)


if (
  sum(matched) == 0
) {

  stop(
    "No CopyKAT cells could be matched to metadata."
  )

}


# ------------------------------------------------------------
# 16. CopyKAT by sample
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("COPYKAT PREDICTIONS BY SAMPLE\n")
cat("====================================================\n")

sample_prediction <- validation %>%

  filter(
    !is.na(sample_id)
  ) %>%

  count(
    sample_id,
    copykat.pred
  ) %>%

  group_by(
    sample_id
  ) %>%

  mutate(
    percentage =
      100 * n / sum(n)
  ) %>%

  ungroup()


print(
  sample_prediction
)


# ------------------------------------------------------------
# 17. Save sample summary
# ------------------------------------------------------------

write.csv(
  sample_prediction,
  file.path(
    result_dir,
    "CopyKAT_prediction_by_sample.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 18. Aneuploid summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("ANEUPLOID CELL SUMMARY\n")
cat("====================================================\n")

aneuploid <- validation %>%

  filter(
    copykat.pred == "aneuploid"
  )


cat(
  "Aneuploid cells:",
  nrow(aneuploid),
  "\n"
)


# ------------------------------------------------------------
# 19. Diploid summary
# ------------------------------------------------------------

diploid <- validation %>%

  filter(
    copykat.pred == "diploid"
  )


cat(
  "Diploid cells:",
  nrow(diploid),
  "\n"
)


# ------------------------------------------------------------
# 20. Not-defined summary
# ------------------------------------------------------------

not_defined <- validation %>%

  filter(
    copykat.pred == "not.defined"
  )


cat(
  "Not-defined cells:",
  nrow(not_defined),
  "\n"
)


# ------------------------------------------------------------
# 21. Aneuploid cells by sample
# ------------------------------------------------------------

aneuploid_by_sample <- validation %>%

  filter(
    !is.na(sample_id)
  ) %>%

  group_by(
    sample_id
  ) %>%

  summarise(

    total_cells = n(),

    aneuploid =
      sum(
        copykat.pred ==
          "aneuploid"
      ),

    diploid =
      sum(
        copykat.pred ==
          "diploid"
      ),

    not_defined =
      sum(
        copykat.pred ==
          "not.defined"
      ),

    aneuploid_percent =
      100 *
      aneuploid /
      total_cells,

    .groups = "drop"

  )


cat("\nAneuploid cells by sample:\n")

print(
  aneuploid_by_sample
)


write.csv(
  aneuploid_by_sample,
  file.path(
    result_dir,
    "CopyKAT_aneuploid_by_sample.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 22. Validate tumor-candidate status
# ------------------------------------------------------------

if (
  "tumor_candidate" %in%
    colnames(validation)
) {

  cat("\n")
  cat("====================================================\n")
  cat("TUMOR-CANDIDATE VALIDATION\n")
  cat("====================================================\n")

  cat(
    "Tumor-candidate distribution:\n"
  )

  print(
    table(
      validation$tumor_candidate,
      useNA = "ifany"
    )
  )


  cat(
    "\nAneuploid cells by tumor-candidate status:\n"
  )

  print(
    table(
      validation$copykat.pred,
      validation$tumor_candidate,
      useNA = "ifany"
    )
  )

}


# ------------------------------------------------------------
# 23. Validate epithelial status
# ------------------------------------------------------------

if (
  "epithelial_status" %in%
    colnames(validation)
) {

  cat("\n")
  cat("====================================================\n")
  cat("EPITHELIAL VALIDATION\n")
  cat("====================================================\n")

  print(
    table(
      validation$copykat.pred,
      validation$epithelial_status,
      useNA = "ifany"
    )
  )

}


# ------------------------------------------------------------
# 24. High-confidence malignant candidate definition
# ------------------------------------------------------------

# IMPORTANT:
#
# We are NOT claiming that CopyKAT aneuploid = malignant.
#
# We create a conservative candidate flag requiring:
#
# 1. CopyKAT = aneuploid
# 2. tumor_candidate = TRUE
#
# Epithelial status is additionally reported when available.


validation$malignant_candidate <- FALSE


if (
  "tumor_candidate" %in%
    colnames(validation)
) {

  tumor_flag <- tolower(
    trimws(
      as.character(
        validation$tumor_candidate
      )
    )
  )

  tumor_positive <- tumor_flag %in% c(
    "true",
    "1",
    "yes",
    "tumor",
    "candidate"
  )

  validation$malignant_candidate <-
    validation$copykat.pred ==
      "aneuploid" &
    tumor_positive

}


# ------------------------------------------------------------
# 25. If tumor_candidate is unavailable
# ------------------------------------------------------------

if (
  !("tumor_candidate" %in%
      colnames(validation))
) {

  validation$malignant_candidate <-
    validation$copykat.pred ==
      "aneuploid"

}


# ------------------------------------------------------------
# 26. Malignant candidate summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("MALIGNANT CANDIDATE SUMMARY\n")
cat("====================================================\n")

cat(
  "CopyKAT aneuploid:",
  sum(
    validation$copykat.pred ==
      "aneuploid",
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "High-confidence malignant candidates:",
  sum(
    validation$malignant_candidate,
    na.rm = TRUE
  ),
  "\n"
)


# ------------------------------------------------------------
# 27. Save full validation table
# ------------------------------------------------------------

validation_file <- file.path(
  result_dir,
  "GSE228499_CopyKAT_validation_table.csv"
)

write.csv(
  validation,
  validation_file,
  row.names = FALSE
)

cat(
  "\nSaved validation table:\n",
  validation_file,
  "\n"
)


# ------------------------------------------------------------
# 28. Save malignant candidates
# ------------------------------------------------------------

malignant_file <- file.path(
  result_dir,
  "GSE228499_malignant_candidates.csv"
)

malignant_candidates <- validation %>%

  filter(
    malignant_candidate
  )

write.csv(
  malignant_candidates,
  malignant_file,
  row.names = FALSE
)

cat(
  "Saved malignant candidates:\n",
  malignant_file,
  "\n"
)


# ------------------------------------------------------------
# 29. CopyKAT prediction bar plot
# ------------------------------------------------------------

prediction_plot_data <- validation %>%

  count(
    copykat.pred
  ) %>%

  mutate(
    copykat.pred =
      factor(
        copykat.pred,
        levels = c(
          "diploid",
          "aneuploid",
          "not.defined"
        )
      )
  )


pdf(
  file.path(
    figure_dir,
    "CopyKAT_prediction_summary.pdf"
  ),
  width = 8,
  height = 6
)

print(

  ggplot(
    prediction_plot_data,
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
      title = "CopyKAT CNV Prediction Summary",
      x = "CopyKAT prediction",
      y = "Number of cells"
    ) +

    theme_classic()

)

dev.off()


# ------------------------------------------------------------
# 30. Prediction by sample plot
# ------------------------------------------------------------

pdf(
  file.path(
    figure_dir,
    "CopyKAT_prediction_by_sample.pdf"
  ),
  width = 10,
  height = 6
)

print(

  ggplot(
    sample_prediction,
    aes(
      x = sample_id,
      y = percentage,
      fill = copykat.pred
    )
  ) +

    geom_col(
      position = "stack"
    ) +

    labs(
      title = "CopyKAT Predictions Across Samples",
      x = "Sample",
      y = "Percentage of cells",
      fill = "CopyKAT prediction"
    ) +

    theme_classic()

)

dev.off()


# ------------------------------------------------------------
# 31. Aneuploid percentage plot
# ------------------------------------------------------------

pdf(
  file.path(
    figure_dir,
    "CopyKAT_aneuploid_percentage_by_sample.pdf"
  ),
  width = 10,
  height = 6
)

print(

  ggplot(
    aneuploid_by_sample,
    aes(
      x = sample_id,
      y = aneuploid_percent
    )
  ) +

    geom_col() +

    geom_text(
      aes(
        label =
          paste0(
            round(
              aneuploid_percent,
              1
            ),
            "%"
          )
      ),
      vjust = -0.4
    ) +

    labs(
      title = "Aneuploid Cell Fraction by Sample",
      x = "Sample",
      y = "Aneuploid cells (%)"
    ) +

    theme_classic()

)

dev.off()


# ------------------------------------------------------------
# 32. Save simple text report
# ------------------------------------------------------------

report_file <- file.path(
  result_dir,
  "CopyKAT_validation_summary.txt"
)

report_lines <- c(

  "COPYKAT VALIDATION SUMMARY",

  "==========================",

  paste(
    "Total CopyKAT predictions:",
    nrow(validation)
  ),

  paste(
    "Aneuploid:",
    sum(
      validation$copykat.pred ==
        "aneuploid",
      na.rm = TRUE
    )
  ),

  paste(
    "Diploid:",
    sum(
      validation$copykat.pred ==
        "diploid",
      na.rm = TRUE
    )
  ),

  paste(
    "Not defined:",
    sum(
      validation$copykat.pred ==
        "not.defined",
      na.rm = TRUE
    )
  ),

  paste(
    "High-confidence malignant candidates:",
    sum(
      validation$malignant_candidate,
      na.rm = TRUE
    )
  ),

  "",

  "Interpretation:",

  "CopyKAT aneuploid cells are malignant candidates,",

  "not automatically confirmed malignant cells.",

  "Final malignant-cell annotation should integrate",

  "CNV evidence with epithelial identity and marker evidence."

)


writeLines(
  report_lines,
  report_file
)


# ------------------------------------------------------------
# 33. Final summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("COPYKAT VALIDATION COMPLETED\n")
cat("====================================================\n")

cat(
  "Prediction table:",
  validation_file,
  "\n"
)

cat(
  "Malignant candidates:",
  malignant_file,
  "\n"
)

cat(
  "Figures:",
  figure_dir,
  "\n"
)

cat("\n")
cat(
  "IMPORTANT:\n",
  "Aneuploid cells are malignant candidates, not yet\n",
  "final confirmed malignant cells.\n",
  sep = ""
)

cat("\n")
cat("====================================================\n")
cat("READY FOR CNV + MARKER INTEGRATION\n")
cat("====================================================\n")