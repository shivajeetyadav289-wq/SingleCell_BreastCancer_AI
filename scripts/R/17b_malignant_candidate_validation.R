# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 17b_malignant_candidate_validation.R
#
# Purpose:
# Validate the CopyKAT + marker malignant candidates by
# examining:
#
#   1. CopyKAT prediction by sample
#   2. CopyKAT prediction by cell type
#   3. Epithelial identity
#   4. Tumor-candidate status
#   5. Luminal/basal marker scores
#   6. CNV profile coverage
#   7. Sample-specific aneuploidy
#
# IMPORTANT:
# This script does NOT:
#   - rerun CopyKAT
#   - change CopyKAT predictions
#   - overwrite malignant candidates
#   - create a final malignant Seurat object
#
# It is a validation and QC step.
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

integration_dir <- file.path(
  project_dir,
  "results",
  "cnv",
  "integration"
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

figure_dir <- file.path(
  project_dir,
  "figures",
  "cnv",
  "candidate_validation"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. Input files
# ------------------------------------------------------------

integration_file <- file.path(
  integration_dir,
  "GSE228499_CNV_marker_integrated_results.csv"
)

malignant_file <- file.path(
  integration_dir,
  "GSE228499_final_malignant_candidates.csv"
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
cat("MALIGNANT CANDIDATE VALIDATION\n")
cat("====================================================\n")

cat(
  "Dataset: GSE228499\n"
)

cat(
  "Purpose: CNV + marker quality control\n"
)


# ------------------------------------------------------------
# 5. Check files
# ------------------------------------------------------------

required_files <- c(
  integration_file,
  malignant_file,
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
# 6. Load integrated results
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("LOADING INTEGRATED RESULTS\n")
cat("====================================================\n")

validation <- read.csv(
  integration_file,
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


# ------------------------------------------------------------
# 7. Load malignant candidates
# ------------------------------------------------------------

malignant_candidates <- read.csv(
  malignant_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "Integrated malignant candidates:",
  nrow(malignant_candidates),
  "\n"
)


# ------------------------------------------------------------
# 8. Load CopyKAT prediction
# ------------------------------------------------------------

prediction <- readRDS(
  prediction_file
)

cat(
  "CopyKAT predictions:",
  nrow(prediction),
  "\n"
)


# ------------------------------------------------------------
# 9. Basic validation
# ------------------------------------------------------------

if (
  !"copykat.pred" %in%
    colnames(validation)
) {

  stop(
    "copykat.pred column not found."
  )

}

if (
  !"sample_id" %in%
    colnames(validation)
) {

  stop(
    "sample_id column not found."
  )

}


# ------------------------------------------------------------
# 10. Normalize strings
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
# 11. Overall CopyKAT summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("OVERALL COPYKAT SUMMARY\n")
cat("====================================================\n")

overall_prediction <- validation %>%

  count(
    copykat.pred
  ) %>%

  mutate(
    percentage =
      100 * n / sum(n)
  )

print(
  overall_prediction
)


write.csv(
  overall_prediction,
  file.path(
    integration_dir,
    "17b_overall_CopyKAT_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 12. CopyKAT by sample
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("COPYKAT BY SAMPLE\n")
cat("====================================================\n")

sample_prediction <- validation %>%

  group_by(
    sample_id
  ) %>%

  summarise(

    total = n(),

    diploid =
      sum(
        copykat.pred ==
          "diploid",
        na.rm = TRUE
      ),

    aneuploid =
      sum(
        copykat.pred ==
          "aneuploid",
        na.rm = TRUE
      ),

    not_defined =
      sum(
        copykat.pred ==
          "not.defined",
        na.rm = TRUE
      ),

    aneuploid_percent =
      100 *
      aneuploid /
      total,

    .groups = "drop"

  )


print(
  sample_prediction
)


write.csv(
  sample_prediction,
  file.path(
    integration_dir,
    "17b_CopyKAT_by_sample.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 13. Cell-type validation
# ------------------------------------------------------------

if (
  "cell_type" %in%
    colnames(validation)
) {

  cat("\n")
  cat("====================================================\n")
  cat("COPYKAT BY CELL TYPE\n")
  cat("====================================================\n")

  celltype_prediction <- validation %>%

    group_by(
      cell_type,
      copykat.pred
    ) %>%

    summarise(
      cells = n(),
      .groups = "drop"
    )

  print(
    celltype_prediction
  )

  write.csv(
    celltype_prediction,
    file.path(
      integration_dir,
      "17b_CopyKAT_by_cell_type.csv"
    ),
    row.names = FALSE
  )

}


# ------------------------------------------------------------
# 14. Tumor candidate validation
# ------------------------------------------------------------
if (
  "tumor_candidate_flag" %in%
    colnames(validation)
) {

  # Ensure the flag is logical

  validation$tumor_candidate_flag <-
    as.logical(
      validation$tumor_candidate_flag
    )

  cat("\n")
  cat("====================================================\n")
  cat("TUMOR-CANDIDATE × COPYKAT\n")
  cat("====================================================\n")

  cat(
    "Tumor-candidate status:\n"
  )

  print(
    table(
      validation$tumor_candidate_flag,
      useNA = "ifany"
    )
  )

  cat(
    "\nTumor-candidate × CopyKAT:\n"
  )

  tumor_cnv <- table(
    validation$tumor_candidate_flag,
    validation$copykat.pred,
    useNA = "ifany"
  )

  print(
    tumor_cnv
  )

}


# ------------------------------------------------------------
# 15. Epithelial validation
# ------------------------------------------------------------

if (
  "epithelial_evidence" %in%
    colnames(validation)
) {

  cat("\n")
  cat("====================================================\n")
  cat("EPITHELIAL EVIDENCE × COPYKAT\n")
  cat("====================================================\n")

  epithelial_cnv <- table(
    validation$epithelial_evidence,
    validation$copykat.pred
  )

  print(
    epithelial_cnv
  )

}
# ------------------------------------------------------------
# 15b. Tumor + epithelial population validation
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("TUMOR-CANDIDATE EPITHELIAL POPULATION\n")
cat("====================================================\n")

if (
  "tumor_candidate_flag" %in%
    colnames(validation) &&
  "epithelial_evidence" %in%
    colnames(validation)
) {

  combined_population <- validation %>%

    group_by(
      tumor_candidate_flag,
      epithelial_evidence,
      copykat.pred
    ) %>%

    summarise(
      cells = n(),
      .groups = "drop"
    )

  print(
    combined_population
  )

  write.csv(
    combined_population,
    file.path(
      integration_dir,
      "17b_tumor_epithelial_CopyKAT_summary.csv"
    ),
    row.names = FALSE
  )

}


# ------------------------------------------------------------
# 16. Aneuploid cell marker summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("ANEUPLOID CELL MARKER VALIDATION\n")
cat("====================================================\n")

aneuploid_cells <- validation %>%

  filter(
    copykat.pred ==
      "aneuploid"
  )


cat(
  "Aneuploid cells:",
  nrow(aneuploid_cells),
  "\n"
)


# ------------------------------------------------------------
# 17. Score summaries
# ------------------------------------------------------------

score_columns <- c(
  "epithelial_score",
  "luminal_score",
  "basal_score",
  "malignant_evidence_score"
)

available_score_columns <- intersect(
  score_columns,
  colnames(validation)
)


if (
  length(
    available_score_columns
  ) > 0
) {

  cat("\n")
  cat("Score summary for aneuploid cells:\n")

  for (
    score_name in available_score_columns
  ) {

    values <- suppressWarnings(
      as.numeric(
        aneuploid_cells[[score_name]]
      )
    )

    if (
      all(is.na(values))
    ) {

      next

    }

    cat(
      "\n",
      score_name,
      ":\n",
      sep = ""
    )

    print(
      summary(
        values
      )
    )

  }

}


# ------------------------------------------------------------
# 18. Compare scores by CopyKAT status
# ------------------------------------------------------------

score_comparison <- data.frame()

if (
  "epithelial_score" %in%
    colnames(validation)
) {

  score_comparison <- validation %>%

    group_by(
      copykat.pred
    ) %>%

    summarise(

      cells = n(),

      epithelial_score_mean =
        mean(
          epithelial_score,
          na.rm = TRUE
        ),

      epithelial_score_median =
        median(
          epithelial_score,
          na.rm = TRUE
        ),

      luminal_score_mean =
        if (
          "luminal_score" %in%
            colnames(validation)
        )
          mean(
            luminal_score,
            na.rm = TRUE
          )
        else
          NA_real_,

      basal_score_mean =
        if (
          "basal_score" %in%
            colnames(validation)
        )
          mean(
            basal_score,
            na.rm = TRUE
          )
        else
          NA_real_,

      .groups = "drop"

    )

}


if (
  nrow(score_comparison) > 0
) {

  cat("\n")
  cat("====================================================\n")
  cat("MARKER SCORE COMPARISON\n")
  cat("====================================================\n")

  print(
    score_comparison
  )

  write.csv(
    score_comparison,
    file.path(
      integration_dir,
      "17b_marker_score_comparison.csv"
    ),
    row.names = FALSE
  )

}


# ------------------------------------------------------------
# 19. Sample-specific aneuploidity check
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("SAMPLE-SPECIFIC ANEUPLOIDITY CHECK\n")
cat("====================================================\n")

sample_check <- validation %>%

  group_by(
    sample_id
  ) %>%

  summarise(

    cells = n(),

    aneuploid_cells =
      sum(
        copykat.pred ==
          "aneuploid",
        na.rm = TRUE
      ),

    aneuploid_fraction =
      mean(
        copykat.pred ==
          "aneuploid",
        na.rm = TRUE
      ),

    tumor_candidates =
      if (
        "tumor_candidate_flag" %in%
          colnames(validation)
      )
        sum(
          tumor_candidate_flag,
          na.rm = TRUE
        )
      else
        NA_integer_,

    epithelial_cells =
      if (
        "epithelial_evidence" %in%
          colnames(validation)
      )
        sum(
          epithelial_evidence,
          na.rm = TRUE
        )
      else
        NA_integer_,

    .groups = "drop"

  )


sample_check$aneuploid_fraction <-
  round(
    100 *
      sample_check$aneuploid_fraction,
    2
  )


print(
  sample_check
)


write.csv(
  sample_check,
  file.path(
    integration_dir,
    "17b_sample_specific_aneuploidity.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 20. Plot: aneuploid fraction by sample
# ------------------------------------------------------------

pdf(
  file.path(
    figure_dir,
    "17b_aneuploid_fraction_by_sample.pdf"
  ),
  width = 10,
  height = 6
)

print(

  ggplot(
    sample_check,
    aes(
      x = sample_id,
      y = aneuploid_fraction
    )
  ) +

    geom_col() +

    geom_text(
      aes(
        label =
          paste0(
            aneuploid_fraction,
            "%"
          )
      ),
      vjust = -0.4
    ) +

    labs(
      title =
        "CopyKAT Aneuploid Fraction by Sample",
      x = "Sample",
      y = "Aneuploid cells (%)"
    ) +

    theme_classic()

)

dev.off()


# ------------------------------------------------------------
# 21. Plot: CopyKAT prediction by sample
# ------------------------------------------------------------

plot_sample <- validation %>%

  count(
    sample_id,
    copykat.pred
  ) %>%

  group_by(
    sample_id
  ) %>%

  mutate(
    percentage =
      100 *
      n /
      sum(n)
  ) %>%

  ungroup()


pdf(
  file.path(
    figure_dir,
    "17b_CopyKAT_prediction_by_sample.pdf"
  ),
  width = 11,
  height = 6
)

print(

  ggplot(
    plot_sample,
    aes(
      x = sample_id,
      y = percentage,
      fill = copykat.pred
    )
  ) +

    geom_col() +

    labs(
      title =
        "CopyKAT Prediction Distribution by Sample",
      x = "Sample",
      y = "Cells (%)",
      fill = "CopyKAT"
    ) +

    theme_classic()

)

dev.off()


# ------------------------------------------------------------
# 22. Plot: epithelial score vs CNV status
# ------------------------------------------------------------

if (
  "epithelial_score" %in%
    colnames(validation)
) {

  pdf(
    file.path(
      figure_dir,
      "17b_epithelial_score_by_CopyKAT.pdf"
    ),
    width = 9,
    height = 6
  )

  print(

    ggplot(
      validation,
      aes(
        x = copykat.pred,
        y = epithelial_score
      )
    ) +

      geom_boxplot(
        outlier.size = 0.3
      ) +

      labs(
        title =
          "Epithelial Score by CopyKAT Prediction",
        x = "CopyKAT prediction",
        y = "Epithelial score"
      ) +

      theme_classic()

  )

  dev.off()

}


# ------------------------------------------------------------
# 23. Plot: luminal/basal scores
# ------------------------------------------------------------

if (
  "luminal_score" %in%
    colnames(validation) &&
  "basal_score" %in%
    colnames(validation)
) {

  score_long <- rbind(

    data.frame(
      cell_id =
        validation$cell_id,
      copykat =
        validation$copykat.pred,
      marker_type =
        "Luminal",
      score =
        validation$luminal_score
    ),

    data.frame(
      cell_id =
        validation$cell_id,
      copykat =
        validation$copykat.pred,
      marker_type =
        "Basal",
      score =
        validation$basal_score
    )

  )


  pdf(
    file.path(
      figure_dir,
      "17b_luminal_basal_scores_by_CopyKAT.pdf"
    ),
    width = 10,
    height = 6
  )

  print(

    ggplot(
      score_long,
      aes(
        x = copykat,
        y = score
      )
    ) +

      geom_boxplot(
        outlier.size = 0.2
      ) +

      facet_wrap(
        ~ marker_type,
        scales = "free_y"
      ) +

      labs(
        title =
          "Luminal and Basal Marker Scores by CopyKAT",
        x = "CopyKAT prediction",
        y = "Score"
      ) +

      theme_classic()

  )

  dev.off()

}


# ------------------------------------------------------------
# 24. Save aneuploid validation table
# ------------------------------------------------------------

aneuploid_file <- file.path(
  integration_dir,
  "17b_aneuploid_validation_table.csv"
)

write.csv(
  aneuploid_cells,
  aneuploid_file,
  row.names = FALSE
)


# ------------------------------------------------------------
# 25. Save malignant candidate validation table
# ------------------------------------------------------------

candidate_validation <- validation %>%

  filter(
    final_malignant_candidate
  )


candidate_file <- file.path(
  integration_dir,
  "17b_malignant_candidate_validation.csv"
)

write.csv(
  candidate_validation,
  candidate_file,
  row.names = FALSE
)


# ------------------------------------------------------------
# 26. Identify samples with strong aneuploid signal
# ------------------------------------------------------------

strong_samples <- sample_check %>%

  filter(
    aneuploid_fraction >= 50
  )


cat("\n")
cat("====================================================\n")
cat("SAMPLES WITH >50% ANEUPLOID CELLS\n")
cat("====================================================\n")

if (
  nrow(strong_samples) > 0
) {

  print(
    strong_samples
  )

} else {

  cat(
    "No sample exceeds 50% aneuploid cells.\n"
  )

}


# ------------------------------------------------------------
# 27. Final report
# ------------------------------------------------------------

report_file <- file.path(
  integration_dir,
  "17b_malignant_candidate_validation_summary.txt"
)


report_lines <- c(

  "MALIGNANT CANDIDATE VALIDATION SUMMARY",

  "======================================",

  "",

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

  "",

  paste(
    "Integrated malignant candidates:",
    sum(
      validation$final_malignant_candidate,
      na.rm = TRUE
    )
  ),

  "",

  "Important observations:",

  "CopyKAT aneuploidy is interpreted as CNV evidence.",

  "Aneuploid cells are not automatically considered malignant.",

  "Epithelial and tumor-candidate evidence must support",

  "the malignant-cell interpretation.",

  "",

  "Sample-specific aneuploidity must be investigated",

  "before final malignant-cell labeling."

)


writeLines(
  report_lines,
  report_file
)


# ------------------------------------------------------------
# 28. Final summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("17B VALIDATION COMPLETED\n")
cat("====================================================\n")

cat(
  "Integrated results:",
  integration_file,
  "\n"
)

cat(
  "Aneuploid validation:",
  aneuploid_file,
  "\n"
)

cat(
  "Candidate validation:",
  candidate_file,
  "\n"
)

cat(
  "Summary:",
  report_file,
  "\n"
)

cat("\nFigures:\n")

cat(
  figure_dir,
  "\n"
)

cat("\n")
cat("====================================================\n")
cat("DO NOT FINALIZE MALIGNANT LABELS YET\n")
cat("====================================================\n")

cat(
  "Review BC11, BC12 and BC17 CNV patterns and\n",
  "marker evidence before transferring labels to Seurat.\n",
  sep = ""
)