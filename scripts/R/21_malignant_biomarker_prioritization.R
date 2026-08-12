# ============================================================
# STEP 21: MALIGNANT BIOMARKER PRIORITIZATION
# ============================================================
#
# Project:
# AI-Assisted Discovery of Breast Cancer Biomarkers
# Using scRNA-seq and Machine Learning
#
# Dataset:
# GSE228499
#
# Purpose:
# Prioritize robust candidate biomarkers from the
# CNV-supported malignant epithelial population.
#
# Input:
# Step 20 differential-expression results
# Step 20 malignant characterization object
#
# Main criteria:
#   1. Adjusted P-value
#   2. Log2 fold-change
#   3. Malignant-cell prevalence
#   4. Diploid-cell prevalence
#   5. Malignant specificity
#   6. Sample consistency
#   7. Epithelial evidence
#
# IMPORTANT:
# This is a research prioritization framework.
# Candidates are NOT clinical biomarkers at this stage.
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

result_dir <- file.path(
  project_dir,
  "results",
  "malignant",
  "characterization"
)

figure_dir <- file.path(
  project_dir,
  "figures",
  "malignant",
  "characterization"
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

de_file <- file.path(
  result_dir,
  "20_malignant_vs_diploid_epithelial_DEG.csv"
)

seurat_file <- file.path(
  result_dir,
  "GSE228499_malignant_characterization.rds"
)


# ------------------------------------------------------------
# 4. Header
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("STEP 21: MALIGNANT BIOMARKER PRIORITIZATION\n")
cat("====================================================\n")


# ------------------------------------------------------------
# 5. Validate input
# ------------------------------------------------------------

if (!file.exists(de_file)) {

  stop(
    "\nDE result file not found:\n",
    de_file,
    "\n"
  )

}


if (!file.exists(seurat_file)) {

  stop(
    "\nMalignant Seurat object not found:\n",
    seurat_file,
    "\n"
  )

}


# ------------------------------------------------------------
# 6. Load DE results
# ------------------------------------------------------------

cat("\n")
cat("Loading Step 20 DE results...\n")

de <- read.csv(
  de_file,
  stringsAsFactors = FALSE
)


cat(
  "Genes loaded:",
  nrow(de),
  "\n"
)


# ------------------------------------------------------------
# 7. Validate DE columns
# ------------------------------------------------------------

required_columns <- c(

  "gene",
  "avg_log2FC",
  "pct.1",
  "pct.2",
  "p_val_adj"

)


missing_columns <- setdiff(

  required_columns,

  colnames(de)

)


if (length(missing_columns) > 0) {

  stop(

    "\nMissing required columns:\n",

    paste(
      missing_columns,
      collapse = "\n"
    ),

    "\n"

  )

}


# ------------------------------------------------------------
# 8. Load malignant object
# ------------------------------------------------------------

cat("\n")
cat("Loading malignant characterization object...\n")

malignant_obj <- readRDS(
  seurat_file
)


cat(
  "Malignant cells:",
  ncol(malignant_obj),
  "\n"
)


# ------------------------------------------------------------
# 9. Basic DE filtering
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("INITIAL DE FILTERING\n")
cat("====================================================\n")


de_filtered <-

  de %>%

  filter(

    !is.na(avg_log2FC),

    !is.na(pct.1),

    !is.na(pct.2),

    !is.na(p_val_adj),

    p_val_adj < 0.05

  )


cat(
  "Significant genes:",
  nrow(de_filtered),
  "\n"
)


# ------------------------------------------------------------
# 10. Calculate specificity
# ------------------------------------------------------------

de_filtered <-

  de_filtered %>%

  mutate(

    prevalence_difference =
      pct.1 - pct.2,

    absolute_log2FC =
      abs(avg_log2FC)

  )


# ------------------------------------------------------------
# 11. Define candidate direction
# ------------------------------------------------------------

de_filtered <-

  de_filtered %>%

  mutate(

    direction = case_when(

      avg_log2FC > 0 ~
        "Malignant_enriched",

      avg_log2FC < 0 ~
        "Diploid_enriched",

      TRUE ~
        "No_direction"

    )

  )


# ------------------------------------------------------------
# 12. Keep malignant-enriched genes
# ------------------------------------------------------------

malignant_de <-

  de_filtered %>%

  filter(

    avg_log2FC > 0

  )


cat(
  "Malignant-enriched genes:",
  nrow(malignant_de),
  "\n"
)


# ============================================================
# PART A
# CORE BIOMARKER FILTER
# ============================================================


# ------------------------------------------------------------
# 13. Apply conservative candidate thresholds
# ------------------------------------------------------------
#
# Criteria:
#
# adjusted P < 0.05
# log2FC >= 0.50
# malignant prevalence >= 25%
# prevalence difference >= 15 percentage points
#
# These thresholds are for research prioritization,
# not clinical validation.
# ------------------------------------------------------------

core_candidates <-

  malignant_de %>%

  filter(

    p_val_adj < 0.05,

    avg_log2FC >= 0.50,

    pct.1 >= 0.25,

    prevalence_difference >= 0.15

  )


cat("\n")
cat("Core candidates after filtering:\n")

cat(
  "Genes:",
  nrow(core_candidates),
  "\n"
)


# ------------------------------------------------------------
# 14. More stringent candidates
# ------------------------------------------------------------

high_confidence_candidates <-

  malignant_de %>%

  filter(

    p_val_adj < 0.01,

    avg_log2FC >= 1.0,

    pct.1 >= 0.40,

    prevalence_difference >= 0.20

  )


cat(
  "High-priority candidates:",
  nrow(high_confidence_candidates),
  "\n"
)


# ============================================================
# PART B
# SAMPLE CONSISTENCY
# ============================================================


# ------------------------------------------------------------
# 15. Extract malignant-cell expression by sample
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("SAMPLE CONSISTENCY ANALYSIS\n")
cat("====================================================\n")


sample_ids <- sort(

  unique(
    malignant_obj$sample_id
  )

)


cat(
  "Malignant samples:",
  paste(
    sample_ids,
    collapse = ", "
  ),
  "\n"
)


# ------------------------------------------------------------
# 16. Calculate gene prevalence by sample
# ------------------------------------------------------------

cat(
  "Calculating expression prevalence by sample...\n"
)


sample_consistency <- data.frame()


for (
  gene in core_candidates$gene
) {

  if (
    !gene %in%
    rownames(malignant_obj)
  ) {

    next

  }


  for (
    sample in sample_ids
  ) {

    sample_cells <- colnames(
      malignant_obj
    )[

      malignant_obj$sample_id ==
        sample

    ]


    if (
      length(sample_cells) == 0
    ) {

      next

    }


    expr <- FetchData(

      malignant_obj,

      vars =
        gene,

      cells =
        sample_cells

    )[, 1]


    prevalence <-

      mean(

        expr > 0,

        na.rm = TRUE

      )


    mean_expression <-

      mean(

        expr,

        na.rm = TRUE

      )


    sample_consistency <-

      rbind(

        sample_consistency,

        data.frame(

          gene =
            gene,

          sample_id =
            sample,

          expressing_fraction =
            prevalence,

          mean_expression =
            mean_expression

        )

      )

  }

}


# ------------------------------------------------------------
# 17. Save raw sample consistency
# ------------------------------------------------------------

write.csv(

  sample_consistency,

  file.path(

    result_dir,

    "21_sample_gene_expression_consistency.csv"

  ),

  row.names = FALSE

)


# ------------------------------------------------------------
# 18. Summarize consistency
# ------------------------------------------------------------

sample_consistency_summary <-

  sample_consistency %>%

  group_by(
    gene
  ) %>%

  summarise(

    samples_tested =
      n(),

    samples_prevalence_25pct =
      sum(
        expressing_fraction >= 0.25,
        na.rm = TRUE
      ),

    samples_prevalence_50pct =
      sum(
        expressing_fraction >= 0.50,
        na.rm = TRUE
      ),

    minimum_sample_prevalence =
      min(
        expressing_fraction,
        na.rm = TRUE
      ),

    mean_sample_prevalence =
      mean(
        expressing_fraction,
        na.rm = TRUE
      ),

    maximum_sample_prevalence =
      max(
        expressing_fraction,
        na.rm = TRUE
      ),

    .groups = "drop"

  )


# ------------------------------------------------------------
# 19. Merge consistency into candidate table
# ------------------------------------------------------------

core_candidates <-

  core_candidates %>%

  left_join(

    sample_consistency_summary,

    by = "gene"

  )


# ------------------------------------------------------------
# 20. Define consistency category
# ------------------------------------------------------------

core_candidates <-

  core_candidates %>%

  mutate(

    consistency_category = case_when(

      samples_prevalence_25pct == 3 ~
        "Consistent_3_of_3",

      samples_prevalence_25pct == 2 ~
        "Consistent_2_of_3",

      samples_prevalence_25pct == 1 ~
        "Sample_specific",

      TRUE ~
        "Not_detected"

    )

  )


# ============================================================
# PART C
# EPITHELIAL IDENTITY
# ============================================================


# ------------------------------------------------------------
# 21. Canonical epithelial markers
# ------------------------------------------------------------

epithelial_markers <- c(

  "EPCAM",

  "KRT8",

  "KRT18",

  "KRT19",

  "MUC1",

  "CEACAM6",

  "KRT14",

  "KRT17"

)


available_epithelial_markers <-

  epithelial_markers[

    epithelial_markers %in%

      rownames(
        malignant_obj
      )

  ]


cat("\n")
cat(
  "Available epithelial markers:",
  length(
    available_epithelial_markers
  ),
  "\n"
)


print(
  available_epithelial_markers
)


# ------------------------------------------------------------
# 22. Calculate epithelial marker expression
# ------------------------------------------------------------

epithelial_summary <- data.frame()


for (
  gene in available_epithelial_markers
) {

  expr <- FetchData(

    malignant_obj,

    vars =
      gene

  )[, 1]


  epithelial_summary <-

    rbind(

      epithelial_summary,

      data.frame(

        gene =
          gene,

        mean_expression =
          mean(
            expr,
            na.rm = TRUE
          ),

        expressing_fraction =
          mean(
            expr > 0,
            na.rm = TRUE
          )

      )

    )

}


write.csv(

  epithelial_summary,

  file.path(

    result_dir,

    "21_epithelial_marker_expression.csv"

  ),

  row.names = FALSE

)


# ============================================================
# PART D
# FLAG UNUSUAL / LOW-CONFIDENCE GENES
# ============================================================


# ------------------------------------------------------------
# 23. Genes requiring biological review
# ------------------------------------------------------------
#
# These are NOT automatically removed.
#
# They are flagged because they are unusual for a
# breast epithelial malignant-cell biomarker panel
# and require biological/contextual validation.
# ------------------------------------------------------------

review_genes <- c(

  "ALB",

  "PVALB",

  "KCNJ3",

  "CST5",

  "NEFH",

  "GRIA3",

  "CHAD",

  "AGTR1"

)


core_candidates <-

  core_candidates %>%

  mutate(

    biological_review_flag =

      ifelse(

        gene %in%
          review_genes,

        "Review_required",

        "Standard_candidate"

      )

  )


# ------------------------------------------------------------
# 24. Save flagged genes
# ------------------------------------------------------------

review_table <-

  core_candidates %>%

  filter(

    biological_review_flag ==
      "Review_required"

  )


write.csv(

  review_table,

  file.path(

    result_dir,

    "21_genes_requiring_biological_review.csv"

  ),

  row.names = FALSE

)


# ============================================================
# PART E
# BIOMARKER SCORING
# ============================================================


# ------------------------------------------------------------
# 25. Calculate normalized components
# ------------------------------------------------------------

core_candidates <-

  core_candidates %>%

  mutate(

    fc_score = pmin(

      avg_log2FC /
        2,

      1

    ),

    prevalence_score =
      pmin(

        pct.1,

        1

      ),

    specificity_score =
      pmin(

        pmax(
          prevalence_difference,
          0
        ),

        1

      ),

    significance_score =

      pmin(

        -log10(
          p_val_adj
        ) /

        20,

        1

      ),

    consistency_score =

      case_when(

        samples_prevalence_25pct == 3
        ~ 1.00,

        samples_prevalence_25pct == 2
        ~ 0.67,

        samples_prevalence_25pct == 1
        ~ 0.33,

        TRUE
        ~ 0

      )

  )


# ------------------------------------------------------------
# 26. Calculate final prioritization score
# ------------------------------------------------------------
#
# Weighting:
#
# Effect size       25%
# Malignant pct     20%
# Specificity       20%
# Significance      15%
# Sample consistency 20%
#
# This is a project-specific research score.
# It is NOT a validated clinical score.
# ------------------------------------------------------------

core_candidates <-

  core_candidates %>%

  mutate(

    biomarker_score =

      0.25 * fc_score +

      0.20 * prevalence_score +

      0.20 * specificity_score +

      0.15 * significance_score +

      0.20 * consistency_score

  )


# ------------------------------------------------------------
# 27. Apply review penalty
# ------------------------------------------------------------

core_candidates <-

  core_candidates %>%

  mutate(

    adjusted_biomarker_score =

      ifelse(

        biological_review_flag ==
          "Review_required",

        biomarker_score *
          0.50,

        biomarker_score

      )

  )


# ============================================================
# PART F
# FINAL RANKING
# ============================================================


# ------------------------------------------------------------
# 28. Rank candidates
# ------------------------------------------------------------

prioritized_candidates <-

  core_candidates %>%

  arrange(

    desc(
      adjusted_biomarker_score
    ),

    p_val_adj

  ) %>%

  mutate(

    priority_rank =
      row_number()

  )


# ------------------------------------------------------------
# 29. Save complete ranking
# ------------------------------------------------------------

write.csv(

  prioritized_candidates,

  file.path(

    result_dir,

    "21_prioritized_malignant_biomarkers.csv"

  ),

  row.names = FALSE

)


# ------------------------------------------------------------
# 30. Save top 50
# ------------------------------------------------------------

top_50 <-

  prioritized_candidates %>%

  head(50)


write.csv(

  top_50,

  file.path(

    result_dir,

    "21_top_50_malignant_biomarker_candidates.csv"

  ),

  row.names = FALSE

)


# ------------------------------------------------------------
# 31. Save top 20
# ------------------------------------------------------------

top_20 <-

  prioritized_candidates %>%

  head(20)


write.csv(

  top_20,

  file.path(

    result_dir,

    "21_top_20_malignant_biomarker_candidates.csv"

  ),

  row.names = FALSE

)


# ============================================================
# PART G
# CONSISTENCY FILTERED CANDIDATES
# ============================================================


# ------------------------------------------------------------
# 32. Robust 3-sample candidates
# ------------------------------------------------------------

robust_candidates <-

  prioritized_candidates %>%

  filter(

    samples_prevalence_25pct == 3,

    biological_review_flag ==
      "Standard_candidate"

  ) %>%

  arrange(

    desc(
      adjusted_biomarker_score
    )

  )


write.csv(

  robust_candidates,

  file.path(

    result_dir,

    "21_robust_3sample_biomarker_candidates.csv"

  ),

  row.names = FALSE

)


cat("\n")

cat(

  "Robust 3-sample candidates:",

  nrow(
    robust_candidates
  ),

  "\n"

)


# ============================================================
# PART H
# VISUALIZATION
# ============================================================


# ------------------------------------------------------------
# 33. Top candidate score plot
# ------------------------------------------------------------

plot_data <-

  top_20 %>%

  mutate(

    gene =
      reorder(
        gene,
        adjusted_biomarker_score
      )

  )


pdf(

  file.path(

    figure_dir,

    "21_top_biomarker_scores.pdf"

  ),

  width = 10,

  height = 8

)


print(

  ggplot(

    plot_data,

    aes(

      x =
        adjusted_biomarker_score,

      y =
        gene

    )

  ) +

    geom_col() +

    labs(

      title =
        "Top Prioritized Malignant Biomarker Candidates",

      x =
        "Adjusted biomarker score",

      y =
        "Gene"

    ) +

    theme_classic()

)


dev.off()


# ------------------------------------------------------------
# 34. Effect size vs specificity
# ------------------------------------------------------------

scatter_data <-

  prioritized_candidates %>%

  head(100)


pdf(

  file.path(

    figure_dir,

    "21_biomarker_effect_vs_specificity.pdf"

  ),

  width = 9,

  height = 7

)


print(

  ggplot(

    scatter_data,

    aes(

      x =
        prevalence_difference,

      y =
        avg_log2FC

    )

  ) +

    geom_point(

      alpha = 0.7

    ) +

    labs(

      title =
        "Biomarker Effect Size vs Malignant Specificity",

      x =
        "Prevalence difference (malignant - diploid)",

      y =
        "Average log2 fold change"

    ) +

    theme_classic()

)


dev.off()


# ============================================================
# PART I
# SUMMARY REPORT
# ============================================================


# ------------------------------------------------------------
# 35. Candidate category counts
# ------------------------------------------------------------

review_count <-

  sum(

    prioritized_candidates$biological_review_flag ==
      "Review_required"

  )


standard_count <-

  sum(

    prioritized_candidates$biological_review_flag ==
      "Standard_candidate"

  )


consistent_3sample <-

  sum(

    prioritized_candidates$samples_prevalence_25pct ==
      3

  )


# ------------------------------------------------------------
# 36. Write summary
# ------------------------------------------------------------

summary_file <- file.path(

  result_dir,

  "21_biomarker_prioritization_summary.txt"

)


summary_text <- c(

  "MALIGNANT BIOMARKER PRIORITIZATION SUMMARY",

  "===========================================",

  "",

  paste(
    "Significant DE genes:",
    nrow(de_filtered)
  ),

  "",

  paste(
    "Malignant-enriched genes:",
    nrow(malignant_de)
  ),

  "",

  paste(
    "Core biomarker candidates:",
    nrow(core_candidates)
  ),

  "",

  paste(
    "High-priority candidates:",
    nrow(high_confidence_candidates)
  ),

  "",

  paste(
    "3-sample consistent candidates:",
    consistent_3sample
  ),

  "",

  paste(
    "Standard candidates:",
    standard_count
  ),

  "",

  paste(
    "Genes requiring biological review:",
    review_count
  ),

  "",

  "Primary biomarker criteria:",

  "Adjusted P-value < 0.05",

  "Average log2FC >= 0.50",

  "Malignant prevalence >= 25%",

  "Prevalence difference >= 15 percentage points",

  "",

  "Sample consistency:",

  "Expression prevalence evaluated independently",

  "across BC11, BC12 and BC17.",

  "",

  "Scoring:",

  "Effect size = 25%",

  "Malignant prevalence = 20%",

  "Specificity = 20%",

  "Statistical significance = 15%",

  "Sample consistency = 20%",

  "",

  "Important:",

  "This prioritization score is a research",

  "framework and is not a validated clinical",

  "biomarker score.",

  "",

  "Genes flagged for biological review are",

  "not automatically excluded. They require",

  "additional biological and technical validation."

)


writeLines(

  summary_text,

  summary_file

)


# ============================================================
# FINAL OUTPUT
# ============================================================

cat("\n")
cat("====================================================\n")
cat("STEP 21 COMPLETED\n")
cat("====================================================\n")

cat(
  "Significant DE genes:",
  nrow(de_filtered),
  "\n"
)

cat(
  "Malignant-enriched genes:",
  nrow(malignant_de),
  "\n"
)

cat(
  "Core candidates:",
  nrow(core_candidates),
  "\n"
)

cat(
  "High-priority candidates:",
  nrow(high_confidence_candidates),
  "\n"
)

cat(
  "3-sample consistent candidates:",
  consistent_3sample,
  "\n"
)

cat("\n")
cat("Top 20 candidates:\n")

print(

  top_20[

    ,

    c(

      "priority_rank",

      "gene",

      "avg_log2FC",

      "pct.1",

      "pct.2",

      "prevalence_difference",

      "samples_prevalence_25pct",

      "adjusted_biomarker_score",

      "biological_review_flag"

    )

  ]

)


cat("\n")
cat("Results saved to:\n")

cat(
  result_dir,
  "\n"
)

cat("\n")
cat("====================================================\n")
cat("READY FOR PATHWAY + BIOLOGICAL VALIDATION\n")
cat("====================================================\n")