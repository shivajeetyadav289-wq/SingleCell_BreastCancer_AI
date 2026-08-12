# ============================================================
# STEP 25
# FINAL CANDIDATE BIOLOGICAL REVIEW
# ============================================================
#
# Purpose:
#   Perform final biological review using the COMPLETE
#   Step-24 integrated ranking (986 candidates).
#
# Important:
#   This is a research prioritization framework.
#   It is NOT a clinical biomarker validation.
#
# ============================================================


cat("\n")
cat("====================================================\n")
cat("STEP 25: FINAL CANDIDATE BIOLOGICAL REVIEW\n")
cat("====================================================\n")


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})


# ------------------------------------------------------------
# 2. Paths
# ------------------------------------------------------------

input_file <- file.path(
  "results",
  "malignant",
  "integrated_prioritization",
  "24_integrated_malignant_biomarker_ranking.csv"
)

output_dir <- file.path(
  "results",
  "malignant",
  "final_candidates"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. Read complete Step-24 ranking
# ------------------------------------------------------------

cat("\nReading complete Step-24 ranking...\n")

x <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "Step-24 candidates:",
  nrow(x),
  "\n"
)


# ------------------------------------------------------------
# 4. Basic validation
# ------------------------------------------------------------

required_columns <- c(
  "gene",
  "final_integrated_score",
  "adjusted_biomarker_score",
  "avg_log2FC",
  "pct.1",
  "pct.2",
  "prevalence_difference",
  "mean_sample_prevalence",
  "minimum_sample_prevalence",
  "maximum_sample_prevalence",
  "consistency_category",
  "sample_consistency_score",
  "specificity_score",
  "prevalence_score",
  "significance_score",
  "high_priority",
  "robust_3sample",
  "final_biological_category",
  "biological_class",
  "GO_BP_support",
  "GO_CC_support",
  "GO_MF_support",
  "KEGG_support",
  "pathway_support",
  "epithelial_relevance",
  "proliferation_relevance",
  "mitochondrial_relevance",
  "stress_relevance",
  "metabolic_flag",
  "biological_review_flag"
)


missing_columns <- setdiff(
  required_columns,
  colnames(x)
)


if (length(missing_columns) > 0) {

  stop(
    paste(
      "Missing required columns:",
      paste(missing_columns, collapse = ", ")
    )
  )

}


# ------------------------------------------------------------
# 5. Standardize logical columns
# ------------------------------------------------------------

logical_columns <- c(
  "high_priority",
  "robust_3sample",
  "GO_BP_support",
  "GO_CC_support",
  "GO_MF_support",
  "KEGG_support",
  "pathway_support",
  "epithelial_relevance",
  "proliferation_relevance",
  "mitochondrial_relevance",
  "stress_relevance",
  "metabolic_flag"
)


for (col in logical_columns) {

  if (col %in% colnames(x)) {

    x[[col]] <- as.logical(
      x[[col]]
    )

  }

}


# ------------------------------------------------------------
# 6. Remove duplicated genes
# ------------------------------------------------------------

duplicate_genes <- duplicated(
  x$gene
)

cat(
  "Duplicated genes:",
  sum(duplicate_genes),
  "\n"
)


x <- x[
  !duplicate_genes,
]


# ------------------------------------------------------------
# 7. Strong candidate definition
# ------------------------------------------------------------
#
# Strong candidate:
#
#   adjusted biomarker score >= 0.50
#   absolute log2FC >= 0.50
#   malignant prevalence >= 25%
#   prevalence difference >= 15 percentage points
#   AND
#   consistent in at least 2 of 3 samples
#
# ------------------------------------------------------------

x$strong_candidate <- (
  x$adjusted_biomarker_score >= 0.50 &
  x$absolute_log2FC >= 0.50 &
  x$mean_sample_prevalence >= 0.25 &
  x$prevalence_difference >= 0.15 &
  x$consistency_category %in%
    c(
      "Consistent_2_of_3",
      "Consistent_3_of_3"
    )
)


cat(
  "Strong candidates:",
  sum(x$strong_candidate, na.rm = TRUE),
  "\n"
)


# ------------------------------------------------------------
# 8. Strong 3-sample candidates
# ------------------------------------------------------------

x$strong_3sample <- (
  x$strong_candidate &
  x$consistency_category ==
    "Consistent_3_of_3"
)


cat(
  "Strong 3-sample candidates:",
  sum(x$strong_3sample, na.rm = TRUE),
  "\n"
)


# ------------------------------------------------------------
# 9. Biological categories
# ------------------------------------------------------------

x$review_category <- "Other_Research_Candidate"


# Tumor / epithelial has highest biological priority
x$review_category[
  x$epithelial_relevance == TRUE
] <- "Tumor_Epithelial"


# Proliferation
x$review_category[
  x$proliferation_relevance == TRUE &
  x$review_category ==
    "Other_Research_Candidate"
] <- "Proliferation"


# Mitochondrial / OXPHOS
x$review_category[
  x$mitochondrial_relevance == TRUE &
  x$review_category ==
    "Other_Research_Candidate"
] <- "Mitochondrial_OXPHOS"


# Stress
x$review_category[
  x$stress_relevance == TRUE &
  x$review_category ==
    "Other_Research_Candidate"
] <- "Stress_Response"


# ------------------------------------------------------------
# 10. Biological review status
# ------------------------------------------------------------

x$final_review_status <- "Standard"


x$final_review_status[
  x$biological_review_flag ==
    "Review_required"
] <- "Biological_review_required"


# ------------------------------------------------------------
# 11. Biological support score
# ------------------------------------------------------------

x$pathway_support_score <- rowSums(
  cbind(
    x$GO_BP_support,
    x$GO_CC_support,
    x$GO_MF_support,
    x$KEGG_support,
    x$pathway_support
  ),
  na.rm = TRUE
)


# ------------------------------------------------------------
# 12. Biological relevance score
# ------------------------------------------------------------

x$biological_relevance_score <- 0


# Tumor/epithelial evidence
x$biological_relevance_score <-
  x$biological_relevance_score +
  ifelse(
    x$epithelial_relevance,
    0.25,
    0
  )


# Proliferation
x$biological_relevance_score <-
  x$biological_relevance_score +
  ifelse(
    x$proliferation_relevance,
    0.15,
    0
  )


# Pathway support
x$biological_relevance_score <-
  x$biological_relevance_score +
  pmin(
    x$pathway_support_score / 5,
    1
  ) * 0.20


# Sample consistency
x$biological_relevance_score <-
  x$biological_relevance_score +
  x$sample_consistency_score * 0.20


# Biomarker score
x$biological_relevance_score <-
  x$biological_relevance_score +
  x$adjusted_biomarker_score * 0.20


# ------------------------------------------------------------
# 13. Penalize sample-specific candidates
# ------------------------------------------------------------

x$consistency_penalty <- 0

x$consistency_penalty[
  x$consistency_category ==
    "Sample_specific"
] <- 0.15


# ------------------------------------------------------------
# 14. Penalize biological review candidates
# ------------------------------------------------------------

x$review_penalty <- 0

x$review_penalty[
  x$biological_review_flag ==
    "Review_required"
] <- 0.10


# ------------------------------------------------------------
# 15. Final integrated biological score
# ------------------------------------------------------------

x$final_biological_score <- pmax(
  0,
  x$biological_relevance_score -
    x$consistency_penalty -
    x$review_penalty
)


# ------------------------------------------------------------
# 16. Tumor/epithelial priority score
# ------------------------------------------------------------
#
# Tumor/epithelial candidates are ranked separately.
#
# This prevents OXPHOS genes from being interpreted
# automatically as tumor-specific biomarkers.
#
# ------------------------------------------------------------

x$tumor_priority_score <- (
  x$adjusted_biomarker_score * 0.35 +
  x$sample_consistency_score * 0.25 +
  x$specificity_score * 0.20 +
  x$prevalence_score * 0.20
)


# ------------------------------------------------------------
# 17. OXPHOS priority score
# ------------------------------------------------------------

x$OXPHOS_priority_score <- (
  x$adjusted_biomarker_score * 0.30 +
  x$sample_consistency_score * 0.25 +
  x$prevalence_score * 0.20 +
  x$significance_score * 0.15 +
  pmin(
    x$pathway_support_score / 5,
    1
  ) * 0.10
)


# ------------------------------------------------------------
# 18. Final candidate ranking
# ------------------------------------------------------------

x <- x %>%
  arrange(
    desc(strong_3sample),
    desc(strong_candidate),
    desc(final_biological_score),
    desc(adjusted_biomarker_score)
  )


x$final_panel_rank <- seq_len(
  nrow(x)
)


# ------------------------------------------------------------
# 19. Final research candidate selection
# ------------------------------------------------------------
#
# We do NOT use only the previous top 50.
#
# Selection is performed from the entire 986 candidate
# dataset.
#
# ------------------------------------------------------------

strong_pool <- x %>%
  filter(
    strong_candidate == TRUE
  )


cat(
  "Strong candidate pool:",
  nrow(strong_pool),
  "\n"
)


# ------------------------------------------------------------
# 20. Final research panel
# ------------------------------------------------------------
#
# Prioritize:
#
#   1. Strong 3-sample tumor/epithelial
#   2. Strong 3-sample other candidates
#   3. Strong 3-sample OXPHOS
#   4. Strong 2-of-3 candidates
#
# Maximum panel size = 50
#
# ------------------------------------------------------------

tumor_pool <- strong_pool %>%
  filter(
    review_category ==
      "Tumor_Epithelial",
    strong_3sample == TRUE
  ) %>%
  arrange(
    desc(tumor_priority_score),
    desc(final_biological_score)
  )


OXPHOS_pool <- strong_pool %>%
  filter(
    review_category ==
      "Mitochondrial_OXPHOS",
    strong_3sample == TRUE
  ) %>%
  arrange(
    desc(OXPHOS_priority_score),
    desc(final_biological_score)
  )


other_pool <- strong_pool %>%
  filter(
    review_category ==
      "Other_Research_Candidate",
    strong_3sample == TRUE
  ) %>%
  arrange(
    desc(final_biological_score)
  )


proliferation_pool <- strong_pool %>%
  filter(
    review_category ==
      "Proliferation",
    strong_3sample == TRUE
  ) %>%
  arrange(
    desc(final_biological_score)
  )


# ------------------------------------------------------------
# 21. Allocate final panel
# ------------------------------------------------------------

panel_parts <- list()


# Tumor/epithelial:
# keep all strong 3-sample candidates
if (nrow(tumor_pool) > 0) {

  panel_parts[[length(panel_parts) + 1]] <-
    tumor_pool

}


# Proliferation:
if (nrow(proliferation_pool) > 0) {

  panel_parts[[length(panel_parts) + 1]] <-
    proliferation_pool

}


# OXPHOS:
if (nrow(OXPHOS_pool) > 0) {

  panel_parts[[length(panel_parts) + 1]] <-
    OXPHOS_pool

}


# Other candidates:
if (nrow(other_pool) > 0) {

  panel_parts[[length(panel_parts) + 1]] <-
    other_pool

}


if (length(panel_parts) > 0) {

  panel <- bind_rows(
    panel_parts
  )

} else {

  panel <- x[
    FALSE,
  ]

}


# ------------------------------------------------------------
# 22. Remove duplicate genes
# ------------------------------------------------------------

panel <- panel %>%
  distinct(
    gene,
    .keep_all = TRUE
  )


# ------------------------------------------------------------
# 23. Limit to top 50
# ------------------------------------------------------------

panel <- panel %>%
  arrange(
    desc(
      strong_3sample
    ),
    desc(
      final_biological_score
    )
  ) %>%
  slice_head(
    n = 50
  )


panel$final_panel_rank <- seq_len(
  nrow(panel)
)


# ------------------------------------------------------------
# 24. Final review class
# ------------------------------------------------------------

panel$final_review_class <- panel$review_category


panel$final_tier <- "Research_candidate"


panel$final_tier[
  panel$review_category ==
    "Tumor_Epithelial"
] <- "Tier_1_Tumor_Epithelial"


panel$final_tier[
  panel$review_category ==
    "Proliferation"
] <- "Tier_1_Proliferation"


panel$final_tier[
  panel$review_category ==
    "Mitochondrial_OXPHOS"
] <- "Tier_2_OXPHOS"


panel$final_tier[
  panel$review_category ==
    "Other_Research_Candidate"
] <- "Tier_3_Other"


# ------------------------------------------------------------
# 25. Save complete ranking
# ------------------------------------------------------------

write.csv(
  x,
  file.path(
    output_dir,
    "25_complete_final_candidate_ranking.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 26. Save final panel
# ------------------------------------------------------------

write.csv(
  panel,
  file.path(
    output_dir,
    "25_final_research_candidate_panel.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 27. Save top 50
# ------------------------------------------------------------

write.csv(
  panel,
  file.path(
    output_dir,
    "25_top50_final_research_candidates.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 28. Tumor / epithelial candidates
# ------------------------------------------------------------

tumor_final <- panel %>%
  filter(
    review_category ==
      "Tumor_Epithelial"
  )


write.csv(
  tumor_final,
  file.path(
    output_dir,
    "25_final_tumor_epithelial_candidates.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 29. OXPHOS candidates
# ------------------------------------------------------------

OXPHOS_final <- panel %>%
  filter(
    review_category ==
      "Mitochondrial_OXPHOS"
  )


write.csv(
  OXPHOS_final,
  file.path(
    output_dir,
    "25_final_OXPHOS_candidates.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 30. Proliferation candidates
# ------------------------------------------------------------

proliferation_final <- panel %>%
  filter(
    review_category ==
      "Proliferation"
  )


write.csv(
  proliferation_final,
  file.path(
    output_dir,
    "25_final_proliferation_candidates.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 31. Category summary
# ------------------------------------------------------------

category_summary <- panel %>%
  count(
    review_category,
    name = "genes"
  )


write.csv(
  category_summary,
  file.path(
    output_dir,
    "25_final_panel_category_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 32. Strong candidate summary
# ------------------------------------------------------------

strong_summary <- data.frame(
  total_step24_candidates = nrow(x),

  strong_candidates =
    sum(
      x$strong_candidate,
      na.rm = TRUE
    ),

  strong_3sample_candidates =
    sum(
      x$strong_3sample,
      na.rm = TRUE
    ),

  consistent_3_of_3 =
    sum(
      x$consistency_category ==
        "Consistent_3_of_3",
      na.rm = TRUE
    ),

  consistent_2_of_3 =
    sum(
      x$consistency_category ==
        "Consistent_2_of_3",
      na.rm = TRUE
    ),

  sample_specific =
    sum(
      x$consistency_category ==
        "Sample_specific",
      na.rm = TRUE
    ),

  final_panel_size =
    nrow(panel),

  tumor_epithelial_panel =
    sum(
      panel$review_category ==
        "Tumor_Epithelial"
    ),

  OXPHOS_panel =
    sum(
      panel$review_category ==
        "Mitochondrial_OXPHOS"
    ),

  proliferation_panel =
    sum(
      panel$review_category ==
        "Proliferation"
    ),

  other_panel =
    sum(
      panel$review_category ==
        "Other_Research_Candidate"
    )
)


write.csv(
  strong_summary,
  file.path(
    output_dir,
    "25_final_candidate_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 33. Final report
# ------------------------------------------------------------

report_file <- file.path(
  output_dir,
  "25_final_candidate_biological_review_summary.txt"
)


sink(
  report_file
)


cat("FINAL CANDIDATE BIOLOGICAL REVIEW SUMMARY\n")
cat("=========================================\n\n")


cat(
  "Step 24 integrated candidates:",
  nrow(x),
  "\n"
)


cat(
  "Strong candidates:",
  sum(
    x$strong_candidate,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Strong 3-sample candidates:",
  sum(
    x$strong_3sample,
    na.rm = TRUE
  ),
  "\n\n"
)


cat("Consistency distribution:\n")
print(
  table(
    x$consistency_category
  )
)


cat("\nBiological categories in complete Step-24 dataset:\n")
print(
  table(
    x$final_biological_category
  )
)


cat("\nFinal panel size:")
cat(
  nrow(panel),
  "\n\n"
)


cat("FINAL PANEL CATEGORY DISTRIBUTION:\n")
print(
  table(
    panel$review_category
  )
)


cat("\nFINAL RESEARCH PANEL:\n\n")


print(
  panel[
    ,
    c(
      "final_panel_rank",
      "gene",
      "final_biological_score",
      "final_tier",
      "avg_log2FC",
      "pct.1",
      "mean_sample_prevalence",
      "consistency_category",
      "review_category"
    )
  ]
)


cat("\n\nIMPORTANT INTERPRETATION:\n")

cat(
  "The final panel is a research-oriented shortlist ",
  "derived from the complete Step-24 integrated ",
  "candidate ranking.\n\n"
)


cat(
  "Tumor/epithelial candidates are prioritized ",
  "separately from mitochondrial/OXPHOS candidates ",
  "because metabolic programs are not automatically ",
  "tumor-specific.\n\n"
)


cat(
  "Mitochondrial/OXPHOS candidates represent ",
  "biologically informative metabolic programs and ",
  "require additional validation before being ",
  "considered cancer-specific biomarkers.\n\n"
)


cat(
  "Sample consistency was evaluated across BC11, ",
  "BC12 and BC17.\n\n"
)


cat(
  "All classifications are research annotations and ",
  "are not clinical diagnostic classifications.\n"
)


sink()


# ------------------------------------------------------------
# 34. Console output
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("STEP 25 COMPLETED\n")
cat("====================================================\n\n")


cat(
  "Complete Step-24 candidates:",
  nrow(x),
  "\n"
)


cat(
  "Strong candidates:",
  sum(
    x$strong_candidate,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Strong 3-sample candidates:",
  sum(
    x$strong_3sample,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Final panel:",
  nrow(panel),
  "\n\n"
)


cat("Final panel categories:\n")

print(
  table(
    panel$review_category
  )
)


cat("\nFinal panel genes:\n")

print(
  panel$gene
)


cat("\nOutput directory:\n")
cat(
  output_dir,
  "\n"
)


cat("\n====================================================\n")