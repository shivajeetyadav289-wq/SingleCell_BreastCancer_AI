# ============================================================
# STEP 24
# INTEGRATED MALIGNANT BIOMARKER PRIORITIZATION
# ============================================================
#
# Project:
# SingleCell_BreastCancer_AI
#
# Purpose:
# Integrate Steps 21, 22 and 23 into a final research-level
# biomarker candidate ranking.
#
# This script DOES NOT rerun:
#   - CopyKAT
#   - differential expression
#   - Seurat clustering
#   - pathway enrichment
#
# Evidence integrated:
#   1. Differential expression
#   2. Malignant prevalence
#   3. Specificity
#   4. Statistical significance
#   5. Sample consistency
#   6. Existing biomarker score
#   7. Biological category
#   8. Pathway support
#   9. Biological review flags
#
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

suppressPackageStartupMessages({

  library(dplyr)
  library(ggplot2)

})


# ------------------------------------------------------------
# 2. Directories
# ------------------------------------------------------------

project_dir <- path.expand(
  "~/Projects/SingleCell_BreastCancer_AI"
)

characterization_dir <- file.path(
  project_dir,
  "results",
  "malignant",
  "characterization"
)

pathway_dir <- file.path(
  project_dir,
  "results",
  "malignant",
  "pathway_enrichment"
)

biological_dir <- file.path(
  project_dir,
  "results",
  "malignant",
  "biological_validation"
)

output_dir <- file.path(
  project_dir,
  "results",
  "malignant",
  "integrated_prioritization"
)

figure_dir <- file.path(
  project_dir,
  "figures",
  "malignant",
  "integrated_prioritization"
)


dir.create(
  output_dir,
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
cat("STEP 24: INTEGRATED BIOMARKER PRIORITIZATION\n")
cat("====================================================\n")


# ============================================================
# PART A
# LOAD STEP 21
# ============================================================

step21_file <- file.path(
  characterization_dir,
  "21_prioritized_malignant_biomarkers.csv"
)


if (!file.exists(step21_file)) {

  stop(
    "Step 21 file not found:\n",
    step21_file
  )

}


step21 <- read.csv(
  step21_file,
  stringsAsFactors = FALSE
)


cat("\n")
cat(
  "Step 21 candidates:",
  nrow(step21),
  "\n"
)


# ============================================================
# PART B
# LOAD STEP 23 BIOLOGICAL ANNOTATION
# ============================================================

step23_file <- file.path(
  biological_dir,
  "23_all_candidate_biological_annotation.csv"
)


if (!file.exists(step23_file)) {

  stop(
    "Step 23 biological annotation not found:\n",
    step23_file
  )

}


step23 <- read.csv(
  step23_file,
  stringsAsFactors = FALSE
)


cat(
  "Step 23 annotated candidates:",
  nrow(step23),
  "\n"
)


# ============================================================
# PART C
# MERGE STEP 21 + STEP 23
# ============================================================

integrated <- step21 %>%

  left_join(
    step23 %>%
      select(
        gene,
        high_priority,
        robust_3sample,
        biological_class,
        final_biological_category
      ),
    by = "gene"
  )


# ------------------------------------------------------------
# Missing biological annotation
# ------------------------------------------------------------

integrated$high_priority[
  is.na(integrated$high_priority)
] <- FALSE


integrated$robust_3sample[
  is.na(integrated$robust_3sample)
] <- FALSE


integrated$biological_class[
  is.na(integrated$biological_class)
] <- "Not_annotated"


integrated$final_biological_category[
  is.na(integrated$final_biological_category)
] <- "Other_Candidate"


# ============================================================
# PART D
# PATHWAY SUPPORT
# ============================================================

go_bp_file <- file.path(
  pathway_dir,
  "22_robust_3sample_GO_BP.csv"
)


go_cc_file <- file.path(
  pathway_dir,
  "22_robust_3sample_GO_CC.csv"
)


go_mf_file <- file.path(
  pathway_dir,
  "22_robust_3sample_GO_MF.csv"
)


kegg_file <- file.path(
  pathway_dir,
  "22_robust_3sample_KEGG.csv"
)


# ------------------------------------------------------------
# Function to extract genes from clusterProfiler geneID
# ------------------------------------------------------------

extract_genes_from_enrichment <- function(
  file
) {

  if (!file.exists(file)) {

    return(
      character(0)
    )

  }


  x <- read.csv(
    file,
    stringsAsFactors = FALSE
  )


  if (
    !"geneID" %in%
    colnames(x)
  ) {

    return(
      character(0)
    )

  }


  gene_strings <- x$geneID

  gene_strings <- gene_strings[
    !is.na(gene_strings)
  ]


  if (
    length(gene_strings) == 0
  ) {

    return(
      character(0)
    )

  }


  genes <- unlist(

    strsplit(
      gene_strings,
      "/",
      fixed = TRUE
    )

  )


  unique(
    genes
  )

}


# ------------------------------------------------------------
# Extract pathway-supported genes
# ------------------------------------------------------------

go_bp_genes <- extract_genes_from_enrichment(
  go_bp_file
)

go_cc_genes <- extract_genes_from_enrichment(
  go_cc_file
)

go_mf_genes <- extract_genes_from_enrichment(
  go_mf_file
)

kegg_genes <- extract_genes_from_enrichment(
  kegg_file
)


all_pathway_genes <- unique(
  c(
    go_bp_genes,
    go_cc_genes,
    go_mf_genes,
    kegg_genes
  )
)


cat("\n")
cat("====================================================\n")
cat("PATHWAY SUPPORT\n")
cat("====================================================\n")

cat(
  "GO-BP supported genes:",
  length(go_bp_genes),
  "\n"
)

cat(
  "GO-CC supported genes:",
  length(go_cc_genes),
  "\n"
)

cat(
  "GO-MF supported genes:",
  length(go_mf_genes),
  "\n"
)

cat(
  "KEGG supported genes:",
  length(kegg_genes),
  "\n"
)

cat(
  "Unique pathway-supported genes:",
  length(all_pathway_genes),
  "\n"
)


# ============================================================
# PART E
# PATHWAY SUPPORT FLAGS
# ============================================================

integrated$GO_BP_support <- (

  integrated$gene %in%
    go_bp_genes

)


integrated$GO_CC_support <- (

  integrated$gene %in%
    go_cc_genes

)


integrated$GO_MF_support <- (

  integrated$gene %in%
    go_mf_genes

)


integrated$KEGG_support <- (

  integrated$gene %in%
    kegg_genes

)


integrated$pathway_support <- (

  integrated$gene %in%
    all_pathway_genes

)


# ============================================================
# PART F
# BIOLOGICAL CATEGORY FLAGS
# ============================================================

integrated$epithelial_relevance <- (

  integrated$final_biological_category ==
    "Epithelial_Tumor_Relevant"

)


integrated$proliferation_relevance <- (

  integrated$final_biological_category ==
    "Proliferation_Relevant"

)


integrated$mitochondrial_relevance <- (

  integrated$final_biological_category ==
    "Mitochondrial_Metabolic"

)


integrated$stress_relevance <- (

  integrated$final_biological_category ==
    "Stress_Response"

)


integrated$housekeeping_relevance <- (

  integrated$final_biological_category ==
    "Housekeeping_BroadCellular"

)


# ============================================================
# PART G
# SAMPLE CONSISTENCY SCORE
# ============================================================

integrated$sample_consistency_score <- case_when(

  integrated$consistency_category ==
    "Consistent_3_of_3"
  ~ 1.00,

  integrated$consistency_category ==
    "Consistent_2_of_3"
  ~ 0.67,

  integrated$consistency_category ==
    "Sample_specific"
  ~ 0.33,

  TRUE
  ~ 0.00

)


# ============================================================
# PART H
# EXISTING BIOMARKER SCORE
# ============================================================

# Normalize the existing Step 21 score
# so it contributes consistently.

score_range <- range(
  integrated$adjusted_biomarker_score,
  na.rm = TRUE
)


if (
  diff(score_range) == 0
) {

  integrated$normalized_biomarker_score <- 1

} else {

  integrated$normalized_biomarker_score <- (

    integrated$adjusted_biomarker_score -
      score_range[1]

  ) / diff(score_range)

}


# ============================================================
# PART I
# BIOLOGICAL SUPPORT SCORE
# ============================================================

integrated$biological_support_score <- 0


# Epithelial/tumor relevance
integrated$biological_support_score <- (

  integrated$biological_support_score +

    ifelse(
      integrated$epithelial_relevance,
      0.25,
      0
    )

)


# Pathway support
integrated$biological_support_score <- (

  integrated$biological_support_score +

    ifelse(
      integrated$pathway_support,
      0.20,
      0
    )

)


# Proliferation relevance
integrated$biological_support_score <- (

  integrated$biological_support_score +

    ifelse(
      integrated$proliferation_relevance,
      0.15,
      0
    )

)


# Mitochondrial relevance
#
# This receives positive biological support but
# is deliberately lower than direct epithelial relevance.

integrated$biological_support_score <- (

  integrated$biological_support_score +

    ifelse(
      integrated$mitochondrial_relevance,
      0.10,
      0
    )

)


# Stress response
integrated$biological_support_score <- (

  integrated$biological_support_score +

    ifelse(
      integrated$stress_relevance,
      0.03,
      0
    )

)


# ============================================================
# PART J
# CONFIDENCE SCORE
# ============================================================

integrated$confidence_score <- (

  0.30 *
    integrated$normalized_biomarker_score

  +

  0.20 *
    integrated$sample_consistency_score

  +

  0.20 *
    integrated$prevalence_score

  +

  0.15 *
    integrated$specificity_score

  +

  0.15 *
    integrated$biological_support_score

)


# ============================================================
# PART K
# REVIEW PENALTY
# ============================================================

integrated$review_penalty <- ifelse(

  integrated$biological_review_flag ==
    "Review_required",

  0.15,

  0

)


integrated$final_integrated_score <- (

  integrated$confidence_score -
    integrated$review_penalty

)


# ============================================================
# PART L
# BIOLOGICAL TIER
# ============================================================

integrated$final_tier <- case_when(

  integrated$biological_review_flag ==
    "Review_required"
  ~
  "Tier_4_Biological_Review",

  integrated$epithelial_relevance &
    integrated$robust_3sample &
    integrated$final_integrated_score >= 0.65
  ~
  "Tier_1_Strong_Tumor_Relevant",

  integrated$robust_3sample &
    integrated$final_integrated_score >= 0.65
  ~
  "Tier_2_Strong_Robust_Candidate",

  integrated$robust_3sample &
    integrated$final_integrated_score >= 0.50
  ~
  "Tier_3_Exploratory_Robust",

  integrated$consistency_category ==
    "Consistent_2_of_3"
  ~
  "Tier_3_Moderate_Consistency",

  TRUE
  ~
  "Tier_4_Exploratory"

)


# ============================================================
# PART M
# METABOLIC FLAG
# ============================================================

integrated$metabolic_flag <- case_when(

  integrated$mitochondrial_relevance &
    integrated$pathway_support
  ~
  "OXPHOS_Pathway_Supported",

  integrated$mitochondrial_relevance
  ~
  "OXPHOS_Associated",

  TRUE
  ~
  "Not_OXPHOS_Classified"

)


# ============================================================
# PART N
# FINAL RANK
# ============================================================

integrated <- integrated %>%

  arrange(

    desc(final_integrated_score),

    desc(adjusted_biomarker_score),

    desc(absolute_log2FC)

  )


integrated$integrated_rank <- seq_len(
  nrow(integrated)
)


# ============================================================
# PART O
# SAVE COMPLETE TABLE
# ============================================================

final_columns <- c(

  "integrated_rank",
  "gene",

  "final_integrated_score",
  "final_tier",

  "adjusted_biomarker_score",
  "biomarker_score",

  "avg_log2FC",
  "absolute_log2FC",

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


final_table <- integrated %>%

  select(
    any_of(final_columns)
  )


write.csv(

  final_table,

  file.path(
    output_dir,
    "24_integrated_malignant_biomarker_ranking.csv"
  ),

  row.names = FALSE

)


# ============================================================
# PART P
# TOP 50
# ============================================================

top50 <- final_table %>%

  slice_head(
    n = 50
  )


write.csv(

  top50,

  file.path(
    output_dir,
    "24_top50_integrated_biomarker_candidates.csv"
  ),

  row.names = FALSE

)


# ============================================================
# PART Q
# TIER SUMMARY
# ============================================================

tier_summary <- final_table %>%

  count(

    final_tier,

    name = "genes"

  ) %>%

  arrange(
    desc(genes)
  )


cat("\n")
cat("====================================================\n")
cat("FINAL BIOMARKER TIERS\n")
cat("====================================================\n")

print(
  tier_summary
)


write.csv(

  tier_summary,

  file.path(
    output_dir,
    "24_biomarker_tier_summary.csv"
  ),

  row.names = FALSE

)


# ============================================================
# PART R
# BIOLOGICAL CATEGORY SUMMARY
# ============================================================

category_summary <- final_table %>%

  count(

    final_biological_category,

    name = "genes"

  ) %>%

  arrange(
    desc(genes)
  )


cat("\n")
cat("====================================================\n")
cat("BIOLOGICAL CATEGORY SUMMARY\n")
cat("====================================================\n")

print(
  category_summary
)


write.csv(

  category_summary,

  file.path(
    output_dir,
    "24_biological_category_summary.csv"
  ),

  row.names = FALSE

)


# ============================================================
# PART S
# HIGH-PRIORITY FINAL CANDIDATES
# ============================================================

high_priority_final <- final_table %>%

  filter(
    high_priority
  ) %>%

  arrange(
    integrated_rank
  )


write.csv(

  high_priority_final,

  file.path(
    output_dir,
    "24_high_priority_integrated_candidates.csv"
  ),

  row.names = FALSE

)


# ============================================================
# PART T
# ROBUST FINAL CANDIDATES
# ============================================================

robust_final <- final_table %>%

  filter(
    robust_3sample
  ) %>%

  arrange(
    integrated_rank
  )


write.csv(

  robust_final,

  file.path(
    output_dir,
    "24_robust_3sample_integrated_candidates.csv"
  ),

  row.names = FALSE

)


# ============================================================
# PART U
# TOP OXPHOS CANDIDATES
# ============================================================

top_oxphos <- final_table %>%

  filter(

    mitochondrial_relevance

  ) %>%

  arrange(
    integrated_rank
  )


write.csv(

  top_oxphos,

  file.path(
    output_dir,
    "24_OXPHOS_integrated_candidates.csv"
  ),

  row.names = FALSE

)


# ============================================================
# PART V
# TOP EPITHELIAL/TUMOR CANDIDATES
# ============================================================

top_epithelial <- final_table %>%

  filter(

    epithelial_relevance

  ) %>%

  arrange(
    integrated_rank
  )


write.csv(

  top_epithelial,

  file.path(
    output_dir,
    "24_epithelial_tumor_integrated_candidates.csv"
  ),

  row.names = FALSE

)


# ============================================================
# PART W
# TOP PROLIFERATION CANDIDATES
# ============================================================

top_proliferation <- final_table %>%

  filter(

    proliferation_relevance

  ) %>%

  arrange(
    integrated_rank
  )


write.csv(

  top_proliferation,

  file.path(
    output_dir,
    "24_proliferation_integrated_candidates.csv"
  ),

  row.names = FALSE

)


# ============================================================
# PART X
# TOP 20 DISPLAY
# ============================================================

cat("\n")
cat("====================================================\n")
cat("TOP 20 INTEGRATED CANDIDATES\n")
cat("====================================================\n")


print(

  final_table %>%

    select(

      integrated_rank,
      gene,
      final_integrated_score,
      final_tier,
      avg_log2FC,
      mean_sample_prevalence,
      consistency_category,
      final_biological_category,
      pathway_support,
      metabolic_flag

    ) %>%

    slice_head(
      n = 20
    )

)


# ============================================================
# PART Y
# SCORE DISTRIBUTION
# ============================================================

p1 <- ggplot(

  final_table,

  aes(
    x = final_integrated_score
  )

) +

  geom_histogram(
    bins = 30
  ) +

  labs(

    title =
      "Distribution of Integrated Biomarker Scores",

    x =
      "Integrated score",

    y =
      "Number of genes"

  ) +

  theme_classic()


ggsave(

  file.path(
    figure_dir,
    "24_integrated_score_distribution.pdf"
  ),

  p1,

  width = 9,

  height = 6

)


# ============================================================
# PART Z
# TIER VISUALIZATION
# ============================================================

p2 <- ggplot(

  tier_summary,

  aes(

    x =
      reorder(
        final_tier,
        genes
      ),

    y =
      genes

  )

) +

  geom_col() +

  coord_flip() +

  labs(

    title =
      "Integrated Biomarker Candidate Tiers",

    x =
      "Candidate tier",

    y =
      "Number of genes"

  ) +

  theme_classic()


ggsave(

  file.path(
    figure_dir,
    "24_biomarker_tiers.pdf"
  ),

  p2,

  width = 10,

  height = 6

)


# ============================================================
# PART AA
# SUMMARY FILE
# ============================================================

summary_file <- file.path(

  output_dir,

  "24_integrated_biomarker_prioritization_summary.txt"

)


summary_lines <- c(

  "INTEGRATED MALIGNANT BIOMARKER PRIORITIZATION SUMMARY",

  "====================================================",

  "",

  paste(
    "Step 21 prioritized candidates:",
    nrow(step21)
  ),

  paste(
    "Step 22 high-priority genes:",
    sum(step21$gene %in% integrated$gene[integrated$high_priority])
  ),

  paste(
    "Step 22 robust 3-sample genes:",
    sum(integrated$robust_3sample)
  ),

  paste(
    "Pathway-supported candidates:",
    sum(integrated$pathway_support)
  ),

  paste(
    "Epithelial/tumor-relevant candidates:",
    sum(integrated$epithelial_relevance)
  ),

  paste(
    "Mitochondrial/OXPHOS candidates:",
    sum(integrated$mitochondrial_relevance)
  ),

  paste(
    "Proliferation candidates:",
    sum(integrated$proliferation_relevance)
  ),

  paste(
    "Stress-response candidates:",
    sum(integrated$stress_relevance)
  ),

  paste(
    "Biological-review candidates:",
    sum(
      integrated$biological_review_flag ==
        "Review_required"
    )
  ),

  "",

  "Final tier distribution:",

  paste(

    capture.output(
      print(
        tier_summary
      )
    ),

    collapse = "\n"

  ),

  "",

  "Biological category distribution:",

  paste(

    capture.output(
      print(
        category_summary
      )
    ),

    collapse = "\n"

  ),

  "",

  "Interpretation:",

  "The integrated ranking combines existing Step 21",

  "biomarker scores with sample consistency and",

  "Step 23 biological annotations and Step 22",

  "pathway support.",

  "",

  "Pathway support is treated as biological evidence",

  "rather than proof of biomarker specificity.",

  "",

  "Mitochondrial/OXPHOS-associated genes are retained",

  "as biologically relevant candidates but are not",

  "automatically promoted above epithelial/tumor",

  "relevant candidates.",

  "",

  "Candidates requiring biological review are",

  "penalized in the integrated score and retained",

  "for transparency.",

  "",

  "All rankings are research prioritization results",

  "and are not validated clinical biomarkers."

)


writeLines(

  summary_lines,

  summary_file

)


# ============================================================
# FINAL
# ============================================================

cat("\n")
cat("====================================================\n")
cat("STEP 24 COMPLETED SUCCESSFULLY\n")
cat("====================================================\n")

cat(
  "Integrated candidates:",
  nrow(final_table),
  "\n"
)

cat(
  "Pathway-supported:",
  sum(integrated$pathway_support),
  "\n"
)

cat(
  "Epithelial/tumor relevant:",
  sum(integrated$epithelial_relevance),
  "\n"
)

cat(
  "Mitochondrial/OXPHOS:",
  sum(integrated$mitochondrial_relevance),
  "\n"
)

cat(
  "Robust 3-sample:",
  sum(integrated$robust_3sample),
  "\n"
)

cat("\n")
cat("Results:\n")
cat(output_dir, "\n")

cat("\n")
cat("====================================================\n")
cat("STEP 24 READY FOR REVIEW\n")
cat("====================================================\n")