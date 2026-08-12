# ============================================================
# STEP 23
# MALIGNANT CANDIDATE BIOLOGICAL VALIDATION
# ============================================================
#
# Project:
# SingleCell_BreastCancer_AI
#
# Dataset:
# GSE228499
#
# Purpose:
# Biological classification and validation of the candidate
# genes generated in Steps 21 and 22.
#
# IMPORTANT:
# Step 22 gene lists are treated as the authoritative
# candidate sets.
#
# This script does NOT:
#   - rerun Seurat
#   - rerun CopyKAT
#   - rerun differential expression
#   - change candidate thresholds
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
# 2. Project directories
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

output_dir <- file.path(
  project_dir,
  "results",
  "malignant",
  "biological_validation"
)

figure_dir <- file.path(
  project_dir,
  "figures",
  "malignant",
  "biological_validation"
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
cat("STEP 23: MALIGNANT CANDIDATE BIOLOGICAL VALIDATION\n")
cat("====================================================\n")


# ============================================================
# PART A
# LOAD STEP 21 PRIORITIZATION TABLE
# ============================================================

prioritized_file <- file.path(
  characterization_dir,
  "21_prioritized_malignant_biomarkers.csv"
)


if (!file.exists(prioritized_file)) {

  stop(
    "Missing Step 21 prioritized biomarker table:\n",
    prioritized_file
  )

}


prioritized <- read.csv(
  prioritized_file,
  stringsAsFactors = FALSE
)


cat("\n")
cat(
  "Step 21 prioritized genes:",
  nrow(prioritized),
  "\n"
)


# ============================================================
# PART B
# LOAD AUTHORITATIVE STEP 22 GENE LISTS
# ============================================================

high_priority_file <- file.path(
  pathway_dir,
  "22_high_priority_gene_list.txt"
)

robust_file <- file.path(
  pathway_dir,
  "22_robust_3sample_gene_list.txt"
)


if (!file.exists(high_priority_file)) {

  stop(
    "Missing Step 22 high-priority gene list:\n",
    high_priority_file
  )

}


if (!file.exists(robust_file)) {

  stop(
    "Missing Step 22 robust gene list:\n",
    robust_file
  )

}


high_priority_genes <- readLines(
  high_priority_file
)

high_priority_genes <- unique(
  trimws(high_priority_genes)
)

high_priority_genes <- high_priority_genes[
  high_priority_genes != "" &
  !is.na(high_priority_genes)
]


robust_genes <- readLines(
  robust_file
)

robust_genes <- unique(
  trimws(robust_genes)
)

robust_genes <- robust_genes[
  robust_genes != "" &
  !is.na(robust_genes)
]


cat("\n")
cat("====================================================\n")
cat("AUTHORITATIVE CANDIDATE SETS\n")
cat("====================================================\n")

cat(
  "High-priority genes:",
  length(high_priority_genes),
  "\n"
)

cat(
  "Robust 3-sample genes:",
  length(robust_genes),
  "\n"
)


# ============================================================
# PART C
# RECONCILE STEP 21 AND STEP 22
# ============================================================


step21_genes <- unique(
  prioritized$gene
)


missing_from_step22 <- setdiff(
  step21_genes,
  c(
    high_priority_genes,
    robust_genes
  )
)


high_priority_in_step21 <- intersect(
  high_priority_genes,
  step21_genes
)


robust_in_step21 <- intersect(
  robust_genes,
  step21_genes
)


cat("\n")
cat("Step 21 genes represented in Step 22:\n")

cat(
  "High-priority:",
  length(high_priority_in_step21),
  "\n"
)

cat(
  "Robust:",
  length(robust_in_step21),
  "\n"
)


write.csv(

  data.frame(
    gene = missing_from_step22
  ),

  file.path(
    output_dir,
    "23_genes_not_in_step22_candidate_sets.csv"
  ),

  row.names = FALSE

)


# ============================================================
# PART D
# BUILD ANNOTATION TABLE
# ============================================================


candidate_annotation <- data.frame(

  gene = unique(
    c(
      high_priority_genes,
      robust_genes
    )
  ),

  stringsAsFactors = FALSE

)


candidate_annotation$high_priority <- (

  candidate_annotation$gene %in%
    high_priority_genes

)


candidate_annotation$robust_3sample <- (

  candidate_annotation$gene %in%
    robust_genes

)


# ============================================================
# PART E
# Mitochondrial / OXPHOS GENES
# ============================================================

mitochondrial_genes <- c(

  "ATP5F1A",
  "ATP5F1B",
  "ATP5F1E",
  "ATP5MC1",
  "ATP5MC2",
  "ATP5MC3",
  "ATP5ME",
  "ATP5MPL",
  "ATP5PD",
  "ATP5PB",

  "COX4I1",
  "COX5A",
  "COX5B",
  "COX6A1",
  "COX6C",
  "COX7A2L",
  "COX7B",
  "COX8A",

  "CYC1",

  "NDUFA1",
  "NDUFA4",
  "NDUFA6",
  "NDUFA7",
  "NDUFA8",
  "NDUFA10",
  "NDUFA12",

  "NDUFAB1",

  "NDUFB1",
  "NDUFB4",
  "NDUFB5",
  "NDUFB6",
  "NDUFB7",
  "NDUFB8",
  "NDUFB9",
  "NDUFB10",

  "NDUFS1",
  "NDUFS2",
  "NDUFS3",
  "NDUFS6",
  "NDUFS7",

  "NDUFV1",

  "SDHB",

  "UQCC3",
  "UQCR10",
  "UQCRH",

  "ETFB",
  "ETFRF1",

  "MDH1",
  "MDH2",

  "IDH1",
  "IDH2",

  "ACO2",

  "GHITM",
  "CHCHD10"

)


# ============================================================
# PART F
# PROLIFERATION / CELL-CYCLE GENES
# ============================================================

proliferation_genes <- c(

  "MKI67",
  "TOP2A",
  "PCNA",
  "MCM2",
  "MCM3",
  "MCM4",
  "MCM5",
  "MCM6",
  "MCM7",
  "CCNA2",
  "CCNB1",
  "CCNB2",
  "CDC20",
  "CDC25C",
  "CDK1",
  "CDK2",
  "CDK4",
  "CDK6",
  "UBE2C",
  "BUB1",
  "BUB1B",
  "CENPF",
  "TYMS",
  "TK1",
  "HMGB2",
  "STMN1",
  "TUBA1B",
  "TUBB"

)


# ============================================================
# PART G
# EPITHELIAL / TUMOR-ASSOCIATED GENES
# ============================================================

epithelial_tumor_genes <- c(

  "EPCAM",
  "KRT8",
  "KRT18",
  "KRT19",

  "KRT5",
  "KRT14",
  "KRT17",

  "MUC1",
  "CEACAM6",

  "ERBB2",
  "ESR1",
  "PGR",

  "KRT23",
  "KRT80",

  "LAMC2",
  "MMP9",
  "MMP11",

  "SOX9",
  "AREG",

  "CD24",
  "EPN3",

  "TM4SF1",
  "PHLDA1",

  "PRSS23",
  "FCGRT"

)


# ============================================================
# PART H
# STRESS / IMMEDIATE-EARLY RESPONSE
# ============================================================

stress_genes <- c(

  "FOS",
  "JUN",
  "JUNB",
  "ATF3",
  "DUSP1",
  "DUSP2",
  "HSPA1A",
  "HSPA1B",
  "HSP90AA1",
  "DNAJB1"

)


# ============================================================
# PART I
# HOUSEKEEPING / BROAD CELLULAR
# ============================================================

housekeeping_genes <- c(

  "ACTB",
  "GAPDH",
  "B2M",
  "RPLP0",
  "RPL34",
  "RPL11",
  "RPS3",
  "RPS25",
  "RPL26",
  "MALAT1",
  "UBC",
  "EEF1A1"

)


# ============================================================
# PART J
# BIOLOGICAL CLASSIFICATION
# ============================================================

classify_gene <- function(gene) {

  classes <- character(0)


  if (
    gene %in%
    mitochondrial_genes
  ) {

    classes <- c(
      classes,
      "Mitochondrial_OXPHOS"
    )

  }


  if (
    gene %in%
    proliferation_genes
  ) {

    classes <- c(
      classes,
      "Proliferation_CellCycle"
    )

  }


  if (
    gene %in%
    epithelial_tumor_genes
  ) {

    classes <- c(
      classes,
      "Epithelial_TumorAssociated"
    )

  }


  if (
    gene %in%
    stress_genes
  ) {

    classes <- c(
      classes,
      "Stress_Response"
    )

  }


  if (
    gene %in%
    housekeeping_genes
  ) {

    classes <- c(
      classes,
      "Housekeeping_BroadCellular"
    )

  }


  if (
    length(classes) == 0
  ) {

    classes <- "Other_Candidate"

  }


  paste(
    classes,
    collapse = ";"
  )

}


candidate_annotation$biological_class <- vapply(

  candidate_annotation$gene,

  classify_gene,

  character(1)

)


# ============================================================
# PART K
# FINAL CATEGORY
# ============================================================

candidate_annotation$final_biological_category <- case_when(

  grepl(
    "Epithelial_TumorAssociated",
    candidate_annotation$biological_class
  )
  ~
  "Epithelial_Tumor_Relevant",

  grepl(
    "Proliferation_CellCycle",
    candidate_annotation$biological_class
  )
  ~
  "Proliferation_Relevant",

  grepl(
    "Mitochondrial_OXPHOS",
    candidate_annotation$biological_class
  )
  ~
  "Mitochondrial_Metabolic",

  grepl(
    "Stress_Response",
    candidate_annotation$biological_class
  )
  ~
  "Stress_Response",

  grepl(
    "Housekeeping_BroadCellular",
    candidate_annotation$biological_class
  )
  ~
  "Housekeeping_BroadCellular",

  TRUE
  ~
  "Other_Candidate"

)


# ============================================================
# PART L
# SAVE COMPLETE ANNOTATION
# ============================================================

write.csv(

  candidate_annotation,

  file.path(
    output_dir,
    "23_all_candidate_biological_annotation.csv"
  ),

  row.names = FALSE

)


# ============================================================
# PART M
# HIGH-PRIORITY ANNOTATION
# ============================================================

high_priority_annotation <- candidate_annotation %>%

  filter(
    high_priority
  )


write.csv(

  high_priority_annotation,

  file.path(
    output_dir,
    "23_high_priority_biological_annotation.csv"
  ),

  row.names = FALSE

)


# ============================================================
# PART N
# ROBUST ANNOTATION
# ============================================================

robust_annotation <- candidate_annotation %>%

  filter(
    robust_3sample
  )


write.csv(

  robust_annotation,

  file.path(
    output_dir,
    "23_robust_3sample_biological_annotation.csv"
  ),

  row.names = FALSE

)


# ============================================================
# PART O
# CATEGORY SUMMARY
# ============================================================

high_priority_category_summary <- high_priority_annotation %>%

  count(

    final_biological_category,

    name = "genes"

  ) %>%

  arrange(
    desc(genes)
  )


robust_category_summary <- robust_annotation %>%

  count(

    final_biological_category,

    name = "genes"

  ) %>%

  arrange(
    desc(genes)
  )


cat("\n")
cat("====================================================\n")
cat("HIGH-PRIORITY BIOLOGICAL CATEGORIES\n")
cat("====================================================\n")

print(
  high_priority_category_summary
)


cat("\n")
cat("====================================================\n")
cat("ROBUST 3-SAMPLE BIOLOGICAL CATEGORIES\n")
cat("====================================================\n")

print(
  robust_category_summary
)


write.csv(

  high_priority_category_summary,

  file.path(
    output_dir,
    "23_high_priority_category_summary.csv"
  ),

  row.names = FALSE

)


write.csv(

  robust_category_summary,

  file.path(
    output_dir,
    "23_robust_category_summary.csv"
  ),

  row.names = FALSE

)


# ============================================================
# PART P
# SPECIFIC CANDIDATE GROUPS
# ============================================================

oxphos_candidates <- candidate_annotation %>%

  filter(

    grepl(
      "Mitochondrial_OXPHOS",
      biological_class
    )

  )


proliferation_candidates <- candidate_annotation %>%

  filter(

    grepl(
      "Proliferation_CellCycle",
      biological_class
    )

  )


epithelial_candidates <- candidate_annotation %>%

  filter(

    grepl(
      "Epithelial_TumorAssociated",
      biological_class
    )

  )


stress_candidates <- candidate_annotation %>%

  filter(

    grepl(
      "Stress_Response",
      biological_class
    )

  )


housekeeping_candidates <- candidate_annotation %>%

  filter(

    grepl(
      "Housekeeping_BroadCellular",
      biological_class
    )

  )


write.csv(

  oxphos_candidates,

  file.path(
    output_dir,
    "23_OXPHOS_candidates.csv"
  ),

  row.names = FALSE

)


write.csv(

  proliferation_candidates,

  file.path(
    output_dir,
    "23_proliferation_candidates.csv"
  ),

  row.names = FALSE

)


write.csv(

  epithelial_candidates,

  file.path(
    output_dir,
    "23_epithelial_tumor_candidates.csv"
  ),

  row.names = FALSE

)


write.csv(

  stress_candidates,

  file.path(
    output_dir,
    "23_stress_response_candidates.csv"
  ),

  row.names = FALSE

)


write.csv(

  housekeeping_candidates,

  file.path(
    output_dir,
    "23_housekeeping_candidates.csv"
  ),

  row.names = FALSE

)


# ============================================================
# PART Q
# STEP 22 PATHWAY VALIDATION
# ============================================================

go_bp_file <- file.path(

  pathway_dir,

  "22_robust_3sample_GO_BP.csv"

)


mitochondrial_terms <- data.frame()


if (
  file.exists(go_bp_file)
) {

  go_bp <- read.csv(

    go_bp_file,

    stringsAsFactors = FALSE

  )


  if (
    "Description" %in%
    colnames(go_bp)
  ) {

    mitochondrial_terms <- go_bp %>%

      filter(

        grepl(

          "respiration|oxidative phosphorylation|electron transport|ATP synthesis|ATP metabolic|energy",

          Description,

          ignore.case = TRUE

        )

      )

  }

}


cat("\n")
cat(
  "Mitochondrial/energy GO-BP terms:",
  nrow(mitochondrial_terms),
  "\n"
)


write.csv(

  mitochondrial_terms,

  file.path(
    output_dir,
    "23_mitochondrial_energy_GO_BP_terms.csv"
  ),

  row.names = FALSE

)


# ============================================================
# PART R
# BIOLOGICAL REVIEW FLAGS
# ============================================================

review_candidates <- candidate_annotation %>%

  filter(

    grepl(

      "Stress_Response|Housekeeping_BroadCellular",

      biological_class

    )

  )


write.csv(

  review_candidates,

  file.path(
    output_dir,
    "23_candidates_requiring_biological_review.csv"
  ),

  row.names = FALSE

)


# ============================================================
# PART S
# VISUALIZATION
# ============================================================

plot_data <- high_priority_category_summary


p <- ggplot(

  plot_data,

  aes(

    x =
      reorder(
        final_biological_category,
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
      "Biological Classification of High-Priority Candidates",

    x =
      "Biological category",

    y =
      "Number of genes"

  ) +

  theme_classic()


ggsave(

  file.path(
    figure_dir,
    "23_high_priority_biological_categories.pdf"
  ),

  p,

  width = 9,

  height = 6

)


# ============================================================
# PART T
# SUMMARY
# ============================================================

summary_file <- file.path(

  output_dir,

  "23_malignant_candidate_biological_validation_summary.txt"

)


summary_lines <- c(

  "MALIGNANT CANDIDATE BIOLOGICAL VALIDATION SUMMARY",

  "====================================================",

  "",

  paste(
    "Step 21 prioritized genes:",
    nrow(prioritized)
  ),

  paste(
    "Step 22 high-priority genes:",
    length(high_priority_genes)
  ),

  paste(
    "Step 22 robust 3-sample genes:",
    length(robust_genes)
  ),

  paste(
    "Step 21 genes absent from Step 22 candidate sets:",
    length(missing_from_step22)
  ),

  "",

  "High-priority biological categories:",

  paste(

    capture.output(
      print(
        high_priority_category_summary
      )
    ),

    collapse = "\n"

  ),

  "",

  "Robust 3-sample biological categories:",

  paste(

    capture.output(
      print(
        robust_category_summary
      )
    ),

    collapse = "\n"

  ),

  "",

  paste(
    "Mitochondrial/OXPHOS candidates:",
    nrow(oxphos_candidates)
  ),

  paste(
    "Proliferation candidates:",
    nrow(proliferation_candidates)
  ),

  paste(
    "Epithelial/tumor-associated candidates:",
    nrow(epithelial_candidates)
  ),

  paste(
    "Stress-response candidates:",
    nrow(stress_candidates)
  ),

  paste(
    "Housekeeping candidates:",
    nrow(housekeeping_candidates)
  ),

  paste(
    "Candidates requiring biological review:",
    nrow(review_candidates)
  ),

  "",

  paste(
    "Significant mitochondrial/energy GO-BP terms:",
    nrow(mitochondrial_terms)
  ),

  "",

  "Interpretation:",

  "The robust 3-sample candidate set contains a",

  "strong mitochondrial/oxidative-phosphorylation",

  "biological signal based on GO enrichment.",

  "",

  "This represents an association and does not",

  "establish a cancer-specific mechanism.",

  "",

  "Mitochondrial/metabolic candidates should be",

  "considered separately from epithelial/tumor",

  "associated candidates during final biomarker",

  "selection.",

  "",

  "Stress and housekeeping genes are retained",

  "for transparency but flagged for biological",

  "review.",

  "",

  "All classifications are research annotations",

  "and are not clinical diagnostic classifications."

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
cat("STEP 23 COMPLETED SUCCESSFULLY\n")
cat("====================================================\n")

cat(
  "High-priority genes:",
  length(high_priority_genes),
  "\n"
)

cat(
  "Robust 3-sample genes:",
  length(robust_genes),
  "\n"
)

cat(
  "OXPHOS candidates:",
  nrow(oxphos_candidates),
  "\n"
)

cat(
  "Proliferation candidates:",
  nrow(proliferation_candidates),
  "\n"
)

cat(
  "Epithelial/tumor candidates:",
  nrow(epithelial_candidates),
  "\n"
)

cat(
  "Biological-review candidates:",
  nrow(review_candidates),
  "\n"
)

cat(
  "Mitochondrial/energy GO-BP terms:",
  nrow(mitochondrial_terms),
  "\n"
)

cat("\n")
cat("Results directory:\n")
cat(output_dir, "\n")

cat("\n")
cat("====================================================\n")
cat("READY FOR NEXT INTEGRATION STEP\n")
cat("====================================================\n")