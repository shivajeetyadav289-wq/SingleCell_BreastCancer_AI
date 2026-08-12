# ============================================================
# STEP 22: MALIGNANT PATHWAY ENRICHMENT
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
# Functional and pathway characterization of prioritized
# malignant-cell biomarker candidates.
#
# Gene sets:
#
#   A. 164 high-priority candidates
#   B. 626 robust 3-sample candidates
#
# Analyses:
#
#   1. GO Biological Process
#   2. GO Molecular Function
#   3. GO Cellular Component
#   4. KEGG
#   5. Reactome
#
# IMPORTANT:
# Pathway enrichment identifies biological processes associated
# with the candidate gene sets. It does NOT establish causality
# or clinical biomarker validity.
# ============================================================


# ------------------------------------------------------------
# 1. Required packages
# ------------------------------------------------------------

suppressPackageStartupMessages({

  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)

})


# ------------------------------------------------------------
# 2. Optional Reactome package
# ------------------------------------------------------------

reactome_available <-

  requireNamespace(
    "ReactomePA",
    quietly = TRUE
  )


# ------------------------------------------------------------
# 3. Project directories
# ------------------------------------------------------------

project_dir <- path.expand(
  "~/Projects/SingleCell_BreastCancer_AI"
)

input_dir <- file.path(
  project_dir,
  "results",
  "malignant",
  "characterization"
)

output_dir <- file.path(
  project_dir,
  "results",
  "malignant",
  "pathway_enrichment"
)

figure_dir <- file.path(
  project_dir,
  "figures",
  "malignant",
  "pathway_enrichment"
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
# 4. Input files
# ------------------------------------------------------------

high_priority_file <- file.path(

  input_dir,

  "21_top_20_malignant_biomarker_candidates.csv"

)

# NOTE:
# The complete high-priority set is needed here.
# Step 21 currently saves the 164 candidates in
# 21_prioritized_malignant_biomarkers.csv.
#
# We therefore reconstruct the 164-gene set from that file.

prioritized_file <- file.path(

  input_dir,

  "21_prioritized_malignant_biomarkers.csv"

)


robust_file <- file.path(

  input_dir,

  "21_robust_3sample_biomarker_candidates.csv"

)


# ------------------------------------------------------------
# 5. Header
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("STEP 22: MALIGNANT PATHWAY ENRICHMENT\n")
cat("====================================================\n")


# ------------------------------------------------------------
# 6. Check files
# ------------------------------------------------------------

if (!file.exists(prioritized_file)) {

  stop(
    "\nMissing Step 21 prioritized candidate file:\n",
    prioritized_file,
    "\n"
  )

}


if (!file.exists(robust_file)) {

  stop(
    "\nMissing Step 21 robust candidate file:\n",
    robust_file,
    "\n"
  )

}


# ------------------------------------------------------------
# 7. Load candidate tables
# ------------------------------------------------------------

cat("\n")
cat("Loading Step 21 candidate tables...\n")


prioritized <- read.csv(

  prioritized_file,

  stringsAsFactors = FALSE

)


robust <- read.csv(

  robust_file,

  stringsAsFactors = FALSE

)


cat(
  "Prioritized genes loaded:",
  nrow(prioritized),
  "\n"
)


cat(
  "Robust 3-sample genes loaded:",
  nrow(robust),
  "\n"
)


# ------------------------------------------------------------
# 8. Select high-priority candidates
# ------------------------------------------------------------
#
# Step 21 defines:
#
# p_val_adj < 0.01
# avg_log2FC >= 1
# pct.1 >= 0.40
# prevalence difference >= 0.20
#
# We reconstruct this exact high-priority set from the
# complete prioritized table.
# ------------------------------------------------------------

high_priority <-

  prioritized %>%

  filter(

    p_val_adj < 0.01,

    avg_log2FC >= 1.0,

    pct.1 >= 0.40,

    prevalence_difference >= 0.20

  )


# Remove genes explicitly flagged for biological review.

high_priority_clean <-

  high_priority %>%

  filter(

    biological_review_flag !=
      "Review_required"

  )


# ------------------------------------------------------------
# 9. Extract gene symbols
# ------------------------------------------------------------

high_priority_genes <-

  unique(

    high_priority_clean$gene

  )


robust_genes <-

  unique(

    robust$gene

  )


cat("\n")
cat("====================================================\n")
cat("GENE SETS\n")
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


# ------------------------------------------------------------
# 10. Save gene lists
# ------------------------------------------------------------

writeLines(

  high_priority_genes,

  file.path(

    output_dir,

    "22_high_priority_gene_list.txt"

  )

)


writeLines(

  robust_genes,

  file.path(

    output_dir,

    "22_robust_3sample_gene_list.txt"

  )

)


# ============================================================
# PART A
# GENE SYMBOL → ENTREZ ID
# ============================================================


# ------------------------------------------------------------
# 11. Convert gene symbols
# ------------------------------------------------------------

convert_symbols <- function(

  genes

) {

  genes <- unique(

    genes[

      !is.na(genes) &
      genes != ""

    ]

  )


  mapping <-

    suppressMessages(

      bitr(

        genes,

        fromType =
          "SYMBOL",

        toType =
          c(
            "ENTREZID",
            "ENSEMBL"
          ),

        OrgDb =
          org.Hs.eg.db

      )

    )


  mapping <-

    mapping %>%

    distinct(

      SYMBOL,

      .keep_all = TRUE

    )


  return(mapping)

}


# ------------------------------------------------------------
# 12. Convert both gene sets
# ------------------------------------------------------------

cat("\n")
cat("Converting gene symbols to Entrez IDs...\n")


high_mapping <-

  convert_symbols(

    high_priority_genes

  )


robust_mapping <-

  convert_symbols(

    robust_genes

  )


cat(
  "High-priority genes mapped:",
  nrow(high_mapping),
  "\n"
)


cat(
  "Robust genes mapped:",
  nrow(robust_mapping),
  "\n"
)


# ------------------------------------------------------------
# 13. Save mapping
# ------------------------------------------------------------

write.csv(

  high_mapping,

  file.path(

    output_dir,

    "22_high_priority_gene_ID_mapping.csv"

  ),

  row.names = FALSE

)


write.csv(

  robust_mapping,

  file.path(

    output_dir,

    "22_robust_gene_ID_mapping.csv"

  ),

  row.names = FALSE

)


# ============================================================
# PART B
# BACKGROUND GENE SET
# ============================================================
#
# The background should ideally represent genes detectable
# in the experiment.
#
# Step 20 tested 5,637 genes.
#
# Therefore we use the genes from the Step 20 DE table as the
# experimental background rather than the entire human genome.
# ============================================================


cat("\n")
cat("Preparing experimental background...\n")


background_genes <-

  unique(

    prioritized$gene

  )


background_mapping <-

  convert_symbols(

    background_genes

  )


background_entrez <-

  unique(

    background_mapping$ENTREZID

  )


cat(
  "Background genes mapped:",
  length(background_entrez),
  "\n"
)


# ============================================================
# PART C
# ENRICHMENT FUNCTION
# ============================================================


run_GO_enrichment <- function(

  entrez_genes,

  background,

  dataset_name

) {


  cat("\n")
  cat("---------------------------------------------\n")
  cat(
    "GO enrichment:",
    dataset_name,
    "\n"
  )
  cat("---------------------------------------------\n")


  if (

    length(entrez_genes) < 5

  ) {

    warning(
      "Too few genes for enrichment: ",
      dataset_name
    )

    return(NULL)

  }


  # Biological Process

  ego_BP <-

    enrichGO(

      gene =
        entrez_genes,

      universe =
        background,

      OrgDb =
        org.Hs.eg.db,

      keyType =
        "ENTREZID",

      ont =
        "BP",

      pAdjustMethod =
        "BH",

      pvalueCutoff =
        0.05,

      qvalueCutoff =
        0.20,

      readable =
        TRUE

    )


  # Molecular Function

  ego_MF <-

    enrichGO(

      gene =
        entrez_genes,

      universe =
        background,

      OrgDb =
        org.Hs.eg.db,

      keyType =
        "ENTREZID",

      ont =
        "MF",

      pAdjustMethod =
        "BH",

      pvalueCutoff =
        0.05,

      qvalueCutoff =
        0.20,

      readable =
        TRUE

    )


  # Cellular Component

  ego_CC <-

    enrichGO(

      gene =
        entrez_genes,

      universe =
        background,

      OrgDb =
        org.Hs.eg.db,

      keyType =
        "ENTREZID",

      ont =
        "CC",

      pAdjustMethod =
        "BH",

      pvalueCutoff =
        0.05,

      qvalueCutoff =
        0.20,

      readable =
        TRUE

    )


  # Save results

  if (!is.null(ego_BP)) {

    write.csv(

      as.data.frame(ego_BP),

      file.path(

        output_dir,

        paste0(
          "22_",
          dataset_name,
          "_GO_BP.csv"
        )

      ),

      row.names = FALSE

    )

  }


  if (!is.null(ego_MF)) {

    write.csv(

      as.data.frame(ego_MF),

      file.path(

        output_dir,

        paste0(
          "22_",
          dataset_name,
          "_GO_MF.csv"
        )

      ),

      row.names = FALSE

    )

  }


  if (!is.null(ego_CC)) {

    write.csv(

      as.data.frame(ego_CC),

      file.path(

        output_dir,

        paste0(
          "22_",
          dataset_name,
          "_GO_CC.csv"
        )

      ),

      row.names = FALSE

    )

  }


  return(

    list(

      BP = ego_BP,

      MF = ego_MF,

      CC = ego_CC

    )

  )

}


# ============================================================
# PART D
# KEGG FUNCTION
# ============================================================


run_KEGG <- function(

  entrez_genes,

  background,

  dataset_name

) {


  cat("\n")
  cat(
    "KEGG enrichment:",
    dataset_name,
    "\n"
  )


  if (

    length(entrez_genes) < 5

  ) {

    return(NULL)

  }


  ekegg <-

    enrichKEGG(

      gene =
        entrez_genes,

      universe =
        background,

      organism =
        "hsa",

      keyType =
        "ncbi-geneid",

      pAdjustMethod =
        "BH",

      pvalueCutoff =
        0.05,

      qvalueCutoff =
        0.20

    )


  if (!is.null(ekegg)) {

    write.csv(

      as.data.frame(ekegg),

      file.path(

        output_dir,

        paste0(
          "22_",
          dataset_name,
          "_KEGG.csv"
        )

      ),

      row.names = FALSE

    )

  }


  return(ekegg)

}


# ============================================================
# PART E
# REACTOME FUNCTION
# ============================================================


run_Reactome <- function(

  entrez_genes,

  background,

  dataset_name

) {


  if (!reactome_available) {

    cat(
      "ReactomePA not installed; skipping Reactome.\n"
    )

    return(NULL)

  }


  cat("\n")
  cat(
    "Reactome enrichment:",
    dataset_name,
    "\n"
  )


  result <-

    ReactomePA::enrichPathway(

      gene =
        entrez_genes,

      universe =
        background,

      organism =
        "human",

      pAdjustMethod =
        "BH",

      pvalueCutoff =
        0.05,

      qvalueCutoff =
        0.20,

      readable =
        TRUE

    )


  if (!is.null(result)) {

    write.csv(

      as.data.frame(result),

      file.path(

        output_dir,

        paste0(
          "22_",
          dataset_name,
          "_Reactome.csv"
        )

      ),

      row.names = FALSE

    )

  }


  return(result)

}


# ============================================================
# PART F
# RUN ENRICHMENT
# ============================================================


high_entrez <-

  unique(

    high_mapping$ENTREZID

  )


robust_entrez <-

  unique(

    robust_mapping$ENTREZID

  )


high_GO <-

  run_GO_enrichment(

    high_entrez,

    background_entrez,

    "high_priority"

  )


robust_GO <-

  run_GO_enrichment(

    robust_entrez,

    background_entrez,

    "robust_3sample"

  )


high_KEGG <-

  run_KEGG(

    high_entrez,

    background_entrez,

    "high_priority"

  )


robust_KEGG <-

  run_KEGG(

    robust_entrez,

    background_entrez,

    "robust_3sample"

  )


high_Reactome <-

  run_Reactome(

    high_entrez,

    background_entrez,

    "high_priority"

  )


robust_Reactome <-

  run_Reactome(

    robust_entrez,

    background_entrez,

    "robust_3sample"

  )


# ============================================================
# PART G
# TOP PATHWAY SUMMARY
# ============================================================


extract_top_terms <- function(

  enrichment_object,

  n = 20

) {


  if (

    is.null(enrichment_object)

  ) {

    return(
      data.frame()
    )

  }


  df <-

    as.data.frame(

      enrichment_object

    )


  if (

    nrow(df) == 0

  ) {

    return(
      data.frame()
    )

  }


  df %>%

    arrange(

      p.adjust

    ) %>%

    head(n)

}


# ------------------------------------------------------------
# 14. High-priority GO summaries
# ------------------------------------------------------------

if (!is.null(high_GO)) {


  high_BP_top <-

    extract_top_terms(

      high_GO$BP

    )


  high_MF_top <-

    extract_top_terms(

      high_GO$MF

    )


  high_CC_top <-

    extract_top_terms(

      high_GO$CC

    )


  write.csv(

    high_BP_top,

    file.path(

      output_dir,

      "22_high_priority_top_GO_BP.csv"

    ),

    row.names = FALSE

  )


  write.csv(

    high_MF_top,

    file.path(

      output_dir,

      "22_high_priority_top_GO_MF.csv"

    ),

    row.names = FALSE

  )


  write.csv(

    high_CC_top,

    file.path(

      output_dir,

      "22_high_priority_top_GO_CC.csv"

    ),

    row.names = FALSE

  )

}


# ------------------------------------------------------------
# 15. Robust GO summaries
# ------------------------------------------------------------

if (!is.null(robust_GO)) {


  robust_BP_top <-

    extract_top_terms(

      robust_GO$BP

    )


  robust_MF_top <-

    extract_top_terms(

      robust_GO$MF

    )


  robust_CC_top <-

    extract_top_terms(

      robust_GO$CC

    )


  write.csv(

    robust_BP_top,

    file.path(

      output_dir,

      "22_robust_top_GO_BP.csv"

    ),

    row.names = FALSE

  )


  write.csv(

    robust_MF_top,

    file.path(

      output_dir,

      "22_robust_top_GO_MF.csv"

    ),

    row.names = FALSE

  )


  write.csv(

    robust_CC_top,

    file.path(

      output_dir,

      "22_robust_top_GO_CC.csv"

    ),

    row.names = FALSE

  )

}


# ============================================================
# PART H
# DOTPLOTS
# ============================================================


make_dotplot <- function(

  enrichment_object,

  filename,

  title,

  n = 15

) {


  if (

    is.null(enrichment_object)

  ) {

    return(NULL)

  }


  df <-

    as.data.frame(

      enrichment_object

    )


  if (

    nrow(df) == 0

  ) {

    cat(
      "No significant terms for:",
      title,
      "\n"
    )

    return(NULL)

  }


  n_plot <-

    min(
      n,
      nrow(df)
    )


  p <-

    dotplot(

      enrichment_object,

      showCategory =
        n_plot

    ) +

    ggtitle(

      title

    ) +

    theme_classic()


  ggsave(

    filename,

    p,

    width =
      10,

    height =
      7

  )

}


# ------------------------------------------------------------
# 16. High-priority plots
# ------------------------------------------------------------

if (!is.null(high_GO)) {


  make_dotplot(

    high_GO$BP,

    file.path(

      figure_dir,

      "22_high_priority_GO_BP_dotplot.pdf"

    ),

    "High-Priority Candidates: GO Biological Process"

  )


  make_dotplot(

    high_GO$MF,

    file.path(

      figure_dir,

      "22_high_priority_GO_MF_dotplot.pdf"

    ),

    "High-Priority Candidates: GO Molecular Function"

  )


  make_dotplot(

    high_GO$CC,

    file.path(

      figure_dir,

      "22_high_priority_GO_CC_dotplot.pdf"

    ),

    "High-Priority Candidates: GO Cellular Component"

  )

}


# ------------------------------------------------------------
# 17. Robust plots
# ------------------------------------------------------------

if (!is.null(robust_GO)) {


  make_dotplot(

    robust_GO$BP,

    file.path(

      figure_dir,

      "22_robust_GO_BP_dotplot.pdf"

    ),

    "Robust 3-Sample Candidates: GO Biological Process"

  )


  make_dotplot(

    robust_GO$MF,

    file.path(

      figure_dir,

      "22_robust_GO_MF_dotplot.pdf"

    ),

    "Robust 3-Sample Candidates: GO Molecular Function"

  )


  make_dotplot(

    robust_GO$CC,

    file.path(

      figure_dir,

      "22_robust_GO_CC_dotplot.pdf"

    ),

    "Robust 3-Sample Candidates: GO Cellular Component"

  )

}


# ============================================================
# PART I
# KEGG DOTPLOTS
# ============================================================


if (!is.null(high_KEGG)) {

  make_dotplot(

    high_KEGG,

    file.path(

      figure_dir,

      "22_high_priority_KEGG_dotplot.pdf"

    ),

    "High-Priority Candidates: KEGG Pathways"

  )

}


if (!is.null(robust_KEGG)) {

  make_dotplot(

    robust_KEGG,

    file.path(

      figure_dir,

      "22_robust_KEGG_dotplot.pdf"

    ),

    "Robust 3-Sample Candidates: KEGG Pathways"

  )

}


# ============================================================
# PART J
# PATHWAY OVERLAP
# ============================================================


get_term_ids <- function(

  enrichment_object

) {


  if (

    is.null(enrichment_object)

  ) {

    return(
      character(0)
    )

  }


  df <-

    as.data.frame(

      enrichment_object

    )


  if (

    nrow(df) == 0

  ) {

    return(
      character(0)
    )

  }


  df %>%

    filter(

      p.adjust < 0.05

    ) %>%

    pull(

      Description

    ) %>%

    unique()

}


high_BP_terms <-

  get_term_ids(

    if (
      !is.null(high_GO)
    )
      high_GO$BP
    else
      NULL

  )


robust_BP_terms <-

  get_term_ids(

    if (
      !is.null(robust_GO)
    )
      robust_GO$BP
    else
      NULL

  )


shared_BP_terms <-

  intersect(

    high_BP_terms,

    robust_BP_terms

  )


write.csv(

  data.frame(

    shared_GO_BP_terms =
      shared_BP_terms

  ),

  file.path(

    output_dir,

    "22_shared_GO_BP_terms.csv"

  ),

  row.names = FALSE

)


# ============================================================
# PART K
# SUMMARY
# ============================================================


count_terms <- function(

  enrichment_object

) {


  if (

    is.null(enrichment_object)

  ) {

    return(0)

  }


  df <-

    as.data.frame(

      enrichment_object

    )


  if (

    nrow(df) == 0

  ) {

    return(0)

  }


  sum(

    df$p.adjust < 0.05,

    na.rm = TRUE

  )

}


high_BP_count <-

  count_terms(

    if (
      !is.null(high_GO)
    )
      high_GO$BP
    else
      NULL

  )


high_MF_count <-

  count_terms(

    if (
      !is.null(high_GO)
    )
      high_GO$MF
    else
      NULL

  )


high_CC_count <-

  count_terms(

    if (
      !is.null(high_GO)
    )
      high_GO$CC
    else
      NULL

  )


robust_BP_count <-

  count_terms(

    if (
      !is.null(robust_GO)
    )
      robust_GO$BP
    else
      NULL

  )


robust_MF_count <-

  count_terms(

    if (
      !is.null(robust_GO)
    )
      robust_GO$MF
    else
      NULL

  )


robust_CC_count <-

  count_terms(

    if (
      !is.null(robust_GO)
    )
      robust_GO$CC
    else
      NULL

  )


# ------------------------------------------------------------
# 18. Summary report
# ------------------------------------------------------------

summary_file <- file.path(

  output_dir,

  "22_pathway_enrichment_summary.txt"

)


summary_text <- c(

  "MALIGNANT PATHWAY ENRICHMENT SUMMARY",

  "====================================",

  "",

  paste(
    "High-priority genes:",
    length(high_priority_genes)
  ),

  paste(
    "High-priority genes mapped to Entrez:",
    length(high_entrez)
  ),

  "",

  paste(
    "Robust 3-sample genes:",
    length(robust_genes)
  ),

  paste(
    "Robust genes mapped to Entrez:",
    length(robust_entrez)
  ),

  "",

  "GO Biological Process:",

  paste(
    "High-priority significant terms:",
    high_BP_count
  ),

  paste(
    "Robust significant terms:",
    robust_BP_count
  ),

  "",

  "GO Molecular Function:",

  paste(
    "High-priority significant terms:",
    high_MF_count
  ),

  paste(
    "Robust significant terms:",
    robust_MF_count
  ),

  "",

  "GO Cellular Component:",

  paste(
    "High-priority significant terms:",
    high_CC_count
  ),

  paste(
    "Robust significant terms:",
    robust_CC_count
  ),

  "",

  paste(
    "Shared significant GO-BP terms:",
    length(shared_BP_terms)
  ),

  "",

  paste(
    "ReactomePA available:",
    reactome_available
  ),

  "",

  "Interpretation:",

  "Pathway enrichment identifies biological",

  "processes associated with prioritized",

  "malignant-cell candidates.",

  "",

  "Enrichment does not establish causality",

  "and does not constitute clinical validation.",

  "",

  "The robust 3-sample gene set is particularly",

  "important because these candidates show",

  "expression prevalence across BC11, BC12",

  "and BC17."

)


writeLines(

  summary_text,

  summary_file

)


# ============================================================
# FINAL
# ============================================================


cat("\n")
cat("====================================================\n")
cat("STEP 22 COMPLETED\n")
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
  "High-priority GO-BP terms:",
  high_BP_count,
  "\n"
)


cat(
  "Robust GO-BP terms:",
  robust_BP_count,
  "\n"
)


cat(
  "Shared significant GO-BP terms:",
  length(shared_BP_terms),
  "\n"
)


cat("\n")
cat("Results directory:\n")
cat(output_dir, "\n")


cat("\n")
cat("====================================================\n")
cat("READY FOR BIOLOGICAL VALIDATION\n")
cat("====================================================\n")