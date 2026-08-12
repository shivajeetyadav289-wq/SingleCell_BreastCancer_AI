#!/usr/bin/env Rscript

# ============================================================
# STEP 27
# EXTERNAL EXPRESSION VALIDATION
# ============================================================
#
# Purpose:
#   Validate expression of the Step 25 final candidate panel
#   in the independent GSE176078 breast-cancer scRNA-seq dataset.
#
# Input:
#   GSE176078 sparse MatrixMarket matrix
#   GSE176078 metadata
#   Step 25 final research candidate panel
#
# Strategy:
#   1. Read the 50 Step-25 candidates
#   2. Match candidates to external genes
#   3. Restrict analysis to Cancer Epithelial cells
#   4. Extract only candidate-gene rows from sparse matrix
#   5. Calculate cell-level expression prevalence
#   6. Calculate sample-level expression prevalence
#   7. Rank externally supported candidates
#
# ============================================================

options(stringsAsFactors = FALSE)

cat("\n====================================================\n")
cat("STEP 27: EXTERNAL EXPRESSION VALIDATION\n")
cat("====================================================\n\n")

# ------------------------------------------------------------
# 1. INPUT PATHS
# ------------------------------------------------------------

candidate_file <-
  "results/malignant/final_candidates/25_final_research_candidate_panel.csv"

base_dir <-
  "data/external/GSE176078/Wu_etal_2021_BRCA_scRNASeq"

# Correct actual directory name if needed
if (!dir.exists(base_dir)) {
  base_dir <-
    "data/external/GSE176078/Wu_etal_2021_BRCA_scRNASeq"
}

matrix_file <- file.path(base_dir, "count_matrix_sparse.mtx")
genes_file  <- file.path(base_dir, "count_matrix_genes.tsv")
barcodes_file <- file.path(base_dir, "count_matrix_barcodes.tsv")
metadata_file <- file.path(base_dir, "metadata.csv")

output_dir <- "results/malignant/external_validation"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2. CHECK FILES
# ------------------------------------------------------------

required_files <- c(
  candidate_file,
  matrix_file,
  genes_file,
  barcodes_file,
  metadata_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  cat("ERROR: Missing files:\n")
  print(missing_files)
  quit(status = 1)
}

cat("Input files found.\n\n")

# ------------------------------------------------------------
# 3. LOAD REQUIRED PACKAGES
# ------------------------------------------------------------

if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Package 'Matrix' is required but not installed.")
}

library(Matrix)

# ------------------------------------------------------------
# 4. LOAD STEP 25 CANDIDATES
# ------------------------------------------------------------

cat("Loading Step 25 candidate panel...\n")

candidates <- read.csv(
  candidate_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!"gene" %in% colnames(candidates)) {
  stop("Column 'gene' not found in Step 25 candidate file.")
}

candidate_genes <- unique(candidates$gene)

cat("Step 25 candidates:", length(candidate_genes), "\n\n")

# ------------------------------------------------------------
# 5. LOAD GENE AND BARCODE INFORMATION
# ------------------------------------------------------------

cat("Loading gene and barcode information...\n")

genes_df <- read.delim(
  genes_file,
  header = FALSE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

barcodes_df <- read.delim(
  barcodes_file,
  header = FALSE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

external_genes <- genes_df[[1]]
external_barcodes <- barcodes_df[[1]]

cat("External genes:", length(external_genes), "\n")
cat("External barcodes:", length(external_barcodes), "\n\n")

# ------------------------------------------------------------
# 6. LOAD METADATA
# ------------------------------------------------------------

cat("Loading metadata...\n")

meta <- read.csv(
  metadata_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (nrow(meta) != length(external_barcodes)) {
  stop(
    "Metadata/barcode mismatch: ",
    nrow(meta),
    " metadata rows vs ",
    length(external_barcodes),
    " barcodes."
  )
}

# First metadata column contains barcode
metadata_barcodes <- meta[[1]]

if (!all(metadata_barcodes == external_barcodes)) {

  cat("Barcode order differs. Matching metadata by barcode...\n")

  idx <- match(external_barcodes, metadata_barcodes)

  if (any(is.na(idx))) {
    stop("Some matrix barcodes could not be matched to metadata.")
  }

  meta <- meta[idx, , drop = FALSE]
}

cat("Metadata cells:", nrow(meta), "\n\n")

# ------------------------------------------------------------
# 7. SELECT CANCER EPITHELIAL CELLS
# ------------------------------------------------------------

if (!"celltype_major" %in% colnames(meta)) {
  stop("Column 'celltype_major' not found in metadata.")
}

cancer_cells <- meta$celltype_major == "Cancer Epithelial"

if (sum(cancer_cells) == 0) {
  stop("No Cancer Epithelial cells found.")
}

cancer_indices <- which(cancer_cells)

cat("Cancer epithelial cells:", length(cancer_indices), "\n")
cat(
  "Cancer epithelial samples:",
  length(unique(meta$orig.ident[cancer_indices])),
  "\n\n"
)

# ------------------------------------------------------------
# 8. MATCH CANDIDATES TO EXTERNAL GENES
# ------------------------------------------------------------

cat("Matching candidates to external gene list...\n")

gene_index <- match(candidate_genes, external_genes)

supported <- !is.na(gene_index)

supported_genes <- candidate_genes[supported]
missing_genes <- candidate_genes[!supported]

cat("Candidates present:", sum(supported), "\n")
cat("Candidates absent:", sum(!supported), "\n\n")

# ------------------------------------------------------------
# 9. READ SPARSE MATRIX
# ------------------------------------------------------------

cat("Reading sparse expression matrix...\n")
cat("This may take some time because the matrix contains\n")
cat("approximately 178 million non-zero entries.\n\n")

full_matrix <- readMM(matrix_file)

cat("Matrix loaded.\n")
cat(
  "Dimensions:",
  nrow(full_matrix),
  "genes x",
  ncol(full_matrix),
  "cells\n\n"
)

# ------------------------------------------------------------
# 10. CONVERT TO SPARSE COMPRESSED FORMAT
# ------------------------------------------------------------

full_matrix <- as(full_matrix, "dgCMatrix")

# ------------------------------------------------------------
# 11. EXTRACT ONLY CANCER EPITHELIAL CELLS
# ------------------------------------------------------------

cat("Restricting matrix to Cancer Epithelial cells...\n")

cancer_matrix <- full_matrix[, cancer_indices, drop = FALSE]

rm(full_matrix)
gc()

cat(
  "Cancer matrix:",
  nrow(cancer_matrix),
  "genes x",
  ncol(cancer_matrix),
  "cells\n\n"
)

# ------------------------------------------------------------
# 12. EXTRACT ONLY CANDIDATE GENES
# ------------------------------------------------------------

candidate_gene_indices <- gene_index[supported]

candidate_matrix <-
  cancer_matrix[candidate_gene_indices, , drop = FALSE]

rm(cancer_matrix)
gc()

rownames(candidate_matrix) <- supported_genes

cat(
  "Candidate matrix:",
  nrow(candidate_matrix),
  "genes x",
  ncol(candidate_matrix),
  "Cancer Epithelial cells\n\n"
)

# ------------------------------------------------------------
# 13. CELL-LEVEL EXPRESSION VALIDATION
# ------------------------------------------------------------

cat("Calculating cell-level expression...\n")

cell_count <- ncol(candidate_matrix)

expression_counts <- Matrix::rowSums(candidate_matrix > 0)

expression_prevalence <-
  as.numeric(expression_counts / cell_count)

mean_expression <-
  as.numeric(Matrix::rowMeans(candidate_matrix))

# ------------------------------------------------------------
# 14. SAMPLE-LEVEL EXPRESSION PREVALENCE
# ------------------------------------------------------------

cat("Calculating sample-level expression consistency...\n")

cancer_meta <- meta[cancer_indices, , drop = FALSE]

samples <- unique(cancer_meta$orig.ident)

sample_results <- vector("list", length(supported_genes))

for (i in seq_along(supported_genes)) {

  gene <- supported_genes[i]

  gene_expression <-
    as.numeric(candidate_matrix[i, ] > 0)

  sample_prev <- numeric(length(samples))

  for (j in seq_along(samples)) {

    sample_name <- samples[j]

    sample_cells <-
      which(cancer_meta$orig.ident == sample_name)

    if (length(sample_cells) == 0) {
      sample_prev[j] <- NA
    } else {

      sample_prev[j] <-
        mean(gene_expression[sample_cells])
    }
  }

  sample_results[[i]] <- data.frame(
    gene = gene,
    samples_tested = length(samples),
    samples_detected_25pct =
      sum(sample_prev >= 0.25, na.rm = TRUE),
    samples_detected_50pct =
      sum(sample_prev >= 0.50, na.rm = TRUE),
    minimum_sample_prevalence =
      min(sample_prev, na.rm = TRUE),
    mean_sample_prevalence =
      mean(sample_prev, na.rm = TRUE),
    maximum_sample_prevalence =
      max(sample_prev, na.rm = TRUE)
  )
}

sample_summary <- do.call(
  rbind,
  sample_results
)

# ------------------------------------------------------------
# 15. COMBINE RESULTS
# ------------------------------------------------------------

expression_results <- data.frame(
  gene = supported_genes,
  external_cells_expressing =
    as.numeric(expression_counts),
  external_expression_prevalence =
    expression_prevalence,
  external_mean_expression =
    mean_expression,
  stringsAsFactors = FALSE
)

results <- merge(
  expression_results,
  sample_summary,
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)

# ------------------------------------------------------------
# 16. ADD STEP 25 INFORMATION
# ------------------------------------------------------------

step25_info <- candidates

keep_columns <- intersect(
  c(
    "gene",
    "final_panel_rank",
    "final_integrated_score",
    "final_review_class",
    "final_biological_category",
    "biological_class",
    "strong_candidate",
    "strong_3sample"
  ),
  colnames(step25_info)
)

step25_info <- step25_info[, keep_columns, drop = FALSE]

results <- merge(
  step25_info,
  results,
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)

# ------------------------------------------------------------
# 17. VALIDATION STATUS
# ------------------------------------------------------------

results$external_validation_status <- "Not_supported"

results$external_validation_status[
  results$external_expression_prevalence >= 0.05
] <- "Externally_expressed"

results$external_validation_status[
  results$external_expression_prevalence >= 0.25 &
  results$samples_detected_25pct >= 3
] <- "Strong_external_expression_support"

# ------------------------------------------------------------
# 18. SORT RESULTS
# ------------------------------------------------------------

results <- results[
  order(
    -results$external_expression_prevalence,
    -results$mean_sample_prevalence
  ),
]

# ------------------------------------------------------------
# 19. SAVE SUPPORTED CANDIDATES
# ------------------------------------------------------------

supported_file <- file.path(
  output_dir,
  "27_external_expression_validation.csv"
)

write.csv(
  results,
  supported_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 20. SAMPLE CONSISTENCY FILE
# ------------------------------------------------------------

sample_file <- file.path(
  output_dir,
  "27_external_sample_consistency.csv"
)

write.csv(
  sample_summary[
    order(
      -sample_summary$mean_sample_prevalence
    ),
  ],
  sample_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 21. TOP CANDIDATES
# ------------------------------------------------------------

top_candidates <- head(results, 20)

top_file <- file.path(
  output_dir,
  "27_external_top_candidates.csv"
)

write.csv(
  top_candidates,
  top_file,
  row.names = FALSE
)

# ------------------------------------------------------------
# 22. VALIDATED GENE LIST
# ------------------------------------------------------------

validated_genes <-
  results$gene[
    results$external_expression_prevalence >= 0.25 &
    results$samples_detected_25pct >= 3
  ]

validated_file <- file.path(
  output_dir,
  "27_external_validated_gene_list.txt"
)

writeLines(
  validated_genes,
  validated_file
)

# ------------------------------------------------------------
# 23. MISSING CANDIDATES
# ------------------------------------------------------------

missing_file <- file.path(
  output_dir,
  "27_candidates_not_found_externally.txt"
)

writeLines(
  missing_genes,
  missing_file
)

# ------------------------------------------------------------
# 24. SUMMARY
# ------------------------------------------------------------

strong_support <-
  sum(
    results$external_validation_status ==
      "Strong_external_expression_support"
  )

expressed_support <-
  sum(
    results$external_validation_status ==
      "Externally_expressed"
  )

summary_file <- file.path(
  output_dir,
  "27_external_expression_validation_summary.txt"
)

summary_text <- c(
  "====================================================",
  "STEP 27 EXTERNAL EXPRESSION VALIDATION SUMMARY",
  "====================================================",
  "",
  paste("Step 25 candidates:", length(candidate_genes)),
  paste("Candidates present in GSE176078:", length(supported_genes)),
  paste("Candidates absent:", length(missing_genes)),
  paste("Cancer epithelial cells:", length(cancer_indices)),
  paste(
    "External cancer epithelial samples:",
    length(samples)
  ),
  "",
  paste(
    "Externally expressed candidates:",
    expressed_support
  ),
  paste(
    "Strong external expression support:",
    strong_support
  ),
  "",
  "Strong external expression support criteria:",
  "Expression prevalence >= 25% in cancer epithelial cells",
  "AND expression prevalence >= 25% in at least 3 samples",
  "",
  "Interpretation:",
  "External expression support indicates reproducible",
  "expression in an independent breast-cancer dataset.",
  "It does not establish clinical biomarker validity.",
  "It does not establish causality or diagnostic utility.",
  "",
  "The external dataset is used as an independent",
  "biological/expression validation resource.",
  "",
  "====================================================",
  "STEP 27 COMPLETED",
  "===================================================="
)

writeLines(
  summary_text,
  summary_file
)

# ------------------------------------------------------------
# 25. PRINT SUMMARY
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("STEP 27 COMPLETED\n")
cat("====================================================\n\n")

cat(
  "Step 25 candidates:",
  length(candidate_genes),
  "\n"
)

cat(
  "Candidates present in GSE176078:",
  length(supported_genes),
  "\n"
)

cat(
  "Candidates absent:",
  length(missing_genes),
  "\n"
)

cat(
  "Cancer epithelial cells:",
  length(cancer_indices),
  "\n"
)

cat(
  "External samples:",
  length(samples),
  "\n\n"
)

cat(
  "Externally expressed candidates:",
  expressed_support,
  "\n"
)

cat(
  "Strong external expression support:",
  strong_support,
  "\n\n"
)

cat("Top candidates:\n")
print(
  top_candidates[
    ,
    intersect(
      c(
        "gene",
        "external_expression_prevalence",
        "external_mean_expression",
        "mean_sample_prevalence",
        "samples_detected_25pct",
        "external_validation_status"
      ),
      colnames(top_candidates)
    )
  ],
  row.names = FALSE
)

cat("\nResults saved to:\n")
cat(output_dir, "\n\n")
cat("Done.\n")