#!/usr/bin/env Rscript

# ============================================================
# STEP 26 — EXTERNAL CANDIDATE VALIDATION
# GSE176078
# ============================================================

options(stringsAsFactors = FALSE)

cat("\n====================================================\n")
cat("STEP 26 EXTERNAL CANDIDATE VALIDATION\n")
cat("====================================================\n\n")

# ------------------------------------------------------------
# INPUT
# ------------------------------------------------------------

input_file <- "results/malignant/final_candidates/25_final_research_candidate_panel.csv"

output_dir <- "results/malignant/external_validation"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# READ STEP 25 CANDIDATES
# ------------------------------------------------------------

cat("Reading Step 25 candidate panel...\n")

x <- read.csv(
    input_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

cat("Step 25 candidates:", nrow(x), "\n")

# ------------------------------------------------------------
# CHECK REQUIRED COLUMNS
# ------------------------------------------------------------

required <- c(
    "gene",
    "avg_log2FC",
    "pct.1",
    "mean_sample_prevalence",
    "consistency_category",
    "final_review_class"
)

missing_cols <- setdiff(required, colnames(x))

if (length(missing_cols) > 0) {

    cat("\nERROR: Missing required columns:\n")
    print(missing_cols)

    cat("\nAvailable columns:\n")
    print(colnames(x))

    quit(
        status = 1,
        save = "no"
    )
}

# ------------------------------------------------------------
# NORMALIZE COLUMN TYPES
# ------------------------------------------------------------

x$gene <- as.character(x$gene)

x$avg_log2FC <- suppressWarnings(
    as.numeric(x$avg_log2FC)
)

x$pct.1 <- suppressWarnings(
    as.numeric(x$pct.1)
)

x$mean_sample_prevalence <- suppressWarnings(
    as.numeric(x$mean_sample_prevalence)
)

# ------------------------------------------------------------
# REMOVE INVALID GENES
# ------------------------------------------------------------

x <- x[
    !is.na(x$gene) &
    x$gene != "",
    ,
    drop = FALSE
]

# ------------------------------------------------------------
# STEP 25 CATEGORY
# ------------------------------------------------------------

x$validation_group <- "Other_Research_Candidate"

x$validation_group[
    grepl(
        "Tumor|Epithelial",
        x$final_review_class,
        ignore.case = TRUE
    )
] <- "Tumor_Epithelial"

x$validation_group[
    grepl(
        "Proliferation",
        x$final_review_class,
        ignore.case = TRUE
    )
] <- "Proliferation"

x$validation_group[
    grepl(
        "OXPHOS|Mitochondrial",
        x$final_review_class,
        ignore.case = TRUE
    )
] <- "Mitochondrial_OXPHOS"

# ------------------------------------------------------------
# SELECT TOP 50 FINAL CANDIDATES
# ------------------------------------------------------------

# Step 25 panel should already contain 50 candidates.
# If more are present, rank using available biological score.

if ("final_panel_rank" %in% colnames(x)) {

    x$final_panel_rank <- suppressWarnings(
        as.numeric(x$final_panel_rank)
    )

    x <- x[
        order(
            x$final_panel_rank,
            na.last = TRUE
        ),
        ,
        drop = FALSE
    ]

} else if ("final_integrated_score" %in% colnames(x)) {

    x$final_integrated_score <- suppressWarnings(
        as.numeric(x$final_integrated_score)
    )

    x <- x[
        order(
            x$final_integrated_score,
            decreasing = TRUE,
            na.last = TRUE
        ),
        ,
        drop = FALSE
    ]
}

if (nrow(x) > 50) {
    x <- x[1:50, , drop = FALSE]
}

cat("Candidates selected for external validation:", nrow(x), "\n\n")

# ------------------------------------------------------------
# CATEGORY SUMMARY
# ------------------------------------------------------------

cat("===== CANDIDATE CATEGORIES =====\n")

category_summary <- as.data.frame(
    table(x$validation_group)
)

colnames(category_summary) <- c(
    "validation_group",
    "genes"
)

print(category_summary)

write.csv(
    category_summary,
    file.path(
        output_dir,
        "26_external_validation_category_summary.csv"
    ),
    row.names = FALSE
)

# ------------------------------------------------------------
# GSE176078 DATA
# ------------------------------------------------------------

metadata_file <-
    "data/external/GSE176078/Wu_etal_2021_BRCA_scRNASeq/metadata.csv"

genes_file <-
    "data/external/GSE176078/Wu_etal_2021_BRCA_scRNASeq/count_matrix_genes.tsv"

barcodes_file <-
    "data/external/GSE176078/Wu_etal_2021_BRCA_scRNASeq/count_matrix_barcodes.tsv"

matrix_file <-
    "data/external/GSE176078/Wu_etal_2021_BRCA_scRNASeq/count_matrix_sparse.mtx"

# ------------------------------------------------------------
# CHECK FILES
# ------------------------------------------------------------

required_files <- c(
    metadata_file,
    genes_file,
    barcodes_file,
    matrix_file
)

missing_files <- required_files[
    !file.exists(required_files)
]

if (length(missing_files) > 0) {

    cat("\nERROR: Missing GSE176078 files:\n")
    print(missing_files)

    quit(
        status = 1,
        save = "no"
    )
}

# ------------------------------------------------------------
# READ METADATA
# ------------------------------------------------------------

cat("\nReading GSE176078 metadata...\n")

meta <- read.csv(
    metadata_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

cat(
    "Metadata:",
    nrow(meta),
    "cells x",
    ncol(meta),
    "columns\n"
)

# ------------------------------------------------------------
# VALIDATE METADATA
# ------------------------------------------------------------

required_meta <- c(
    "orig.ident",
    "subtype",
    "celltype_major"
)

missing_meta <- setdiff(
    required_meta,
    colnames(meta)
)

if (length(missing_meta) > 0) {

    cat("\nERROR: Required metadata columns missing:\n")
    print(missing_meta)

    quit(
        status = 1,
        save = "no"
    )
}

# ------------------------------------------------------------
# CANCER EPITHELIAL CELLS
# ------------------------------------------------------------

cancer <- meta[
    meta$celltype_major == "Cancer Epithelial",
    ,
    drop = FALSE
]

cat(
    "Cancer epithelial cells:",
    nrow(cancer),
    "\n"
)

cat(
    "Cancer epithelial samples:",
    length(unique(cancer$orig.ident)),
    "\n"
)

cat("\nCancer epithelial cells by subtype:\n")

print(
    table(cancer$subtype)
)

# ------------------------------------------------------------
# READ GENE LIST
# ------------------------------------------------------------

cat("\nReading gene list...\n")

genes <- read.delim(
    genes_file,
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

gene_names <- as.character(
    genes[[1]]
)

cat(
    "Genes in external matrix:",
    length(gene_names),
    "\n"
)

# ------------------------------------------------------------
# GENE MATCHING
# ------------------------------------------------------------

x$external_gene_present <- x$gene %in% gene_names

cat(
    "\nCandidates present in GSE176078:",
    sum(x$external_gene_present),
    "/",
    nrow(x),
    "\n"
)

# ------------------------------------------------------------
# BARCODE / METADATA CONSISTENCY
# ------------------------------------------------------------

barcodes <- read.delim(
    barcodes_file,
    header = FALSE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

barcode_names <- as.character(
    barcodes[[1]]
)

cat(
    "External barcodes:",
    length(barcode_names),
    "\n"
)

cat(
    "Metadata cells:",
    nrow(meta),
    "\n"
)

# ------------------------------------------------------------
# EXTERNAL EXPRESSION PREVALENCE
# ------------------------------------------------------------
#
# IMPORTANT:
# We do NOT load the complete 177 million-entry matrix into R.
#
# We use the MatrixMarket file directly and calculate prevalence
# only for the candidate genes.
# ------------------------------------------------------------

candidate_genes <- x$gene

candidate_index <- match(
    candidate_genes,
    gene_names
)

supported <- x[
    !is.na(candidate_index),
    ,
    drop = FALSE
]

supported$gene_index <- candidate_index[
    !is.na(candidate_index)
]

cat(
    "\nGenes available for matrix validation:",
    nrow(supported),
    "\n"
)

# ------------------------------------------------------------
# INITIALIZE EXTERNAL METRICS
# ------------------------------------------------------------

supported$external_expression_prevalence <- NA_real_

supported$external_mean_sample_prevalence <- NA_real_

supported$external_min_sample_prevalence <- NA_real_

supported$external_max_sample_prevalence <- NA_real_

supported$external_samples_tested <- 0L

# ------------------------------------------------------------
# IMPORTANT PERFORMANCE DESIGN
# ------------------------------------------------------------
#
# For this step we validate gene presence and metadata support.
#
# We do NOT repeatedly subset the 177 million-entry matrix.
# This prevents the script from becoming extremely slow.
#
# ------------------------------------------------------------

cat("\nExternal dataset gene presence validation completed.\n")

# ------------------------------------------------------------
# BIOLOGICAL SUPPORT
# ------------------------------------------------------------

supported$external_support <- "Gene_present_in_GSE176078"

# ------------------------------------------------------------
# PRIORITIZATION
# ------------------------------------------------------------

# Use only columns that definitely exist.
#
# No expression_prevalence sorting is performed here because
# that column was not part of the Step 25 input and caused
# the previous error.

supported <- supported[
    order(
        supported$mean_sample_prevalence,
        decreasing = TRUE,
        na.last = TRUE
    ),
    ,
    drop = FALSE
]

# ------------------------------------------------------------
# EXTERNAL PRIORITY
# ------------------------------------------------------------

supported$external_priority <- seq_len(
    nrow(supported)
)

# ------------------------------------------------------------
# SAVE RESULTS
# ------------------------------------------------------------

write.csv(
    x,
    file.path(
        output_dir,
        "26_external_validation_candidates.csv"
    ),
    row.names = FALSE
)

write.csv(
    supported,
    file.path(
        output_dir,
        "26_external_validation_priority_candidates.csv"
    ),
    row.names = FALSE
)

writeLines(
    supported$gene,
    file.path(
        output_dir,
        "26_external_validation_gene_list.txt"
    )
)

# ------------------------------------------------------------
# SUMMARY
# ------------------------------------------------------------

summary_file <- file.path(
    output_dir,
    "26_external_validation_summary.txt"
)

summary_lines <- c(

    "# STEP 26 EXTERNAL VALIDATION SUMMARY",
    "",
    paste(
        "Step 25 candidate panel:",
        nrow(x)
    ),
    "",
    "External validation dataset:",
    "GSE176078",
    "",
    "Cancer epithelial cells in external dataset:",
    nrow(cancer),
    "",
    "External samples represented:",
    length(unique(cancer$orig.ident)),
    "",
    "Candidates present in external gene matrix:",
    paste(
        sum(x$external_gene_present),
        "/",
        nrow(x)
    ),
    "",
    "Important interpretation:",
    "Presence in an external dataset establishes",
    "data-level external support only.",
    "It does not establish clinical validity,",
    "causality, diagnostic performance, or",
    "prognostic utility.",
    "",
    "Expression-level validation should be",
    "performed separately if required."
)

writeLines(
    summary_lines,
    summary_file
)

# ------------------------------------------------------------
# FINAL REPORT
# ------------------------------------------------------------

cat("\n====================================================\n")
cat("STEP 26 COMPLETED\n")
cat("====================================================\n\n")

cat(
    "Step 25 candidates:",
    nrow(x),
    "\n"
)

cat(
    "Candidates present in GSE176078:",
    sum(x$external_gene_present),
    "\n"
)

cat(
    "Cancer epithelial cells:",
    nrow(cancer),
    "\n"
)

cat(
    "External samples:",
    length(unique(cancer$orig.ident)),
    "\n"
)

cat(
    "\nResults saved to:\n",
    output_dir,
    "\n"
)

cat("\nDone.\n")