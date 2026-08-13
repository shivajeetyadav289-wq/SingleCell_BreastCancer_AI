# ============================================================
# STEP 28
# EXTRACT HIGH-CONFIDENCE ML DATASET
# GSE228499 Single-Cell Breast Cancer AI/ML Project
# ============================================================

suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
})

cat("====================================================\n")
cat("STEP 28 - ML DATASET EXTRACTION\n")
cat("====================================================\n\n")

input_file <- "data/processed/GSE228499_final_malignant_annotated.rds"

output_dir <- "results/ai_ml/28_ml_dataset"

if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

# ------------------------------------------------------------
# Load object
# ------------------------------------------------------------

cat("Loading Seurat object...\n")

obj <- readRDS(input_file)

cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")
cat("Genes:", nrow(obj), "\n")
cat("Cells:", ncol(obj), "\n\n")

meta <- obj@meta.data

# ------------------------------------------------------------
# Define ML labels
# ------------------------------------------------------------

positive <- meta$malignant_status ==
    "CNV_supported_malignant_candidate"

negative <- meta$cnv_supported_status ==
    "CopyKAT_diploid" &
    meta$malignant_status ==
    "Non_malignant_or_unresolved"

keep <- positive | negative

cat("===== LABEL SUMMARY =====\n")

cat(
    "Positive malignant cells:",
    sum(positive, na.rm = TRUE),
    "\n"
)

cat(
    "Negative high-confidence diploid cells:",
    sum(negative, na.rm = TRUE),
    "\n"
)

cat(
    "Total ML cells:",
    sum(keep, na.rm = TRUE),
    "\n\n"
)

# ------------------------------------------------------------
# Check sample distribution
# ------------------------------------------------------------

ml_meta <- meta[keep, , drop = FALSE]

ml_meta$ML_label <- ifelse(
    positive[keep],
    1,
    0
)

ml_meta$ML_class <- ifelse(
    ml_meta$ML_label == 1,
    "malignant",
    "diploid_non_malignant"
)

cat("===== ML LABEL DISTRIBUTION =====\n")
print(table(ml_meta$ML_class))

cat("\n===== SAMPLE DISTRIBUTION =====\n")
print(
    table(
        ml_meta$patient_id,
        ml_meta$ML_class
    )
)

# ------------------------------------------------------------
# Extract cell IDs
# ------------------------------------------------------------

cell_ids <- rownames(ml_meta)

cat("\nSelected cells:", length(cell_ids), "\n")

# ------------------------------------------------------------
# Extract normalized expression
# ------------------------------------------------------------

cat("\n===== EXPRESSION EXTRACTION =====\n")

cat("Default assay:", DefaultAssay(obj), "\n")

rna <- obj[["RNA"]]

# Seurat v5 has sample-specific data layers.
# We extract each sample separately and combine the
# selected cells into one sparse matrix.

data_layers <- paste0(
    "data.",
    unique(ml_meta$patient_id)
)

cat("Expected data layers:\n")
print(data_layers)

available_layers <- Layers(rna)

cat("\nAvailable layers:\n")
print(available_layers)

missing_layers <- setdiff(
    data_layers,
    available_layers
)

if (length(missing_layers) > 0) {
    stop(
        "Missing expected data layers: ",
        paste(missing_layers, collapse = ", ")
    )
}

# ------------------------------------------------------------
# Extract each sample
# ------------------------------------------------------------

mat_list <- list()

for (sample in unique(ml_meta$patient_id)) {

    cat("\nProcessing sample:", sample, "\n")

    layer_name <- paste0("data.", sample)

    sample_cells <- cell_ids[
        ml_meta$patient_id == sample
    ]

    cat(
        "Selected cells:",
        length(sample_cells),
        "\n"
    )

    layer_mat <- LayerData(
        rna,
        layer = layer_name
    )

    # Keep only selected cells belonging to this sample
    common_cells <- intersect(
        colnames(layer_mat),
        sample_cells
    )

    if (length(common_cells) == 0) {
        stop(
            "No matching cells found for sample: ",
            sample
        )
    }

    layer_mat <- layer_mat[
        ,
        common_cells,
        drop = FALSE
    ]

    mat_list[[sample]] <- layer_mat

    rm(layer_mat)

    gc()
}

# ------------------------------------------------------------
# Combine sample matrices
# ------------------------------------------------------------

cat("\n===== COMBINING MATRICES =====\n")

expr <- do.call(
    cbind,
    mat_list
)

cat(
    "Expression dimensions:",
    nrow(expr),
    "genes x",
    ncol(expr),
    "cells\n"
)

# ------------------------------------------------------------
# Ensure cell order matches metadata
# ------------------------------------------------------------

expr <- expr[
    ,
    cell_ids,
    drop = FALSE
]

if (!identical(
    colnames(expr),
    rownames(ml_meta)
)) {

    stop(
        "Expression matrix and metadata cell order do not match."
    )
}

cat("Cell order verified.\n")

# ------------------------------------------------------------
# Save sparse matrix
# ------------------------------------------------------------

cat("\n===== SAVING DATASET =====\n")

writeMM(
    expr,
    file.path(
        output_dir,
        "28_expression_matrix.mtx"
    )
)

system2(
    "gzip",
    c(
        "-f",
        file.path(
            output_dir,
            "28_expression_matrix.mtx"
        )
    )
)

# Genes
write.table(
    data.frame(
        gene = rownames(expr)
    ),
    file.path(
        output_dir,
        "28_genes.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
)

# Cells
write.table(
    data.frame(
        cell = colnames(expr)
    ),
    file.path(
        output_dir,
        "28_cells.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
)

# Metadata
write.csv(
    ml_meta,
    file.path(
        output_dir,
        "28_ml_metadata.csv"
    ),
    row.names = TRUE
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

summary <- data.frame(
    metric = c(
        "genes",
        "cells",
        "malignant_cells",
        "diploid_non_malignant_cells",
        "samples"
    ),
    value = c(
        nrow(expr),
        ncol(expr),
        sum(ml_meta$ML_label == 1),
        sum(ml_meta$ML_label == 0),
        length(unique(ml_meta$patient_id))
    )
)

write.csv(
    summary,
    file.path(
        output_dir,
        "28_ml_dataset_summary.csv"
    ),
    row.names = FALSE
)

cat("\n====================================================\n")
cat("STEP 28 COMPLETED\n")
cat("====================================================\n")

print(summary)

cat("\nOutput directory:\n")
cat(output_dir, "\n")