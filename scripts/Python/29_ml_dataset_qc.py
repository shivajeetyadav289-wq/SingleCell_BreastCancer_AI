# ============================================================
# STEP 29
# ML DATASET QUALITY CONTROL
# Single-Cell Breast Cancer AI/ML Project
# ============================================================

import os
import numpy as np
import pandas as pd

from scipy.io import mmread
from scipy import sparse


# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

DATA_DIR = "results/ai_ml/28_ml_dataset"
OUTPUT_DIR = "results/ai_ml/29_dataset_qc"

os.makedirs(OUTPUT_DIR, exist_ok=True)


print("=" * 60)
print("STEP 29 - ML DATASET QUALITY CONTROL")
print("=" * 60)


# ------------------------------------------------------------
# Load expression matrix
# ------------------------------------------------------------

matrix_file = os.path.join(
    DATA_DIR,
    "28_expression_matrix.mtx.gz"
)

print("\nLoading expression matrix...")

X = mmread(matrix_file)

print("Original matrix type:", type(X))
print("Original dimensions:", X.shape)


# Convert to CSR sparse matrix

X = X.tocsr()

print("Sparse format:", type(X))
print("Genes:", X.shape[0])
print("Cells:", X.shape[1])

print("Non-zero entries:", X.nnz)

total_elements = X.shape[0] * X.shape[1]

sparsity = 1 - (X.nnz / total_elements)

print(
    "Sparsity:",
    f"{sparsity:.4%}"
)


# ------------------------------------------------------------
# Load genes
# ------------------------------------------------------------

genes_file = os.path.join(
    DATA_DIR,
    "28_genes.tsv"
)

genes = pd.read_csv(
    genes_file,
    sep="\t",
    header=None,
    names=["gene"]
)

print("\n===== GENE CHECK =====")

print("Number of genes:", len(genes))

print(
    "Unique genes:",
    genes["gene"].nunique()
)

print("\nFirst 10 genes:")

print(
    genes.head(10).to_string(index=False)
)


# ------------------------------------------------------------
# Load cells
# ------------------------------------------------------------

cells_file = os.path.join(
    DATA_DIR,
    "28_cells.tsv"
)

cells = pd.read_csv(
    cells_file,
    sep="\t",
    header=None,
    names=["cell"]
)

print("\n===== CELL CHECK =====")

print("Number of cells:", len(cells))

print(
    "Unique cells:",
    cells["cell"].nunique()
)

print("\nFirst 10 cells:")

print(
    cells.head(10).to_string(index=False)
)


# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

metadata_file = os.path.join(
    DATA_DIR,
    "28_ml_metadata.csv"
)

meta = pd.read_csv(
    metadata_file,
    index_col=0
)

print("\n===== METADATA CHECK =====")

print("Metadata rows:", len(meta))
print("Metadata columns:", len(meta.columns))

print("\nMetadata columns:")

for col in meta.columns:
    print("-", col)


# ------------------------------------------------------------
# Dimension checks
# ------------------------------------------------------------

print("\n===== DIMENSION CHECK =====")

assert X.shape[0] == len(genes), (
    "Gene count does not match expression matrix."
)

assert X.shape[1] == len(cells), (
    "Cell count does not match expression matrix."
)

assert X.shape[1] == len(meta), (
    "Cell count does not match metadata."
)

print("✓ Gene dimensions match")
print("✓ Cell dimensions match")
print("✓ Metadata dimensions match")


# ------------------------------------------------------------
# Cell ID consistency
# ------------------------------------------------------------

print("\n===== CELL ID CHECK =====")

matrix_cells = cells["cell"].astype(str).values
metadata_cells = meta.index.astype(str).values

if np.array_equal(matrix_cells, metadata_cells):

    print("✓ Expression cell order matches metadata")

else:

    print(
        "WARNING: Expression cell order does not "
        "match metadata."
    )

    common = np.intersect1d(
        matrix_cells,
        metadata_cells
    )

    print(
        "Common cells:",
        len(common)
    )


# ------------------------------------------------------------
# Label distribution
# ------------------------------------------------------------

print("\n===== ML LABEL DISTRIBUTION =====")

if "ML_label" in meta.columns:

    print(
        meta["ML_label"]
        .value_counts(dropna=False)
        .sort_index()
    )

if "ML_class" in meta.columns:

    print(
        meta["ML_class"]
        .value_counts(dropna=False)
    )


# ------------------------------------------------------------
# Sample distribution
# ------------------------------------------------------------

print("\n===== SAMPLE DISTRIBUTION =====")

if "patient_id" in meta.columns:

    sample_table = pd.crosstab(
        meta["patient_id"],
        meta["ML_class"]
    )

    print(sample_table)


# ------------------------------------------------------------
# Missing values
# ------------------------------------------------------------

print("\n===== MISSING VALUES =====")

missing = meta.isna().sum()

missing = missing[
    missing > 0
]

if len(missing) == 0:

    print("✓ No missing metadata values")

else:

    print(missing)


# ------------------------------------------------------------
# Expression statistics
# ------------------------------------------------------------

print("\n===== EXPRESSION MATRIX =====")

print(
    "Minimum:",
    X.min()
)

print(
    "Maximum:",
    X.max()
)

print(
    "Non-zero values:",
    X.nnz
)

memory_mb = (
    X.data.nbytes +
    X.indices.nbytes +
    X.indptr.nbytes
) / (1024 ** 2)

print(
    "Matrix memory estimate:",
    f"{memory_mb:.2f} MB"
)


# ------------------------------------------------------------
# Cells with expression
# ------------------------------------------------------------

cell_nonzero = np.diff(X.indptr)

print("\n===== CELL EXPRESSION =====")

print(
    "Minimum detected genes:",
    cell_nonzero.min()
)

print(
    "Maximum detected genes:",
    cell_nonzero.max()
)

print(
    "Median detected genes:",
    np.median(cell_nonzero)
)


# ------------------------------------------------------------
# Genes with expression
# ------------------------------------------------------------

gene_nonzero = np.diff(X.tocsc().indptr)

print("\n===== GENE EXPRESSION =====")

print(
    "Genes detected in at least one cell:",
    np.sum(gene_nonzero > 0)
)

print(
    "Genes detected in >=10 cells:",
    np.sum(gene_nonzero >= 10)
)

print(
    "Genes detected in >=1% cells:",
    np.sum(gene_nonzero >= X.shape[1] * 0.01)
)


# ------------------------------------------------------------
# Save QC summary
# ------------------------------------------------------------

summary = {

    "genes": X.shape[0],

    "cells": X.shape[1],

    "nonzero_entries": X.nnz,

    "sparsity": sparsity,

    "malignant_cells":
        int((meta["ML_label"] == 1).sum()),

    "diploid_non_malignant_cells":
        int((meta["ML_label"] == 0).sum()),

    "samples":
        meta["patient_id"].nunique(),

    "metadata_missing_values":
        int(meta.isna().sum().sum())

}

summary_df = pd.DataFrame(
    [summary]
)

summary_file = os.path.join(
    OUTPUT_DIR,
    "29_ml_dataset_qc_summary.csv"
)

summary_df.to_csv(
    summary_file,
    index=False
)


# ------------------------------------------------------------
# Save sample distribution
# ------------------------------------------------------------

if "patient_id" in meta.columns:

    sample_table.to_csv(
        os.path.join(
            OUTPUT_DIR,
            "29_sample_class_distribution.csv"
        )
    )


# ------------------------------------------------------------
# Final
# ------------------------------------------------------------

print("\n" + "=" * 60)
print("STEP 29 COMPLETED")
print("=" * 60)

print("\nQC summary saved to:")

print(summary_file)