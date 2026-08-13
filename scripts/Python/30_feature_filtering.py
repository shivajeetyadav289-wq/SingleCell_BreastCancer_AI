# ============================================================
# STEP 30
# FEATURE FILTERING AND FEATURE-SET CONSTRUCTION
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

OUTPUT_DIR = "results/ai_ml/30_feature_selection"

os.makedirs(OUTPUT_DIR, exist_ok=True)


print("=" * 60)
print("STEP 30 - FEATURE FILTERING")
print("=" * 60)


# ------------------------------------------------------------
# Load expression matrix
# ------------------------------------------------------------

matrix_file = os.path.join(
    DATA_DIR,
    "28_expression_matrix.mtx.gz"
)

print("\nLoading expression matrix...")

X = mmread(matrix_file).tocsr()

print(
    "Matrix:",
    X.shape[0],
    "genes x",
    X.shape[1],
    "cells"
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

gene_names = genes["gene"].astype(str).values


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

print(
    "Metadata:",
    meta.shape[0],
    "cells"
)


# ------------------------------------------------------------
# Verify dimensions
# ------------------------------------------------------------

assert X.shape[0] == len(gene_names)
assert X.shape[1] == len(meta)

print("\nDimension checks: PASSED")


# ------------------------------------------------------------
# Gene detection statistics
# ------------------------------------------------------------

print("\nCalculating gene detection statistics...")

# Number of cells with expression > 0
detected_cells = np.asarray(
    (X > 0).sum(axis=1)
).ravel()

detection_fraction = (
    detected_cells / X.shape[1]
)

# Total expression
total_expression = np.asarray(
    X.sum(axis=1)
).ravel()

# Mean expression
mean_expression = (
    total_expression / X.shape[1]
)


gene_stats = pd.DataFrame({

    "gene": gene_names,

    "detected_cells":
        detected_cells,

    "detection_fraction":
        detection_fraction,

    "total_expression":
        total_expression,

    "mean_expression":
        mean_expression
})


# ------------------------------------------------------------
# Save complete gene statistics
# ------------------------------------------------------------

gene_stats = gene_stats.sort_values(
    "detection_fraction",
    ascending=False
)

gene_stats.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "30_gene_expression_statistics.csv"
    ),
    index=False
)


# ------------------------------------------------------------
# Basic expression filtering
# ------------------------------------------------------------

MIN_DETECTION_FRACTION = 0.01

MIN_DETECTED_CELLS = int(
    np.ceil(
        X.shape[1] * MIN_DETECTION_FRACTION
    )
)

print("\n===== FILTERING =====")

print(
    "Detection threshold:",
    f"{MIN_DETECTION_FRACTION:.1%}"
)

print(
    "Minimum detected cells:",
    MIN_DETECTED_CELLS
)


keep = (
    detected_cells >= MIN_DETECTED_CELLS
)

print(
    "Genes before filtering:",
    len(gene_names)
)

print(
    "Genes retained:",
    int(keep.sum())
)

print(
    "Genes removed:",
    int((~keep).sum())
)


# ------------------------------------------------------------
# Filter expression matrix
# ------------------------------------------------------------

X_filtered = X[keep, :].tocsr()

filtered_genes = gene_names[keep]


print(
    "\nFiltered matrix:",
    X_filtered.shape[0],
    "genes x",
    X_filtered.shape[1],
    "cells"
)


# ------------------------------------------------------------
# Save filtered sparse matrix
# ------------------------------------------------------------

filtered_matrix_file = os.path.join(
    OUTPUT_DIR,
    "30_filtered_expression_matrix.mtx"
)

print("\nSaving filtered matrix...")

from scipy.io import mmwrite

mmwrite(
    filtered_matrix_file,
    X_filtered
)


# Compress using system gzip
os.system(
    f"gzip -f '{filtered_matrix_file}'"
)


# ------------------------------------------------------------
# Save filtered genes
# ------------------------------------------------------------

pd.DataFrame({
    "gene": filtered_genes
}).to_csv(
    os.path.join(
        OUTPUT_DIR,
        "30_filtered_genes.tsv"
    ),
    sep="\t",
    index=False,
    header=False
)


# ------------------------------------------------------------
# Create feature summary
# ------------------------------------------------------------

feature_summary = pd.DataFrame({

    "feature_set": [
        "all_genes",
        "expression_filtered_genes"
    ],

    "gene_count": [
        len(gene_names),
        len(filtered_genes)
    ],

    "cells": [
        X.shape[1],
        X_filtered.shape[1]
    ],

    "detection_threshold": [
        0,
        MIN_DETECTION_FRACTION
    ]
})


feature_summary.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "30_feature_set_summary.csv"
    ),
    index=False
)


# ------------------------------------------------------------
# Save metadata unchanged
# ------------------------------------------------------------

meta.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "30_ml_metadata.csv"
    )
)


# ------------------------------------------------------------
# Final report
# ------------------------------------------------------------

print("\n" + "=" * 60)
print("STEP 30 COMPLETED")
print("=" * 60)

print(
    "\nOriginal genes:",
    len(gene_names)
)

print(
    "Expression-filtered genes:",
    len(filtered_genes)
)

print(
    "Cells:",
    X.shape[1]
)

print(
    "\nResults saved to:"
)

print(OUTPUT_DIR)