# ============================================================
# STEP 43
# CELL-TYPE / TUMOR-SPECIFICITY VALIDATION
# ============================================================
#
# Purpose:
#
# Determine whether the final six consensus genes are:
#
#   1. Malignant-associated
#   2. Tumor-epithelial associated
#   3. Broadly expressed
#   4. Weakly specific
#
# Data source:
# GSE228499 discovery dataset
#
# No new ML model is trained.
# No external dataset is downloaded.
#
# ============================================================

import os
import gzip
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# 1. INPUT FILES
# ============================================================

EXPRESSION_FILE = (
    "results/ai_ml/28_ml_dataset/"
    "28_expression_matrix.mtx.gz"
)

GENES_FILE = (
    "results/ai_ml/28_ml_dataset/"
    "28_genes.tsv"
)

CELLS_FILE = (
    "results/ai_ml/28_ml_dataset/"
    "28_cells.tsv"
)

METADATA_FILE = (
    "results/ai_ml/28_ml_dataset/"
    "28_ml_metadata.csv"
)

STEP42_FILE = (
    "results/ai_ml/42_external_validation/"
    "42_six_gene_external_validation.csv"
)

STEP39_FILE = (
    "results/ai_ml/39_consensus_visualization/"
    "39_final_candidate_panel.csv"
)


# ============================================================
# 2. OUTPUT DIRECTORY
# ============================================================

OUTPUT_DIR = (
    "results/ai_ml/43_celltype_tumor_specificity"
)

FIGURE_DIR = os.path.join(
    OUTPUT_DIR,
    "figures"
)

os.makedirs(
    OUTPUT_DIR,
    exist_ok=True
)

os.makedirs(
    FIGURE_DIR,
    exist_ok=True
)


# ============================================================
# 3. FINAL SIX
# ============================================================

FINAL_SIX = [
    "COX8A",
    "MRPS7",
    "PRSS23",
    "ETFB",
    "FCGRT",
    "EPN3"
]


# ============================================================
# 4. START
# ============================================================

print("=" * 70)
print("STEP 43 - CELL-TYPE / TUMOR-SPECIFICITY VALIDATION")
print("=" * 70)


# ============================================================
# 5. FILE CHECK
# ============================================================

print("\n===== INPUT FILE CHECK =====")

required_files = {
    "Expression matrix": EXPRESSION_FILE,
    "Genes": GENES_FILE,
    "Cells": CELLS_FILE,
    "Metadata": METADATA_FILE,
    "Step 42": STEP42_FILE,
    "Step 39": STEP39_FILE
}

for name, path in required_files.items():

    if os.path.exists(path):
        print(f"{name}: FOUND")

    else:
        raise FileNotFoundError(
            f"\nMissing file:\n{name}\n{path}"
        )


# ============================================================
# 6. LOAD GENES
# ============================================================

print("\n===== LOADING GENES =====")

genes = pd.read_csv(
    GENES_FILE,
    sep="\t",
    header=None
)

gene_names = (
    genes.iloc[:, 0]
    .astype(str)
    .str.strip()
    .tolist()
)

print(
    "Genes:",
    len(gene_names)
)


# ============================================================
# 7. LOAD CELLS
# ============================================================

print("\n===== LOADING CELLS =====")

cells = pd.read_csv(
    CELLS_FILE,
    sep="\t",
    header=None
)

cell_names = (
    cells.iloc[:, 0]
    .astype(str)
    .str.strip()
    .tolist()
)

print(
    "Cells:",
    len(cell_names)
)


# ============================================================
# 8. LOAD METADATA
# ============================================================

print("\n===== LOADING METADATA =====")


meta = pd.read_csv(
    METADATA_FILE
)

print(
    "Metadata shape:",
    meta.shape
)

print(
    "\nMetadata columns:"
)

for column in meta.columns:
    print(
        " -",
        column
    )


# ============================================================
# 9. IDENTIFY CELL ID
# ============================================================

print(
    "\n===== IDENTIFYING CELL ID ====="
)

cell_id_column = None

possible_cell_columns = [
    "Unnamed: 0",
    "cell",
    "cell_id",
    "barcode",
    "Cell",
    "Cell_ID"
]

for column in possible_cell_columns:

    if column in meta.columns:

        matches = (
            meta[column]
            .astype(str)
            .isin(
                set(cell_names)
            )
            .sum()
        )

        if matches == len(cell_names):

            cell_id_column = column
            break


if cell_id_column is None:

    # fallback:
    for column in meta.columns:

        matches = (
            meta[column]
            .astype(str)
            .isin(
                set(cell_names)
            )
            .sum()
        )

        if matches >= (
            0.95 * len(cell_names)
        ):

            cell_id_column = column
            break


if cell_id_column is None:

    raise ValueError(
        "Could not identify cell ID column."
    )


print(
    "Cell ID column:",
    cell_id_column
)


meta[
    cell_id_column
] = (
    meta[
        cell_id_column
    ]
    .astype(str)
    .str.strip()
)


# ============================================================
# 10. REORDER METADATA
# ============================================================

meta = (
    meta
    .set_index(
        cell_id_column
    )
    .reindex(
        cell_names
    )
)


if meta.index.isna().any():

    raise ValueError(
        "Metadata contains missing cell IDs."
    )


print(
    "Metadata aligned to expression matrix."
)


# ============================================================
# 11. REQUIRED METADATA CHECK
# ============================================================

print(
    "\n===== REQUIRED METADATA CHECK ====="
)

required_metadata = [
    "ML_class",
    "malignant_status",
    "cell_type",
    "epithelial_status",
    "tumor_candidate",
    "patient_id"
]

for column in required_metadata:

    if column in meta.columns:

        print(
            f"{column}: FOUND"
        )

    else:

        print(
            f"{column}: NOT FOUND"
        )


if "ML_class" not in meta.columns:

    raise ValueError(
        "ML_class column is required."
    )


if "patient_id" not in meta.columns:

    raise ValueError(
        "patient_id column is required."
    )


# ============================================================
# 12. DEFINE CELL GROUPS
# ============================================================

print(
    "\n===== DEFINING CELL GROUPS ====="
)


meta[
    "analysis_class"
] = np.where(

    meta[
        "ML_class"
    ]
    .astype(str)
    .eq(
        "malignant"
    ),

    "Malignant",

    np.where(

        meta[
            "ML_class"
        ]
        .astype(str)
        .eq(
            "diploid_non_malignant"
        ),

        "Diploid_non_malignant",

        "Other"

    )
)


print(
    "\nAnalysis class:"
)

print(
    meta[
        "analysis_class"
    ]
    .value_counts()
)


# ============================================================
# 13. EPITHELIAL GROUP
# ============================================================

if "epithelial_status" in meta.columns:

    meta[
        "epithelial_group"
    ] = (
        meta[
            "epithelial_status"
        ]
        .astype(str)
        .str.lower()
    )

else:

    meta[
        "epithelial_group"
    ] = "unknown"


# ============================================================
# 14. LOAD EXPRESSION MATRIX
# ============================================================

print(
    "\n===== LOADING EXPRESSION MATRIX ====="
)

from scipy.io import mmread


matrix = mmread(
    EXPRESSION_FILE
)


print(
    "Expression matrix shape:",
    matrix.shape
)


# Convert to CSC for efficient gene access
matrix = matrix.tocsc()


# ============================================================
# 15. DIMENSION CHECK
# ============================================================

if matrix.shape[0] != len(
    gene_names
):

    raise ValueError(
        "Gene dimension mismatch."
    )


if matrix.shape[1] != len(
    cell_names
):

    raise ValueError(
        "Cell dimension mismatch."
    )


print(
    "Dimension checks: PASSED"
)


# ============================================================
# 16. GENE INDEX
# ============================================================

gene_to_index = {}

for i, gene in enumerate(
    gene_names
):

    if gene not in gene_to_index:

        gene_to_index[
            gene
        ] = i


missing_genes = [
    gene
    for gene in FINAL_SIX
    if gene not in gene_to_index
]


print(
    "\n===== FINAL SIX GENE CHECK ====="
)

print(
    "Genes requested:",
    len(FINAL_SIX)
)

print(
    "Genes found:",
    len(FINAL_SIX) - len(missing_genes)
)


if missing_genes:

    print(
        "Missing genes:",
        missing_genes
    )


# ============================================================
# 17. HELPER FUNCTION
# ============================================================

def expression_vector(gene):

    idx = gene_to_index.get(
        gene
    )

    if idx is None:

        return np.zeros(
            matrix.shape[1]
        )

    return np.asarray(
        matrix[
            idx,
            :
        ].todense()
    ).ravel()


# ============================================================
# 18. HELPER: DETECTION
# ============================================================

def detection_percent(values):

    if len(values) == 0:

        return np.nan

    return (
        np.sum(
            values > 0
        )
        /
        len(values)
        *
        100
    )


# ============================================================
# 19. HELPER: MEAN EXPRESSION
# ============================================================

def mean_expression(values):

    if len(values) == 0:

        return np.nan

    return float(
        np.mean(values)
    )


# ============================================================
# 20. HELPER: LOG2 FOLD CHANGE
# ============================================================

def log2_fold_change(
    malignant,
    non_malignant
):

    if (
        len(malignant) == 0
        or
        len(non_malignant) == 0
    ):

        return np.nan

    malignant_mean = (
        np.mean(
            malignant
        )
    )

    non_malignant_mean = (
        np.mean(
            non_malignant
        )
    )

    return np.log2(
        (
            malignant_mean + 1e-6
        )
        /
        (
            non_malignant_mean + 1e-6
        )
    )


# ============================================================
# 21. CELL-TYPE DISTRIBUTION
# ============================================================

print(
    "\n===== CELL-TYPE DISTRIBUTION ====="
)

if "cell_type" in meta.columns:

    print(
        meta[
            "cell_type"
        ]
        .fillna(
            "Unknown"
        )
        .value_counts()
        .to_string()
    )


# ============================================================
# 22. ANALYZE SIX GENES
# ============================================================

print(
    "\n===== ANALYZING FINAL SIX GENES ====="
)


results = []


for gene in FINAL_SIX:

    print(
        f"\n--- {gene} ---"
    )


    values = expression_vector(
        gene
    )


    malignant_mask = (
        meta[
            "analysis_class"
        ]
        .values
        == "Malignant"
    )


    diploid_mask = (
        meta[
            "analysis_class"
        ]
        .values
        == "Diploid_non_malignant"
    )


    other_mask = (
        meta[
            "analysis_class"
        ]
        .values
        == "Other"
    )


    malignant_values = (
        values[
            malignant_mask
        ]
    )


    diploid_values = (
        values[
            diploid_mask
        ]
    )


    other_values = (
        values[
            other_mask
        ]
    )


    malignant_detection = (
        detection_percent(
            malignant_values
        )
    )


    diploid_detection = (
        detection_percent(
            diploid_values
        )
    )


    other_detection = (
        detection_percent(
            other_values
        )
    )


    detection_difference = (
        malignant_detection
        -
        diploid_detection
    )


    log2fc = log2_fold_change(
        malignant_values,
        diploid_values
    )


    # --------------------------------------------------------
    # EPITHELIAL ANALYSIS
    # --------------------------------------------------------

    epithelial_detection = np.nan

    non_epithelial_detection = np.nan

    epithelial_difference = np.nan

    if "cell_type" in meta.columns:

        cell_type = (
            meta[
                "cell_type"
            ]
            .fillna(
                "Unknown"
            )
            .astype(str)
            .str.lower()
            .values
        )


        epithelial_mask = np.array([
            (
                "epithelial"
                in x
            )
            or
            (
                "luminal"
                in x
            )
            or
            (
                "tumor"
                in x
            )
            or
            (
                "cancer"
                in x
            )

            for x in cell_type
        ])


        non_epithelial_mask = (
            ~epithelial_mask
        )


        epithelial_values = (
            values[
                epithelial_mask
            ]
        )


        non_epithelial_values = (
            values[
                non_epithelial_mask
            ]
        )


        epithelial_detection = (
            detection_percent(
                epithelial_values
            )
        )


        non_epithelial_detection = (
            detection_percent(
                non_epithelial_values
            )
        )


        epithelial_difference = (
            epithelial_detection
            -
            non_epithelial_detection
        )


    # --------------------------------------------------------
    # MALIGNANT EPITHELIAL
    # --------------------------------------------------------

    malignant_epithelial_detection = np.nan


    if "cell_type" in meta.columns:

        cell_type = (
            meta[
                "cell_type"
            ]
            .fillna(
                "Unknown"
            )
            .astype(str)
            .str.lower()
            .values
        )


        epithelial_mask = np.array([
            (
                "epithelial"
                in x
            )
            or
            (
                "luminal"
                in x
            )
            or
            (
                "tumor"
                in x
            )
            or
            (
                "cancer"
                in x
            )

            for x in cell_type
        ])


        malignant_epithelial_mask = (
            malignant_mask
            &
            epithelial_mask
        )


        malignant_epithelial_values = (
            values[
                malignant_epithelial_mask
            ]
        )


        malignant_epithelial_detection = (
            detection_percent(
                malignant_epithelial_values
            )
        )


    # --------------------------------------------------------
    # SAMPLE CONSISTENCY
    # --------------------------------------------------------

    sample_log2fcs = []

    sample_positive = 0

    samples_tested = 0


    sample_ids = (
        meta[
            "patient_id"
        ]
        .dropna()
        .unique()
    )


    for sample in sample_ids:

        sample_mask = (
            meta[
                "patient_id"
            ]
            .values
            == sample
        )


        sample_malignant = (
            sample_mask
            &
            malignant_mask
        )


        sample_diploid = (
            sample_mask
            &
            diploid_mask
        )


        if (
            sample_malignant.sum()
            == 0
            or
            sample_diploid.sum()
            == 0
        ):

            continue


        sample_malignant_values = (
            values[
                sample_malignant
            ]
        )


        sample_diploid_values = (
            values[
                sample_diploid
            ]
        )


        sample_fc = log2_fold_change(
            sample_malignant_values,
            sample_diploid_values
        )


        if np.isfinite(
            sample_fc
        ):

            samples_tested += 1

            sample_log2fcs.append(
                sample_fc
            )

            if sample_fc > 0:

                sample_positive += 1


    if samples_tested > 0:

        sample_consistency = (
            sample_positive
            /
            samples_tested
        )

        mean_sample_log2fc = (
            np.mean(
                sample_log2fcs
            )
        )

    else:

        sample_consistency = np.nan

        mean_sample_log2fc = np.nan


    # --------------------------------------------------------
    # SPECIFICITY CLASSIFICATION
    # --------------------------------------------------------

    if (
        detection_difference >= 30
        and
        epithelial_difference >= 20
        and
        sample_consistency >= 0.67
    ):

        specificity_class = (
            "Strong_tumor_epithelial_specificity"
        )


    elif (
        detection_difference >= 20
        and
        sample_consistency >= 0.67
    ):

        specificity_class = (
            "Tumor_associated"
        )


    elif (
        detection_difference >= 10
    ):

        specificity_class = (
            "Weak_tumor_association"
        )


    else:

        specificity_class = (
            "Broad_or_unresolved_expression"
        )


    # --------------------------------------------------------
    # STORE
    # --------------------------------------------------------

    results.append({

        "gene": gene,

        "malignant_detection_pct":
            malignant_detection,

        "diploid_detection_pct":
            diploid_detection,

        "other_detection_pct":
            other_detection,

        "malignant_detection_difference_pct":
            detection_difference,

        "malignant_vs_diploid_log2FC":
            log2fc,

        "malignant_mean_expression":
            mean_expression(
                malignant_values
            ),

        "diploid_mean_expression":
            mean_expression(
                diploid_values
            ),

        "epithelial_detection_pct":
            epithelial_detection,

        "non_epithelial_detection_pct":
            non_epithelial_detection,

        "epithelial_detection_difference_pct":
            epithelial_difference,

        "malignant_epithelial_detection_pct":
            malignant_epithelial_detection,

        "samples_tested":
            samples_tested,

        "samples_positive_log2FC":
            sample_positive,

        "sample_consistency":
            sample_consistency,

        "mean_sample_log2FC":
            mean_sample_log2fc,

        "tumor_specificity_class":
            specificity_class

    })


# ============================================================
# 23. RESULTS DATAFRAME
# ============================================================

result_df = pd.DataFrame(
    results
)


# ============================================================
# 24. MERGE STEP 42
# ============================================================

step42_small_columns = [
    "gene",
    "external_support_category",
    "external_validation_score",
    "evidence_layers",
    "final_external_validation_class"
]

available_step42_columns = [
    c
    for c in step42_small_columns
    if c in pd.read_csv(
        STEP42_FILE,
        nrows=1
    ).columns
]


step42 = pd.read_csv(
    STEP42_FILE
)

step42["gene"] = (
    step42["gene"]
    .astype(str)
    .str.strip()
)


step42_small = step42[
    available_step42_columns
].copy()


result_df = result_df.merge(
    step42_small,
    on="gene",
    how="left"
)


# ============================================================
# 25. CREATE FINAL INTERPRETATION
# ============================================================

def final_interpretation(row):

    cls = row[
        "tumor_specificity_class"
    ]


    external = str(
        row.get(
            "external_support_category",
            ""
        )
    )


    if (
        cls
        ==
        "Strong_tumor_epithelial_specificity"
        and
        external
        in [
            "Strong_external_expression_support",
            "External_expression_supported"
        ]
    ):

        return (
            "High-priority tumor-associated candidate"
        )


    if (
        cls
        in [
            "Strong_tumor_epithelial_specificity",
            "Tumor_associated"
        ]
    ):

        return (
            "Tumor-associated candidate"
        )


    if (
        cls
        ==
        "Weak_tumor_association"
    ):

        return (
            "Weak tumor association"
        )


    return (
        "Requires further specificity validation"
    )


result_df[
    "final_specificity_interpretation"
] = result_df.apply(
    final_interpretation,
    axis=1
)


# ============================================================
# 26. SORT
# ============================================================

result_df = result_df.sort_values(
    [
        "malignant_detection_difference_pct",
        "sample_consistency"
    ],
    ascending=[
        False,
        False
    ]
).reset_index(
    drop=True
)


result_df[
    "specificity_rank"
] = (
    result_df.index + 1
)


# ============================================================
# 27. SAVE RESULT
# ============================================================

output_file = os.path.join(
    OUTPUT_DIR,
    "43_celltype_tumor_specificity.csv"
)

result_df.to_csv(
    output_file,
    index=False
)


# ============================================================
# 28. PRINT RESULTS
# ============================================================

print(
    "\n" + "=" * 70
)

print(
    "STEP 43 - TUMOR SPECIFICITY RESULTS"
)

print(
    "=" * 70
)


display_columns = [
    "specificity_rank",
    "gene",
    "malignant_detection_pct",
    "diploid_detection_pct",
    "malignant_detection_difference_pct",
    "malignant_vs_diploid_log2FC",
    "epithelial_detection_pct",
    "malignant_epithelial_detection_pct",
    "samples_tested",
    "samples_positive_log2FC",
    "sample_consistency",
    "tumor_specificity_class",
    "final_specificity_interpretation"
]


display_columns = [
    c
    for c in display_columns
    if c in result_df.columns
]


print(
    result_df[
        display_columns
    ].to_string(
        index=False
    )
)


# ============================================================
# 29. SUMMARY
# ============================================================

print(
    "\n===== SPECIFICITY SUMMARY ====="
)


print(
    result_df[
        "tumor_specificity_class"
    ]
    .value_counts()
    .to_string()
)


# ============================================================
# 30. SAVE SUMMARY
# ============================================================

summary_file = os.path.join(
    OUTPUT_DIR,
    "43_tumor_specificity_summary.txt"
)


with open(
    summary_file,
    "w"
) as f:

    f.write(
        "STEP 43 - CELL-TYPE / TUMOR-SPECIFICITY VALIDATION\n"
    )

    f.write(
        "=" * 70 + "\n\n"
    )

    f.write(
        "Purpose:\n"
    )

    f.write(
        "Determine whether the six consensus candidates "
        "are preferentially associated with malignant "
        "and tumor-epithelial cells.\n\n"
    )

    f.write(
        "RESULTS\n"
    )

    f.write(
        "-" * 70 + "\n"
    )


    for _, row in result_df.iterrows():

        f.write(
            f"\nGene: {row['gene']}\n"
        )

        f.write(
            f"Malignant detection: "
            f"{row['malignant_detection_pct']:.2f}%\n"
        )

        f.write(
            f"Diploid detection: "
            f"{row['diploid_detection_pct']:.2f}%\n"
        )

        f.write(
            f"Detection difference: "
            f"{row['malignant_detection_difference_pct']:.2f} percentage points\n"
        )

        f.write(
            f"Malignant vs diploid log2FC: "
            f"{row['malignant_vs_diploid_log2FC']:.4f}\n"
        )

        f.write(
            f"Epithelial detection: "
            f"{row['epithelial_detection_pct']:.2f}%\n"
        )

        f.write(
            f"Malignant epithelial detection: "
            f"{row['malignant_epithelial_detection_pct']:.2f}%\n"
        )

        f.write(
            f"Sample consistency: "
            f"{row['sample_consistency']:.3f}\n"
        )

        f.write(
            f"Tumor specificity class: "
            f"{row['tumor_specificity_class']}\n"
        )

        f.write(
            f"Final interpretation: "
            f"{row['final_specificity_interpretation']}\n"
        )


    f.write(
        "\n\nIMPORTANT LIMITATION\n"
    )

    f.write(
        "-" * 70 + "\n"
    )

    f.write(
        "Tumor-associated expression does not by itself "
        "establish diagnostic, prognostic, therapeutic, "
        "or clinical biomarker validity.\n"
    )


# ============================================================
# 31. FIGURE 1
# MALIGNANT VS DIPLOID DETECTION
# ============================================================

print(
    "\n===== GENERATING FIGURES ====="
)


fig, ax = plt.subplots(
    figsize=(9, 6)
)


x = np.arange(
    len(result_df)
)

width = 0.36


ax.bar(
    x - width / 2,
    result_df[
        "malignant_detection_pct"
    ],
    width,
    label="Malignant"
)


ax.bar(
    x + width / 2,
    result_df[
        "diploid_detection_pct"
    ],
    width,
    label="Diploid/non-malignant"
)


ax.set_xticks(
    x
)

ax.set_xticklabels(
    result_df[
        "gene"
    ]
)


ax.set_ylabel(
    "Detection (%)"
)

ax.set_xlabel(
    "Gene"
)

ax.set_title(
    "Final Consensus Genes: Malignant vs Diploid Detection"
)

ax.legend()

plt.tight_layout()


figure1 = os.path.join(
    FIGURE_DIR,
    "43_malignant_vs_diploid_detection.png"
)


plt.savefig(
    figure1,
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# ============================================================
# 32. FIGURE 2
# MALIGNANT DETECTION DIFFERENCE
# ============================================================

fig, ax = plt.subplots(
    figsize=(9, 6)
)


ax.bar(
    result_df[
        "gene"
    ],
    result_df[
        "malignant_detection_difference_pct"
    ]
)


ax.axhline(
    0,
    linewidth=1
)


ax.set_ylabel(
    "Detection difference (percentage points)"
)

ax.set_xlabel(
    "Gene"
)

ax.set_title(
    "Malignant Detection Enrichment"
)

plt.tight_layout()


figure2 = os.path.join(
    FIGURE_DIR,
    "43_malignant_detection_difference.png"
)


plt.savefig(
    figure2,
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# ============================================================
# 33. FIGURE 3
# MALIGNANT VS DIPLOID LOG2FC
# ============================================================

fig, ax = plt.subplots(
    figsize=(9, 6)
)


ax.bar(
    result_df[
        "gene"
    ],
    result_df[
        "malignant_vs_diploid_log2FC"
    ]
)


ax.axhline(
    0,
    linewidth=1
)


ax.set_ylabel(
    "log2 fold-change"
)

ax.set_xlabel(
    "Gene"
)

ax.set_title(
    "Malignant vs Diploid Expression"
)

plt.tight_layout()


figure3 = os.path.join(
    FIGURE_DIR,
    "43_malignant_vs_diploid_log2FC.png"
)


plt.savefig(
    figure3,
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# ============================================================
# 34. COMPLETION
# ============================================================

print(
    "\n" + "=" * 70
)

print(
    "STEP 43 COMPLETED"
)

print(
    "=" * 70
)

print(
    "\nMain result:"
)

print(
    output_file
)

print(
    "\nSummary:"
)

print(
    summary_file
)

print(
    "\nFigures:"
)

print(
    figure1
)

print(
    figure2
)

print(
    figure3
)

print(
    "\nResults directory:"
)

print(
    OUTPUT_DIR
)

print(
    "\n" + "=" * 70
)