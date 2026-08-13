# ============================================================
# STEP 38
# CONSENSUS CANDIDATE EXPRESSION VALIDATION
# ============================================================
#
# Purpose:
# Validate the 6 Bioinformatics-ML consensus genes using the
# ACTUAL single-cell expression matrix.
#
# Consensus genes:
#   PRSS23
#   FCGRT
#   ETFB
#   COX8A
#   EPN3
#   MRPS7
#
# Informative samples:
#   BC11
#   BC12
#   BC17
#
# Comparison:
#   malignant
#   diploid_non_malignant
#
# Metrics:
#   - Mean expression
#   - Detection percentage
#   - Log2 fold-change
#   - Detection difference
#   - Sample-level consistency
#
# Important:
# This is expression validation, NOT clinical validation.
# The validation score is exploratory.
#
# ============================================================

import os
import gzip
import numpy as np
import pandas as pd

from scipy import sparse
from scipy.io import mmread


# ============================================================
# 1. CONFIGURATION
# ============================================================

EXPRESSION_MATRIX = (
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

CONSENSUS_FILE = (
    "results/ai_ml/36_bioinformatics_ml_consensus/"
    "36_strong_consensus_candidates.csv"
)

STEP37_FILE = (
    "results/ai_ml/37_consensus_characterization/"
    "37_consensus_candidate_characterization.csv"
)

OUTPUT_DIR = (
    "results/ai_ml/38_consensus_expression_validation"
)

INFORMATIVE_SAMPLES = [
    "BC11",
    "BC12",
    "BC17"
]


# ============================================================
# 2. CREATE OUTPUT DIRECTORY
# ============================================================

os.makedirs(
    OUTPUT_DIR,
    exist_ok=True
)


# ============================================================
# 3. START
# ============================================================

print("=" * 70)
print("STEP 38 - CONSENSUS CANDIDATE EXPRESSION VALIDATION")
print("=" * 70)


# ============================================================
# 4. INPUT FILE CHECK
# ============================================================

print("\n===== INPUT FILE CHECK =====")

required_files = {
    "Expression matrix": EXPRESSION_MATRIX,
    "Genes": GENES_FILE,
    "Cells": CELLS_FILE,
    "Metadata": METADATA_FILE,
    "Consensus candidates": CONSENSUS_FILE
}


for name, path in required_files.items():

    if os.path.exists(path):

        print(
            f"{name}: FOUND"
        )

    else:

        raise FileNotFoundError(
            f"\n{name} file not found:\n{path}"
        )


if os.path.exists(STEP37_FILE):

    print(
        "Step 37 characterization: FOUND"
    )

else:

    print(
        "Step 37 characterization: NOT FOUND"
    )

    print(
        "Continuing without Step 37 merge."
    )


# ============================================================
# 5. LOAD CONSENSUS GENES
# ============================================================

print(
    "\n===== LOADING CONSENSUS GENES ====="
)

consensus = pd.read_csv(
    CONSENSUS_FILE
)


if "gene" not in consensus.columns:

    raise ValueError(
        "Consensus file does not contain "
        "a 'gene' column."
    )


consensus["gene"] = (
    consensus["gene"]
    .astype(str)
    .str.strip()
)


consensus_genes = (
    consensus["gene"]
    .dropna()
    .drop_duplicates()
    .tolist()
)


print(
    "Consensus genes:",
    len(consensus_genes)
)


print(
    "\nGenes:"
)

for gene in consensus_genes:

    print(
        " -",
        gene
    )


# ============================================================
# 6. LOAD GENE NAMES
# ============================================================

print(
    "\n===== LOADING GENE NAMES ====="
)

genes_df = pd.read_csv(
    GENES_FILE,
    sep="\t",
    header=None,
    dtype=str
)


gene_names = (
    genes_df.iloc[:, 0]
    .astype(str)
    .str.strip()
    .tolist()
)


print(
    "Genes in expression matrix:",
    len(gene_names)
)


# ============================================================
# 7. CHECK FOR DUPLICATE GENE NAMES
# ============================================================

duplicate_genes = (
    pd.Series(gene_names)
    .duplicated()
    .sum()
)


print(
    "Duplicate gene names:",
    duplicate_genes
)


# ============================================================
# 8. LOAD CELL NAMES
# ============================================================

print(
    "\n===== LOADING CELL NAMES ====="
)

cells_df = pd.read_csv(
    CELLS_FILE,
    sep="\t",
    header=None,
    dtype=str
)


cell_names = (
    cells_df.iloc[:, 0]
    .astype(str)
    .str.strip()
    .tolist()
)


print(
    "Cells in expression matrix:",
    len(cell_names)
)


# ============================================================
# 9. LOAD METADATA
# ============================================================

print(
    "\n===== LOADING METADATA ====="
)

metadata = pd.read_csv(
    METADATA_FILE
)


print(
    "Metadata dimensions:",
    metadata.shape
)


print(
    "\nMetadata columns:"
)

for column in metadata.columns:

    print(
        " ",
        column
    )


# ============================================================
# 10. IDENTIFY CELL ID COLUMN
# ============================================================
#
# Your current metadata contains:
#
#   Unnamed: 0
#
# This is most likely the original Seurat cell barcode/index.
#
# We first test known column names.
# If none match, we automatically search all columns and
# identify the column with the greatest overlap with the
# expression-matrix cell names.
#
# ============================================================

print(
    "\n===== IDENTIFYING CELL ID COLUMN ====="
)


possible_cell_columns = [

    "cell",

    "cell_id",

    "Cell",

    "Cell_ID",

    "barcode",

    "Barcode",

    "cell_barcode",

    "Cell_Barcode",

    "Unnamed: 0"

]


cell_column = None


expression_cell_set = set(
    cell_names
)


# ------------------------------------------------------------
# First pass: known cell-ID columns
# ------------------------------------------------------------

for column in possible_cell_columns:

    if column in metadata.columns:

        values = (
            metadata[column]
            .astype(str)
            .str.strip()
        )

        match_count = (
            values.isin(
                expression_cell_set
            ).sum()
        )

        if match_count > 0:

            cell_column = column

            print(
                "Detected cell ID column:",
                cell_column
            )

            print(
                "Matching expression cells:",
                match_count
            )

            break


# ------------------------------------------------------------
# Second pass: automatically search all columns
# ------------------------------------------------------------

if cell_column is None:

    print(
        "\nKnown cell-ID columns did not match."
    )

    print(
        "Searching all metadata columns..."
    )


    column_match_counts = {}


    for column in metadata.columns:

        values = (
            metadata[column]
            .astype(str)
            .str.strip()
        )

        match_count = (
            values.isin(
                expression_cell_set
            ).sum()
        )

        column_match_counts[
            column
        ] = match_count


    best_column = max(
        column_match_counts,
        key=column_match_counts.get
    )


    best_match_count = (
        column_match_counts[
            best_column
        ]
    )


    if best_match_count > 0:

        cell_column = best_column

        print(
            "Automatically detected cell ID column:",
            cell_column
        )

        print(
            "Matching expression cells:",
            best_match_count
        )

    else:

        print(
            "\nColumn matching results:"
        )

        for column, count in (
            sorted(
                column_match_counts.items(),
                key=lambda x: x[1],
                reverse=True
            )
        ):

            print(
                f"{column}: {count}"
            )


        raise ValueError(
            "\nCould not identify a metadata "
            "cell ID column."
        )


# ============================================================
# 11. STANDARDIZE CELL ID
# ============================================================

metadata["cell_id"] = (
    metadata[
        cell_column
    ]
    .astype(str)
    .str.strip()
)


# ============================================================
# 12. VERIFY CELL IDs
# ============================================================

metadata_cell_set = set(
    metadata["cell_id"]
)


matched_cells = (
    metadata["cell_id"]
    .isin(
        expression_cell_set
    )
    .sum()
)


print(
    "\n===== CELL ID VALIDATION ====="
)

print(
    "Metadata cells:",
    len(metadata)
)

print(
    "Expression cells:",
    len(cell_names)
)

print(
    "Matching metadata cells:",
    matched_cells
)


if matched_cells == 0:

    raise ValueError(
        "No metadata cell IDs match the "
        "expression matrix."
    )


print(
    "Cell ID validation: PASSED"
)


# ============================================================
# 13. CHECK ML CLASS COLUMN
# ============================================================

print(
    "\n===== CHECKING ML CLASS ====="
)


if "ML_class" not in metadata.columns:

    raise ValueError(
        "\n'ML_class' column not found."
    )


print(
    metadata[
        "ML_class"
    ].value_counts()
)


# ============================================================
# 14. CHECK SAMPLE COLUMN
# ============================================================

print(
    "\n===== IDENTIFYING SAMPLE COLUMN ====="
)


possible_sample_columns = [

    "patient_id",

    "sample_id",

    "orig.ident",

    "sample"

]


sample_column = None


for column in possible_sample_columns:

    if column in metadata.columns:

        sample_column = column

        break


if sample_column is None:

    raise ValueError(
        "\nCould not identify sample column."
    )


print(
    "Sample column:",
    sample_column
)


metadata[
    sample_column
] = (
    metadata[
        sample_column
    ]
    .astype(str)
    .str.strip()
)


# ============================================================
# 15. CHECK SAMPLE DISTRIBUTION
# ============================================================

print(
    "\n===== SAMPLE DISTRIBUTION ====="
)

print(
    pd.crosstab(
        metadata[
            sample_column
        ],
        metadata[
            "ML_class"
        ]
    )
)


# ============================================================
# 16. SELECT INFORMATIVE SAMPLES
# ============================================================

print(
    "\n===== INFORMATIVE SAMPLES ====="
)


metadata = metadata[
    metadata[
        sample_column
    ].isin(
        INFORMATIVE_SAMPLES
    )
].copy()


if len(metadata) == 0:

    raise ValueError(
        "No cells found for BC11, BC12, or BC17."
    )


print(
    "Samples:",
    sorted(
        metadata[
            sample_column
        ].unique()
    )
)


print(
    "\nCells per sample:"
)

print(
    metadata[
        sample_column
    ].value_counts()
)


print(
    "\nClass distribution:"
)

print(
    metadata[
        "ML_class"
    ].value_counts()
)


# ============================================================
# 17. LOAD EXPRESSION MATRIX
# ============================================================

print(
    "\n===== LOADING EXPRESSION MATRIX ====="
)


expression = mmread(
    EXPRESSION_MATRIX
)


if not sparse.issparse(
    expression
):

    expression = sparse.coo_matrix(
        expression
    )


expression = expression.tocsr()


print(
    "Expression matrix shape:",
    expression.shape
)


# ============================================================
# 18. DIMENSION VALIDATION
# ============================================================

n_matrix_genes = (
    expression.shape[0]
)

n_matrix_cells = (
    expression.shape[1]
)


if n_matrix_genes != len(
    gene_names
):

    raise ValueError(
        "\nGene dimension mismatch.\n"
        f"Expression matrix: {n_matrix_genes}\n"
        f"Genes file: {len(gene_names)}"
    )


if n_matrix_cells != len(
    cell_names
):

    raise ValueError(
        "\nCell dimension mismatch.\n"
        f"Expression matrix: {n_matrix_cells}\n"
        f"Cells file: {len(cell_names)}"
    )


print(
    "Dimension checks: PASSED"
)


# ============================================================
# 19. CREATE GENE INDEX
# ============================================================

gene_to_index = {}


for index, gene in enumerate(
    gene_names
):

    # Keep first occurrence
    if gene not in gene_to_index:

        gene_to_index[
            gene
        ] = index


# ============================================================
# 20. CHECK CONSENSUS GENES
# ============================================================

print(
    "\n===== CONSENSUS GENE CHECK ====="
)


present_genes = []

missing_genes = []


for gene in consensus_genes:

    if gene in gene_to_index:

        present_genes.append(
            gene
        )

    else:

        missing_genes.append(
            gene
        )


print(
    "Consensus genes:",
    len(consensus_genes)
)

print(
    "Present in expression matrix:",
    len(present_genes)
)

print(
    "Missing:",
    len(missing_genes)
)


if missing_genes:

    print(
        "\nMissing genes:"
    )

    for gene in missing_genes:

        print(
            " -",
            gene
        )


if len(present_genes) == 0:

    raise ValueError(
        "None of the consensus genes are present "
        "in the expression matrix."
    )


# ============================================================
# 21. CREATE CELL INDEX
# ============================================================

cell_to_index = {
    cell: index
    for index, cell in enumerate(
        cell_names
    )
}


metadata["matrix_index"] = (
    metadata[
        "cell_id"
    ]
    .map(
        cell_to_index
    )
)


missing_matrix_cells = (
    metadata[
        "matrix_index"
    ].isna().sum()
)


print(
    "\nMetadata cells not found in matrix:",
    missing_matrix_cells
)


if missing_matrix_cells > 0:

    metadata = metadata[
        metadata[
            "matrix_index"
        ].notna()
    ].copy()


metadata[
    "matrix_index"
] = (
    metadata[
        "matrix_index"
    ]
    .astype(int)
)


# ============================================================
# 22. SAMPLE-LEVEL EXPRESSION ANALYSIS
# ============================================================

print(
    "\n" + "=" * 70
)

print(
    "SAMPLE-LEVEL EXPRESSION VALIDATION"
)

print(
    "=" * 70
)


sample_results = []


for sample in INFORMATIVE_SAMPLES:

    sample_metadata = metadata[
        metadata[
            sample_column
        ] == sample
    ].copy()


    malignant_metadata = (
        sample_metadata[
            sample_metadata[
                "ML_class"
            ] == "malignant"
        ]
    )


    diploid_metadata = (
        sample_metadata[
            sample_metadata[
                "ML_class"
            ]
            ==
            "diploid_non_malignant"
        ]
    )


    malignant_indices = (
        malignant_metadata[
            "matrix_index"
        ]
        .astype(int)
        .tolist()
    )


    diploid_indices = (
        diploid_metadata[
            "matrix_index"
        ]
        .astype(int)
        .tolist()
    )


    print(
        f"\n--- SAMPLE: {sample} ---"
    )


    print(
        "Malignant cells:",
        len(malignant_indices)
    )


    print(
        "Diploid/non-malignant cells:",
        len(diploid_indices)
    )


    if len(malignant_indices) == 0:

        print(
            "No malignant cells. Skipping."
        )

        continue


    for gene in present_genes:

        gene_index = (
            gene_to_index[
                gene
            ]
        )


        # ----------------------------------------------------
        # Extract malignant expression
        # ----------------------------------------------------

        malignant_values = (
            expression[
                gene_index,
                malignant_indices
            ]
            .toarray()
            .flatten()
        )


        # ----------------------------------------------------
        # Extract diploid expression
        # ----------------------------------------------------

        if len(
            diploid_indices
        ) > 0:

            diploid_values = (
                expression[
                    gene_index,
                    diploid_indices
                ]
                .toarray()
                .flatten()
            )

        else:

            diploid_values = np.array(
                [],
                dtype=float
            )


        # ----------------------------------------------------
        # Mean expression
        # ----------------------------------------------------

        malignant_mean = (
            np.mean(
                malignant_values
            )
            if len(
                malignant_values
            ) > 0
            else np.nan
        )


        diploid_mean = (
            np.mean(
                diploid_values
            )
            if len(
                diploid_values
            ) > 0
            else np.nan
        )


        # ----------------------------------------------------
        # Detection percentage
        # ----------------------------------------------------

        malignant_detection = (
            np.mean(
                malignant_values > 0
            )
            * 100
            if len(
                malignant_values
            ) > 0
            else np.nan
        )


        diploid_detection = (
            np.mean(
                diploid_values > 0
            )
            * 100
            if len(
                diploid_values
            ) > 0
            else np.nan
        )


        # ----------------------------------------------------
        # Log2 fold change
        #
        # Pseudocount = 0.01
        # ----------------------------------------------------

        if (
            not np.isnan(
                malignant_mean
            )
            and
            not np.isnan(
                diploid_mean
            )
        ):

            log2fc = np.log2(
                (
                    malignant_mean
                    + 0.01
                )
                /
                (
                    diploid_mean
                    + 0.01
                )
            )

        else:

            log2fc = np.nan


        # ----------------------------------------------------
        # Detection difference
        # ----------------------------------------------------

        if (
            not np.isnan(
                malignant_detection
            )
            and
            not np.isnan(
                diploid_detection
            )
        ):

            detection_difference = (
                malignant_detection
                -
                diploid_detection
            )

        else:

            detection_difference = np.nan


        # ----------------------------------------------------
        # Store
        # ----------------------------------------------------

        sample_results.append({

            "sample":
                sample,

            "gene":
                gene,

            "malignant_cells":
                len(
                    malignant_values
                ),

            "diploid_cells":
                len(
                    diploid_values
                ),

            "malignant_mean_expression":
                malignant_mean,

            "diploid_mean_expression":
                diploid_mean,

            "log2FC_malignant_vs_diploid":
                log2fc,

            "malignant_detection_pct":
                malignant_detection,

            "diploid_detection_pct":
                diploid_detection,

            "detection_difference_pct":
                detection_difference

        })


# ============================================================
# 23. CREATE SAMPLE RESULT DATAFRAME
# ============================================================

sample_results_df = pd.DataFrame(
    sample_results
)


if len(
    sample_results_df
) == 0:

    raise ValueError(
        "No sample-level results were generated."
    )


# ============================================================
# 24. SAVE SAMPLE-LEVEL RESULTS
# ============================================================

sample_output = os.path.join(
    OUTPUT_DIR,
    "38_sample_level_expression_validation.csv"
)


sample_results_df.to_csv(
    sample_output,
    index=False
)


print(
    "\nSample-level results saved:"
)

print(
    sample_output
)


# ============================================================
# 25. AGGREGATE ACROSS SAMPLES
# ============================================================

print(
    "\n===== AGGREGATING ACROSS SAMPLES ====="
)


summary_rows = []


for gene in present_genes:

    gene_df = (
        sample_results_df[
            sample_results_df[
                "gene"
            ] == gene
        ]
        .copy()
    )


    n_samples = len(
        gene_df
    )


    mean_malignant = (
        gene_df[
            "malignant_mean_expression"
        ].mean()
    )


    mean_diploid = (
        gene_df[
            "diploid_mean_expression"
        ].mean()
    )


    mean_log2fc = (
        gene_df[
            "log2FC_malignant_vs_diploid"
        ].mean()
    )


    mean_malignant_detection = (
        gene_df[
            "malignant_detection_pct"
        ].mean()
    )


    mean_diploid_detection = (
        gene_df[
            "diploid_detection_pct"
        ].mean()
    )


    mean_detection_difference = (
        gene_df[
            "detection_difference_pct"
        ].mean()
    )


    # --------------------------------------------------------
    # Sample consistency
    # --------------------------------------------------------

    positive_log2fc_count = (
        (
            gene_df[
                "log2FC_malignant_vs_diploid"
            ]
            > 0
        )
        .sum()
    )


    log2fc_ge_1_count = (
        (
            gene_df[
                "log2FC_malignant_vs_diploid"
            ]
            >= 1
        )
        .sum()
    )


    positive_detection_count = (
        (
            gene_df[
                "detection_difference_pct"
            ]
            > 0
        )
        .sum()
    )


    # --------------------------------------------------------
    # Consistency category
    # --------------------------------------------------------

    if (
        positive_log2fc_count
        == n_samples
    ):

        consistency = (
            "Consistent_positive_expression"
        )

    elif (
        positive_log2fc_count
        >= 2
    ):

        consistency = (
            "Mostly_positive_expression"
        )

    elif (
        positive_log2fc_count
        == 1
    ):

        consistency = (
            "Sample_specific"
        )

    else:

        consistency = (
            "No_positive_expression_difference"
        )


    # --------------------------------------------------------
    # Consistency score
    # --------------------------------------------------------

    if n_samples > 0:

        consistency_score = (
            positive_log2fc_count
            /
            n_samples
        )

    else:

        consistency_score = np.nan


    summary_rows.append({

        "gene":
            gene,

        "samples_tested":
            n_samples,

        "mean_malignant_expression":
            mean_malignant,

        "mean_diploid_expression":
            mean_diploid,

        "mean_log2FC":
            mean_log2fc,

        "mean_malignant_detection_pct":
            mean_malignant_detection,

        "mean_diploid_detection_pct":
            mean_diploid_detection,

        "mean_detection_difference_pct":
            mean_detection_difference,

        "samples_with_positive_log2FC":
            positive_log2fc_count,

        "samples_with_log2FC_ge_1":
            log2fc_ge_1_count,

        "samples_with_positive_detection_difference":
            positive_detection_count,

        "sample_consistency_score":
            consistency_score,

        "expression_consistency":
            consistency

    })


summary_df = pd.DataFrame(
    summary_rows
)


# ============================================================
# 26. NORMALIZATION FUNCTION
# ============================================================

def minmax_normalize(series):

    values = pd.to_numeric(
        series,
        errors="coerce"
    )

    result = pd.Series(
        np.nan,
        index=series.index,
        dtype=float
    )

    valid = values.notna()

    if valid.sum() == 0:

        return result


    minimum = (
        values[valid].min()
    )

    maximum = (
        values[valid].max()
    )


    if maximum == minimum:

        result.loc[
            valid
        ] = 1.0

    else:

        result.loc[
            valid
        ] = (
            values[valid]
            - minimum
        ) / (
            maximum
            - minimum
        )


    return result


# ============================================================
# 27. CREATE EXPLORATORY EXPRESSION SCORE
# ============================================================

summary_df[
    "log2FC_score"
] = minmax_normalize(
    summary_df[
        "mean_log2FC"
    ]
)


summary_df[
    "detection_difference_score"
] = minmax_normalize(
    summary_df[
        "mean_detection_difference_pct"
    ]
)


#
# Exploratory score:
#
# 50% sample consistency
# 30% mean log2FC
# 20% detection difference
#

summary_df[
    "expression_validation_score"
] = (

    0.50
    *
    summary_df[
        "sample_consistency_score"
    ].fillna(0)

    +

    0.30
    *
    summary_df[
        "log2FC_score"
    ].fillna(0)

    +

    0.20
    *
    summary_df[
        "detection_difference_score"
    ].fillna(0)

)


# ============================================================
# 28. RANK EXPRESSION VALIDATION
# ============================================================

summary_df = (
    summary_df
    .sort_values(
        [
            "expression_validation_score",
            "mean_log2FC",
            "mean_detection_difference_pct"
        ],
        ascending=[
            False,
            False,
            False
        ],
        na_position="last"
    )
    .reset_index(
        drop=True
    )
)


summary_df[
    "expression_validation_rank"
] = (
    summary_df.index + 1
)


# ============================================================
# 29. SAVE SUMMARY
# ============================================================

summary_output = os.path.join(
    OUTPUT_DIR,
    "38_consensus_expression_summary.csv"
)


summary_df.to_csv(
    summary_output,
    index=False
)


# ============================================================
# 30. INTEGRATE WITH STEP 37
# ============================================================

integrated_output = None


if os.path.exists(
    STEP37_FILE
):

    print(
        "\n===== INTEGRATING STEP 37 + STEP 38 ====="
    )


    step37 = pd.read_csv(
        STEP37_FILE
    )


    if "gene" in step37.columns:

        step37["gene"] = (
            step37["gene"]
            .astype(str)
            .str.strip()
        )


        step37_columns = [

            "gene",

            "consensus_rank",

            "exploratory_consensus_score",

            "ML_folds_selected",

            "ML_recurrence",

            "external_expression_prevalence",

            "external_mean_expression",

            "external_support",

            "biological_category"

        ]


        step37_columns = [

            column
            for column in step37_columns
            if column in step37.columns

        ]


        step37_small = (
            step37[
                step37_columns
            ]
            .drop_duplicates(
                subset=["gene"]
            )
        )


        integrated = step37_small.merge(
            summary_df,
            on="gene",
            how="left",
            suffixes=(
                "_step37",
                "_step38"
            )
        )


        integrated = (
            integrated
            .sort_values(
                "expression_validation_score",
                ascending=False,
                na_position="last"
            )
            .reset_index(
                drop=True
            )
        )


        integrated[
            "final_expression_validation_rank"
        ] = (
            integrated.index + 1
        )


        integrated_output = os.path.join(
            OUTPUT_DIR,
            "38_integrated_consensus_validation.csv"
        )


        integrated.to_csv(
            integrated_output,
            index=False
        )


        print(
            "Integrated file saved:"
        )

        print(
            integrated_output
        )


# ============================================================
# 31. SAVE REPORT
# ============================================================

report_output = os.path.join(
    OUTPUT_DIR,
    "38_expression_validation_report.txt"
)


with open(
    report_output,
    "w"
) as report:

    report.write(
        "STEP 38 - CONSENSUS CANDIDATE EXPRESSION VALIDATION\n"
    )

    report.write(
        "=" * 70
        + "\n\n"
    )


    report.write(
        "Purpose:\n"
    )

    report.write(
        "Validate the six Bioinformatics-ML consensus genes "
        "using actual single-cell expression data.\n\n"
    )


    report.write(
        "Informative samples:\n"
    )

    for sample in INFORMATIVE_SAMPLES:

        report.write(
            f"- {sample}\n"
        )


    report.write(
        "\nConsensus genes:\n"
    )

    for gene in present_genes:

        report.write(
            f"- {gene}\n"
        )


    report.write(
        "\n\nRANKED EXPRESSION VALIDATION\n"
    )

    report.write(
        "-" * 70
        + "\n"
    )


    for _, row in summary_df.iterrows():

        report.write(
            f"{int(row['expression_validation_rank'])}. "
            f"{row['gene']} | "
            f"Mean log2FC="
            f"{row['mean_log2FC']:.4f} | "
            f"Positive samples="
            f"{int(row['samples_with_positive_log2FC'])}/"
            f"{int(row['samples_tested'])} | "
            f"Consistency="
            f"{row['expression_consistency']} | "
            f"Score="
            f"{row['expression_validation_score']:.4f}\n"
        )


    report.write(
        "\n\nIMPORTANT LIMITATIONS\n"
    )

    report.write(
        "-" * 70
        + "\n"
    )


    report.write(
        "1. Only three informative samples were available.\n"
    )

    report.write(
        "2. Expression was evaluated using the existing "
        "malignant and diploid_non_malignant labels.\n"
    )

    report.write(
        "3. The expression validation score is exploratory.\n"
    )

    report.write(
        "4. This analysis does not establish clinical "
        "biomarker validity.\n"
    )


# ============================================================
# 32. FINAL RESULTS
# ============================================================

print(
    "\n" + "=" * 70
)

print(
    "STEP 38 COMPLETED"
)

print(
    "=" * 70
)


print(
    "\n===== EXPRESSION VALIDATION SUMMARY ====="
)


display_columns = [

    "expression_validation_rank",

    "gene",

    "mean_log2FC",

    "mean_malignant_detection_pct",

    "mean_diploid_detection_pct",

    "mean_detection_difference_pct",

    "samples_with_positive_log2FC",

    "samples_tested",

    "expression_consistency",

    "expression_validation_score"

]


display_columns = [

    column
    for column in display_columns
    if column in summary_df.columns

]


print(
    summary_df[
        display_columns
    ].to_string(
        index=False
    )
)


print(
    "\n===== OUTPUT FILES ====="
)


print(
    "Sample-level:"
)

print(
    sample_output
)


print(
    "\nSummary:"
)

print(
    summary_output
)


if integrated_output is not None:

    print(
        "\nIntegrated Step 37 + Step 38:"
    )

    print(
        integrated_output
    )


print(
    "\nReport:"
)

print(
    report_output
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