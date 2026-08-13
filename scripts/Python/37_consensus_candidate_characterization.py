# ============================================================
# STEP 37
# CONSENSUS CANDIDATE CHARACTERIZATION AND RANKING
#
# Purpose:
# Integrate the strongest Bioinformatics + ML consensus genes.
#
# Input:
#   Step 25:
#       25_final_research_candidate_panel.csv
#
#   Step 27:
#       27_external_top_candidates.csv
#       26_external_validation_priority_candidates.csv
#
#   Step 35:
#       35_gene_selection_recurrence.csv
#
#   Step 36:
#       36_strong_consensus_candidates.csv
#       36_broad_consensus_candidates.csv
#
# Main objective:
# Characterize and rank the 6 strong consensus candidates.
#
# IMPORTANT:
# This script does NOT invent missing biological or ML values.
# If a value is unavailable in the existing files, it is reported
# as NA.
#
# ============================================================

import os
import re
import pandas as pd
import numpy as np


# ============================================================
# 1. CONFIGURATION
# ============================================================

FINAL_50_FILE = (
    "results/malignant/final_candidates/"
    "25_final_research_candidate_panel.csv"
)

EXTERNAL_TOP_FILE = (
    "results/malignant/external_validation/"
    "27_external_top_candidates.csv"
)

EXTERNAL_PRIORITY_FILE = (
    "results/malignant/external_validation/"
    "26_external_validation_priority_candidates.csv"
)

ML_RECURRENCE_FILE = (
    "results/ai_ml/35_data_driven_features/"
    "35_gene_selection_recurrence.csv"
)

STRONG_CONSENSUS_FILE = (
    "results/ai_ml/36_bioinformatics_ml_consensus/"
    "36_strong_consensus_candidates.csv"
)

BROAD_CONSENSUS_FILE = (
    "results/ai_ml/36_bioinformatics_ml_consensus/"
    "36_broad_consensus_candidates.csv"
)

OUTPUT_DIR = (
    "results/ai_ml/37_consensus_characterization"
)

os.makedirs(
    OUTPUT_DIR,
    exist_ok=True
)


# ============================================================
# 2. HELPER FUNCTIONS
# ============================================================

def clean_gene_column(df, column="gene"):

    if column not in df.columns:
        return df

    df[column] = (
        df[column]
        .astype(str)
        .str.strip()
    )

    return df


def find_column(df, patterns):

    """
    Find the first column whose name matches one
    of the supplied regular-expression patterns.
    """

    for pattern in patterns:

        for column in df.columns:

            if re.search(
                pattern,
                str(column),
                flags=re.IGNORECASE
            ):

                return column

    return None


def numeric_column(df, column):

    if column is None:
        return pd.Series(
            np.nan,
            index=df.index
        )

    return pd.to_numeric(
        df[column],
        errors="coerce"
    )


def safe_normalize(series):

    """
    Min-max normalization.

    If all non-missing values are identical,
    return 1 for available values.
    """

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

    minimum = values[valid].min()
    maximum = values[valid].max()

    if maximum == minimum:

        result.loc[valid] = 1.0

    else:

        result.loc[valid] = (
            values.loc[valid] - minimum
        ) / (
            maximum - minimum
        )

    return result


def yes_no(value):

    if pd.isna(value):
        return "NA"

    if bool(value):
        return "Yes"

    return "No"


# ============================================================
# 3. START
# ============================================================

print("=" * 70)
print("STEP 37 - CONSENSUS CANDIDATE CHARACTERIZATION")
print("=" * 70)


# ============================================================
# 4. CHECK REQUIRED FILES
# ============================================================

print("\n===== INPUT FILE CHECK =====")

required_files = {

    "Final 50":
        FINAL_50_FILE,

    "ML recurrence":
        ML_RECURRENCE_FILE,

    "Strong consensus":
        STRONG_CONSENSUS_FILE,

    "Broad consensus":
        BROAD_CONSENSUS_FILE
}


for name, path in required_files.items():

    if not os.path.exists(path):

        raise FileNotFoundError(
            f"\n{name} file not found:\n{path}"
        )

    print(
        f"{name}: FOUND"
    )


external_top_exists = os.path.exists(
    EXTERNAL_TOP_FILE
)

external_priority_exists = os.path.exists(
    EXTERNAL_PRIORITY_FILE
)


print(
    "External top candidates:",
    "FOUND" if external_top_exists else "NOT FOUND"
)

print(
    "External priority candidates:",
    "FOUND" if external_priority_exists else "NOT FOUND"
)


# ============================================================
# 5. LOAD FINAL 50
# ============================================================

print(
    "\n===== LOADING FINAL 50 ====="
)

final50 = pd.read_csv(
    FINAL_50_FILE
)

final50 = clean_gene_column(
    final50
)

if "gene" not in final50.columns:

    raise ValueError(
        "Final 50 file does not contain a 'gene' column."
    )


final50 = (
    final50
    .dropna(subset=["gene"])
    .drop_duplicates(subset=["gene"])
    .copy()
)


print(
    "Final 50 genes:",
    len(final50)
)


# ============================================================
# 6. LOAD STRONG CONSENSUS
# ============================================================

print(
    "\n===== LOADING STRONG CONSENSUS ====="
)

strong = pd.read_csv(
    STRONG_CONSENSUS_FILE
)

strong = clean_gene_column(
    strong
)

if "gene" not in strong.columns:

    raise ValueError(
        "Strong consensus file does not contain "
        "a 'gene' column."
    )


strong = (
    strong
    .dropna(subset=["gene"])
    .drop_duplicates(subset=["gene"])
    .copy()
)


strong_genes = (
    strong["gene"]
    .tolist()
)


print(
    "Strong consensus genes:",
    len(strong_genes)
)


if len(strong_genes) == 0:

    raise ValueError(
        "No strong consensus genes were found."
    )


print(
    "\nStrong consensus genes:"
)

for gene in strong_genes:

    print(
        " -",
        gene
    )


# ============================================================
# 7. LOAD BROAD CONSENSUS
# ============================================================

print(
    "\n===== LOADING BROAD CONSENSUS ====="
)

broad = pd.read_csv(
    BROAD_CONSENSUS_FILE
)

broad = clean_gene_column(
    broad
)


print(
    "Broad consensus genes:",
    len(broad)
)


# ============================================================
# 8. LOAD ML RECURRENCE
# ============================================================

print(
    "\n===== LOADING ML RECURRENCE ====="
)

ml = pd.read_csv(
    ML_RECURRENCE_FILE
)

ml = clean_gene_column(
    ml
)


if "gene" not in ml.columns:

    raise ValueError(
        "ML recurrence file does not contain "
        "a 'gene' column."
    )


if "folds_selected" not in ml.columns:

    raise ValueError(
        "ML recurrence file does not contain "
        "'folds_selected'."
    )


ml["folds_selected"] = pd.to_numeric(
    ml["folds_selected"],
    errors="coerce"
)


ml = (
    ml
    .dropna(subset=["gene"])
    .drop_duplicates(subset=["gene"])
    .copy()
)


print(
    "ML genes:",
    len(ml)
)


# ============================================================
# 9. LOAD EXTERNAL VALIDATION
# ============================================================

external = None

if external_top_exists:

    print(
        "\n===== LOADING EXTERNAL VALIDATION ====="
    )

    external = pd.read_csv(
        EXTERNAL_TOP_FILE
    )

    external = clean_gene_column(
        external
    )

    print(
        "External top candidates:",
        len(external)
    )

elif external_priority_exists:

    print(
        "\n===== LOADING EXTERNAL VALIDATION PRIORITY ====="
    )

    external = pd.read_csv(
        EXTERNAL_PRIORITY_FILE
    )

    external = clean_gene_column(
        external
    )

    print(
        "External priority candidates:",
        len(external)
    )

else:

    print(
        "\nNo external candidate CSV available."
    )


# ============================================================
# 10. BUILD CONSENSUS TABLE
# ============================================================

print(
    "\n===== BUILDING CONSENSUS TABLE ====="
)

consensus = pd.DataFrame({

    "gene":
        strong_genes

})


# ============================================================
# 11. ML RECURRENCE
# ============================================================

ml_lookup = (
    ml[
        [
            "gene",
            "folds_selected"
        ]
    ]
    .drop_duplicates(
        subset=["gene"]
    )
    .set_index(
        "gene"
    )[
        "folds_selected"
    ]
    .to_dict()
)


consensus[
    "ML_folds_selected"
] = (
    consensus[
        "gene"
    ]
    .map(
        ml_lookup
    )
)


consensus[
    "ML_recurrence"
] = (
    consensus[
        "ML_folds_selected"
    ]
    .apply(
        lambda x:
        "3/3"
        if x == 3
        else (
            "2/3"
            if x == 2
            else (
                "1/3"
                if x == 1
                else "Not_recurrent"
            )
        )
    )
)


# ============================================================
# 12. MERGE FINAL 50 INFORMATION
# ============================================================

final50_for_merge = final50.copy()


# Prevent duplicate column problems

final50_for_merge = (
    final50_for_merge
    .rename(
        columns={
            column:
            f"bio_{column}"
            for column in final50_for_merge.columns
            if column != "gene"
        }
    )
)


consensus = consensus.merge(
    final50_for_merge,
    on="gene",
    how="left"
)


# ============================================================
# 13. IDENTIFY BIOLOGICAL RANK COLUMN
# ============================================================

rank_column = find_column(
    final50,
    [
        r"rank",
        r"priority",
        r"score"
    ]
)


print(
    "\nDetected biological ranking column:",
    rank_column
)


if rank_column is not None:

    consensus[
        "biological_rank_value"
    ] = pd.to_numeric(
        consensus[
            f"bio_{rank_column}"
        ],
        errors="coerce"
    )

else:

    consensus[
        "biological_rank_value"
    ] = np.nan


# ============================================================
# 14. IDENTIFY BIOLOGICAL CATEGORY
# ============================================================

category_column = find_column(
    final50,
    [
        r"category",
        r"class",
        r"type",
        r"program",
        r"biological"
    ]
)


print(
    "Detected biological category column:",
    category_column
)


if category_column is not None:

    consensus[
        "biological_category"
    ] = (
        consensus[
            f"bio_{category_column}"
        ]
        .astype(str)
    )

else:

    consensus[
        "biological_category"
    ] = "Not_available"


# ============================================================
# 15. EXTERNAL VALIDATION
# ============================================================

if external is not None:

    print(
        "\n===== INTEGRATING EXTERNAL VALIDATION ====="
    )

    external = (
        external
        .dropna(subset=["gene"])
        .drop_duplicates(
            subset=["gene"]
        )
        .copy()
    )


    external_for_merge = external.copy()


    external_for_merge = (
        external_for_merge
        .rename(
            columns={
                column:
                f"external_{column}"
                for column in external_for_merge.columns
                if column != "gene"
            }
        )
    )


    consensus = consensus.merge(
        external_for_merge,
        on="gene",
        how="left"
    )


    # --------------------------------------------------------
    # Detect external prevalence
    # --------------------------------------------------------

    prevalence_column = find_column(
        external,
        [
            r"prevalence",
            r"frequency",
            r"fraction",
            r"proportion"
        ]
    )


    print(
        "Detected external prevalence column:",
        prevalence_column
    )


    if prevalence_column is not None:

        consensus[
            "external_expression_prevalence"
        ] = pd.to_numeric(
            consensus[
                f"external_{prevalence_column}"
            ],
            errors="coerce"
        )

    else:

        consensus[
            "external_expression_prevalence"
        ] = np.nan


    # --------------------------------------------------------
    # Detect external mean expression
    # --------------------------------------------------------

    mean_expression_column = find_column(
        external,
        [
            r"mean.*expression",
            r"external_mean",
            r"mean_expression",
            r"avg.*expression"
        ]
    )


    print(
        "Detected external mean expression column:",
        mean_expression_column
    )


    if mean_expression_column is not None:

        consensus[
            "external_mean_expression"
        ] = pd.to_numeric(
            consensus[
                f"external_{mean_expression_column}"
            ],
            errors="coerce"
        )

    else:

        consensus[
            "external_mean_expression"
        ] = np.nan

else:

    consensus[
        "external_expression_prevalence"
    ] = np.nan

    consensus[
        "external_mean_expression"
    ] = np.nan


# ============================================================
# 16. EXTERNAL SUPPORT CLASSIFICATION
# ============================================================

def classify_external(
    prevalence
):

    if pd.isna(prevalence):

        return "Not_available"

    if prevalence >= 0.75:

        return "Strong_external_support"

    elif prevalence >= 0.50:

        return "Moderate_external_support"

    elif prevalence > 0:

        return "Weak_external_support"

    else:

        return "No_external_expression"


consensus[
    "external_support"
] = (
    consensus[
        "external_expression_prevalence"
    ]
    .apply(
        classify_external
    )
)


# ============================================================
# 17. ML SUPPORT SCORE
# ============================================================

consensus[
    "ML_recurrence_score"
] = (
    consensus[
        "ML_folds_selected"
    ]
    / 3.0
)


# ============================================================
# 18. EXTERNAL SUPPORT SCORE
# ============================================================

consensus[
    "external_support_score"
] = (
    consensus[
        "external_expression_prevalence"
    ]
    .clip(
        lower=0,
        upper=1
    )
)


# ============================================================
# 19. BIOLOGICAL PRIORITY SCORE
# ============================================================

if rank_column is not None:

    rank_values = pd.to_numeric(
        consensus[
            "biological_rank_value"
        ],
        errors="coerce"
    )


    # If lower rank = better rank,
    # normalize as inverse rank.

    valid = rank_values.notna()


    consensus[
        "biological_priority_score"
    ] = np.nan


    if valid.sum() > 0:

        maximum_rank = (
            rank_values[valid].max()
        )

        minimum_rank = (
            rank_values[valid].min()
        )


        if maximum_rank != minimum_rank:

            consensus.loc[
                valid,
                "biological_priority_score"
            ] = (
                maximum_rank
                - rank_values[valid]
            ) / (
                maximum_rank
                - minimum_rank
            )

        else:

            consensus.loc[
                valid,
                "biological_priority_score"
            ] = 1.0

else:

    consensus[
        "biological_priority_score"
    ] = np.nan


# ============================================================
# 20. CONSENSUS SCORE
# ============================================================

#
# IMPORTANT:
#
# This is an exploratory prioritization score,
# NOT a statistical validation score.
#
# We use:
#
#   50% ML recurrence
#   30% external expression support
#   20% biological priority
#
# If biological priority is unavailable,
# the available components are reweighted.
#

score_components = []

weights = []


for index in consensus.index:

    ml_score = consensus.loc[
        index,
        "ML_recurrence_score"
    ]

    external_score = consensus.loc[
        index,
        "external_support_score"
    ]

    biological_score = consensus.loc[
        index,
        "biological_priority_score"
    ]


    values = []
    local_weights = []


    if not pd.isna(ml_score):

        values.append(
            ml_score
        )

        local_weights.append(
            0.50
        )


    if not pd.isna(external_score):

        values.append(
            external_score
        )

        local_weights.append(
            0.30
        )


    if not pd.isna(biological_score):

        values.append(
            biological_score
        )

        local_weights.append(
            0.20
        )


    if len(values) == 0:

        score = np.nan

    else:

        score = (
            sum(
                value * weight
                for value, weight
                in zip(
                    values,
                    local_weights
                )
            )
            /
            sum(
                local_weights
            )
        )


    score_components.append(
        score
    )


consensus[
    "exploratory_consensus_score"
] = score_components


# ============================================================
# 21. CONSENSUS LEVEL
# ============================================================

consensus[
    "consensus_level"
] = (
    "Strong_consensus"
)


# ============================================================
# 22. RANK CANDIDATES
# ============================================================

consensus = (
    consensus
    .sort_values(
        [
            "exploratory_consensus_score",
            "ML_folds_selected",
            "external_expression_prevalence"
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


consensus[
    "consensus_rank"
] = (
    consensus.index + 1
)


# ============================================================
# 23. REORDER IMPORTANT COLUMNS
# ============================================================

priority_columns = [

    "consensus_rank",

    "gene",

    "exploratory_consensus_score",

    "ML_folds_selected",

    "ML_recurrence",

    "external_expression_prevalence",

    "external_mean_expression",

    "external_support",

    "biological_priority_score",

    "biological_rank_value",

    "biological_category",

    "consensus_level"

]


existing_priority_columns = [
    column
    for column in priority_columns
    if column in consensus.columns
]


remaining_columns = [
    column
    for column in consensus.columns
    if column not in existing_priority_columns
]


consensus = consensus[
    existing_priority_columns
    + remaining_columns
]


# ============================================================
# 24. SAVE MAIN CONSENSUS TABLE
# ============================================================

main_file = os.path.join(
    OUTPUT_DIR,
    "37_consensus_candidate_characterization.csv"
)


consensus.to_csv(
    main_file,
    index=False
)


# ============================================================
# 25. SAVE COMPACT TABLE
# ============================================================

compact_columns = [

    "consensus_rank",

    "gene",

    "exploratory_consensus_score",

    "ML_folds_selected",

    "ML_recurrence",

    "external_expression_prevalence",

    "external_mean_expression",

    "external_support",

    "biological_category"

]


compact_columns = [
    column
    for column in compact_columns
    if column in consensus.columns
]


compact = consensus[
    compact_columns
].copy()


compact_file = os.path.join(
    OUTPUT_DIR,
    "37_top_consensus_candidates_compact.csv"
)


compact.to_csv(
    compact_file,
    index=False
)


# ============================================================
# 26. SAVE ML SUPPORT TABLE
# ============================================================

ml_support = consensus[
    [
        column
        for column in [
            "gene",
            "ML_folds_selected",
            "ML_recurrence",
            "ML_recurrence_score"
        ]
        if column in consensus.columns
    ]
].copy()


ml_support_file = os.path.join(
    OUTPUT_DIR,
    "37_consensus_ml_support.csv"
)


ml_support.to_csv(
    ml_support_file,
    index=False
)


# ============================================================
# 27. SAVE EXTERNAL SUPPORT TABLE
# ============================================================

external_support = consensus[
    [
        column
        for column in [
            "gene",
            "external_expression_prevalence",
            "external_mean_expression",
            "external_support"
        ]
        if column in consensus.columns
    ]
].copy()


external_support_file = os.path.join(
    OUTPUT_DIR,
    "37_consensus_external_support.csv"
)


external_support.to_csv(
    external_support_file,
    index=False
)


# ============================================================
# 28. SAVE TEXT REPORT
# ============================================================

report_file = os.path.join(
    OUTPUT_DIR,
    "37_consensus_candidate_report.txt"
)


with open(
    report_file,
    "w"
) as report:

    report.write(
        "STEP 37 - CONSENSUS CANDIDATE CHARACTERIZATION\n"
    )

    report.write(
        "=" * 70
        + "\n\n"
    )


    report.write(
        "Purpose:\n"
    )

    report.write(
        "Characterize the strongest genes supported by both "
        "bioinformatics prioritization and recurrent ML "
        "feature selection.\n\n"
    )


    report.write(
        f"Strong consensus candidates: "
        f"{len(consensus)}\n\n"
    )


    report.write(
        "IMPORTANT INTERPRETATION:\n"
    )

    report.write(
        "The consensus score is exploratory and is NOT a "
        "statistical validation score or clinical biomarker score.\n\n"
    )


    report.write(
        "RANKED CONSENSUS CANDIDATES\n"
    )

    report.write(
        "-" * 70
        + "\n"
    )


    for _, row in consensus.iterrows():

        score = row[
            "exploratory_consensus_score"
        ]

        if pd.isna(score):

            score_text = "NA"

        else:

            score_text = f"{score:.4f}"


        report.write(
            f"{int(row['consensus_rank'])}. "
            f"{row['gene']} | "
            f"Score={score_text} | "
            f"ML={row['ML_recurrence']} | "
            f"External={row['external_support']}\n"
        )


    report.write(
        "\n\nSTRONG CONSENSUS GENES\n"
    )

    report.write(
        "-" * 70
        + "\n"
    )


    for gene in strong_genes:

        report.write(
            gene + "\n"
        )


    report.write(
        "\n\nMETHOD LIMITATIONS\n"
    )

    report.write(
        "-" * 70
        + "\n"
    )

    report.write(
        "1. Only three informative samples were available "
        "for LOPO validation.\n"
    )

    report.write(
        "2. The ML model showed substantial inner-CV versus "
        "outer-LOPO performance differences.\n"
    )

    report.write(
        "3. The consensus score is an exploratory ranking "
        "framework.\n"
    )

    report.write(
        "4. Consensus candidates are not clinically validated "
        "biomarkers.\n"
    )


# ============================================================
# 29. PRINT FINAL RESULTS
# ============================================================

print(
    "\n" + "=" * 70
)

print(
    "STEP 37 COMPLETED"
)

print(
    "=" * 70
)


print(
    "\nStrong consensus candidates:",
    len(consensus)
)


print(
    "\n===== RANKED CONSENSUS CANDIDATES ====="
)


display_columns = [
    column
    for column in [
        "consensus_rank",
        "gene",
        "exploratory_consensus_score",
        "ML_recurrence",
        "external_expression_prevalence",
        "external_support"
    ]
    if column in consensus.columns
]


print(
    consensus[
        display_columns
    ].to_string(
        index=False
    )
)


print(
    "\n===== OUTPUT FILES ====="
)


print(
    main_file
)

print(
    compact_file
)

print(
    ml_support_file
)

print(
    external_support_file
)

print(
    report_file
)


print(
    "\nResults saved to:"
)

print(
    OUTPUT_DIR
)

print(
    "\n" + "=" * 70
)