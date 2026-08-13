# ============================================================
# STEP 36
# BIOINFORMATICS–ML CONSENSUS FEATURE ANALYSIS
#
# Purpose:
# Integrate the biological prioritization pipeline with
# data-driven ML feature selection.
#
# Biological feature set:
#   Step 25 final 50 candidates
#
# ML feature set:
#   Step 35 recurrent data-driven genes
#
# ML recurrence:
#   3/3 = selected in all three LOPO folds
#   2/3 = selected in at least two LOPO folds
#   1/3 = selected in one LOPO fold
#
# Main outputs:
#   - Strong consensus candidates
#   - Broad consensus candidates
#   - Final 50 with ML support
#   - ML recurrent genes outside final 50
#   - Full integration table
#   - Summary report
#
# ============================================================


import os
import pandas as pd


# ============================================================
# 1. CONFIGURATION
# ============================================================

FINAL_50_FILE = (
    "results/malignant/final_candidates/"
    "25_final_research_candidate_panel.csv"
)

ML_RECURRENCE_FILE = (
    "results/ai_ml/35_data_driven_features/"
    "35_gene_selection_recurrence.csv"
)

OUTPUT_DIR = (
    "results/ai_ml/36_bioinformatics_ml_consensus"
)


os.makedirs(
    OUTPUT_DIR,
    exist_ok=True
)


# ============================================================
# 2. START
# ============================================================

print("=" * 70)
print("STEP 36 - BIOINFORMATICS–ML CONSENSUS FEATURE ANALYSIS")
print("=" * 70)


# ============================================================
# 3. CHECK INPUT FILES
# ============================================================

print("\n===== INPUT FILE CHECK =====")


if not os.path.exists(FINAL_50_FILE):

    raise FileNotFoundError(
        "\nFinal 50 candidate file not found:\n"
        f"{FINAL_50_FILE}\n"
    )


if not os.path.exists(ML_RECURRENCE_FILE):

    raise FileNotFoundError(
        "\nML recurrence file not found:\n"
        f"{ML_RECURRENCE_FILE}\n"
    )


print(
    "Final 50 file: FOUND"
)

print(
    "ML recurrence file: FOUND"
)


# ============================================================
# 4. LOAD FINAL 50 BIOINFORMATICS CANDIDATES
# ============================================================

print(
    "\n===== LOADING FINAL BIOINFORMATICS CANDIDATES ====="
)


final50 = pd.read_csv(
    FINAL_50_FILE
)


print(
    "Rows:",
    len(final50)
)

print(
    "Columns:",
    len(final50.columns)
)


# ------------------------------------------------------------
# Check gene column
# ------------------------------------------------------------

if "gene" not in final50.columns:

    raise ValueError(
        "\nThe final candidate file does not contain "
        "a 'gene' column.\n"
        "Available columns:\n"
        + "\n".join(
            final50.columns.astype(str)
        )
    )


# ------------------------------------------------------------
# Clean gene names
# ------------------------------------------------------------

final50["gene"] = (
    final50["gene"]
    .astype(str)
    .str.strip()
)


# Remove blank / NA-like values

final50 = final50[
    final50["gene"].notna()
]

final50 = final50[
    final50["gene"] != ""
]

final50 = final50[
    final50["gene"].str.lower() != "nan"
]


# ------------------------------------------------------------
# Remove duplicate genes
# ------------------------------------------------------------

final50 = (
    final50
    .drop_duplicates(
        subset=["gene"]
    )
    .copy()
)


final50_genes = (
    final50["gene"]
    .tolist()
)


final50_set = set(
    final50_genes
)


print(
    "Unique final biological candidates:",
    len(final50_genes)
)


if len(final50_genes) != 50:

    print(
        "\nWARNING:"
    )

    print(
        "Expected 50 final candidates."
    )

    print(
        "Found:",
        len(final50_genes)
    )


# ============================================================
# 5. LOAD ML RECURRENCE RESULTS
# ============================================================

print(
    "\n===== LOADING ML RECURRENCE RESULTS ====="
)


ml = pd.read_csv(
    ML_RECURRENCE_FILE
)


print(
    "Rows:",
    len(ml)
)

print(
    "Columns:",
    len(ml.columns)
)


# ------------------------------------------------------------
# Check required columns
# ------------------------------------------------------------

required_columns = [
    "gene",
    "folds_selected"
]


for column in required_columns:

    if column not in ml.columns:

        raise ValueError(
            f"\nRequired column '{column}' "
            "not found in ML recurrence file.\n"
            "Available columns:\n"
            + "\n".join(
                ml.columns.astype(str)
            )
        )


# ------------------------------------------------------------
# Clean ML data
# ------------------------------------------------------------

ml["gene"] = (
    ml["gene"]
    .astype(str)
    .str.strip()
)


ml["folds_selected"] = pd.to_numeric(
    ml["folds_selected"],
    errors="coerce"
)


ml = ml[
    ml["gene"].notna()
]


ml = ml[
    ml["gene"] != ""
]


ml = ml[
    ml["gene"].str.lower() != "nan"
]


ml = ml[
    ml["folds_selected"].notna()
]


# ------------------------------------------------------------
# Remove duplicate genes
# ------------------------------------------------------------

ml = (
    ml
    .drop_duplicates(
        subset=["gene"]
    )
    .copy()
)


print(
    "Unique ML genes:",
    len(ml)
)


# ============================================================
# 6. CHECK RECURRENCE VALUES
# ============================================================

print(
    "\n===== ML RECURRENCE DISTRIBUTION ====="
)


recurrence_distribution = (
    ml["folds_selected"]
    .value_counts()
    .sort_index(
        ascending=False
    )
)


print(
    recurrence_distribution
)


# ============================================================
# 7. CREATE ML FEATURE GROUPS
# ============================================================

ml_3of3 = set(
    ml.loc[
        ml["folds_selected"] == 3,
        "gene"
    ]
)


ml_2of3 = set(
    ml.loc[
        ml["folds_selected"] >= 2,
        "gene"
    ]
)


ml_1of3 = set(
    ml.loc[
        ml["folds_selected"] == 1,
        "gene"
    ]
)


print(
    "\n===== ML FEATURE GROUPS ====="
)


print(
    "ML selected in 3/3 folds:",
    len(ml_3of3)
)


print(
    "ML selected in >=2/3 folds:",
    len(ml_2of3)
)


print(
    "ML selected in 1/3 folds:",
    len(ml_1of3)
)


# ============================================================
# 8. BIOINFORMATICS × ML OVERLAPS
# ============================================================

strong_consensus = (
    final50_set
    & ml_3of3
)


moderate_consensus = (
    final50_set
    & (
        ml_2of3
        - ml_3of3
    )
)


broad_consensus = (
    final50_set
    & ml_2of3
)


ml_3of3_only = (
    ml_3of3
    - final50_set
)


ml_2of3_only = (
    ml_2of3
    - final50_set
)


final50_not_ml_3of3 = (
    final50_set
    - ml_3of3
)


final50_not_ml_2of3 = (
    final50_set
    - ml_2of3
)


# ============================================================
# 9. PRINT OVERLAP SUMMARY
# ============================================================

print(
    "\n" + "=" * 70
)

print(
    "BIOINFORMATICS–ML OVERLAP"
)

print(
    "=" * 70
)


print(
    "\nFinal 50:",
    len(final50_set)
)


print(
    "ML 3/3:",
    len(ml_3of3)
)


print(
    "ML >=2/3:",
    len(ml_2of3)
)


print(
    "ML 1/3:",
    len(ml_1of3)
)


print(
    "\nFinal 50 ∩ ML 3/3:",
    len(strong_consensus)
)


print(
    "Final 50 ∩ ML 2/3:",
    len(moderate_consensus)
)


print(
    "Final 50 ∩ ML >=2/3:",
    len(broad_consensus)
)


# ============================================================
# 10. STRONG CONSENSUS GENES
# ============================================================

print(
    "\n===== STRONG CONSENSUS ====="
)

print(
    "Final 50 + ML selected in all 3 folds"
)


if len(strong_consensus) == 0:

    print(
        "No strong consensus genes found."
    )

else:

    for gene in sorted(
        strong_consensus
    ):

        print(
            gene
        )


# ============================================================
# 11. MODERATE CONSENSUS GENES
# ============================================================

print(
    "\n===== MODERATE CONSENSUS ====="
)

print(
    "Final 50 + ML selected in exactly 2/3 folds"
)


if len(moderate_consensus) == 0:

    print(
        "No moderate consensus genes found."
    )

else:

    for gene in sorted(
        moderate_consensus
    ):

        print(
            gene
        )


# ============================================================
# 12. BROAD CONSENSUS GENES
# ============================================================

print(
    "\n===== BROAD CONSENSUS ====="
)

print(
    "Final 50 + ML selected in >=2/3 folds"
)


if len(broad_consensus) == 0:

    print(
        "No broad consensus genes found."
    )

else:

    for gene in sorted(
        broad_consensus
    ):

        print(
            gene
        )


# ============================================================
# 13. CREATE ML LOOKUP TABLE
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


# ============================================================
# 14. FINAL 50 WITH ML SUPPORT
# ============================================================

final50_support = final50.copy()


final50_support[
    "ML_folds_selected"
] = (
    final50_support[
        "gene"
    ]
    .map(
        ml_lookup
    )
    .fillna(0)
    .astype(int)
)


# ------------------------------------------------------------
# ML support category
# ------------------------------------------------------------

def classify_ml_support(
    folds
):

    if folds == 3:

        return "Strong_ML_support"

    elif folds == 2:

        return "Moderate_ML_support"

    elif folds == 1:

        return "Weak_ML_support"

    else:

        return "No_ML_recurrence"


final50_support[
    "ML_support"
] = (
    final50_support[
        "ML_folds_selected"
    ]
    .apply(
        classify_ml_support
    )
)


# ------------------------------------------------------------
# Consensus status
# ------------------------------------------------------------

def classify_consensus(
    folds
):

    if folds == 3:

        return "Strong_consensus"

    elif folds == 2:

        return "Moderate_consensus"

    elif folds == 1:

        return "Weak_ML_support"

    else:

        return "Bioinformatics_only"


final50_support[
    "consensus_status"
] = (
    final50_support[
        "ML_folds_selected"
    ]
    .apply(
        classify_consensus
    )
)


# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

final50_support_file = os.path.join(
    OUTPUT_DIR,
    "36_final50_with_ml_support.csv"
)


final50_support.to_csv(
    final50_support_file,
    index=False
)


# ============================================================
# 15. STRONG CONSENSUS TABLE
# ============================================================

strong_rows = []


for gene in sorted(
    strong_consensus
):

    strong_rows.append({

        "gene":
            gene,

        "ML_folds_selected":
            3,

        "bioinformatics_final_50":
            True,

        "consensus_level":
            "Strong_consensus"

    })


strong_df = pd.DataFrame(
    strong_rows
)


strong_file = os.path.join(
    OUTPUT_DIR,
    "36_strong_consensus_candidates.csv"
)


strong_df.to_csv(
    strong_file,
    index=False
)


# ============================================================
# 16. MODERATE CONSENSUS TABLE
# ============================================================

moderate_rows = []


for gene in sorted(
    moderate_consensus
):

    moderate_rows.append({

        "gene":
            gene,

        "ML_folds_selected":
            2,

        "bioinformatics_final_50":
            True,

        "consensus_level":
            "Moderate_consensus"

    })


moderate_df = pd.DataFrame(
    moderate_rows
)


moderate_file = os.path.join(
    OUTPUT_DIR,
    "36_moderate_consensus_candidates.csv"
)


moderate_df.to_csv(
    moderate_file,
    index=False
)


# ============================================================
# 17. BROAD CONSENSUS TABLE
# ============================================================

broad_rows = []


for gene in sorted(
    broad_consensus
):

    recurrence = int(
        ml_lookup.get(
            gene,
            0
        )
    )


    if recurrence == 3:

        level = (
            "Strong_consensus"
        )

    elif recurrence == 2:

        level = (
            "Moderate_consensus"
        )

    else:

        level = (
            "Other"
        )


    broad_rows.append({

        "gene":
            gene,

        "ML_folds_selected":
            recurrence,

        "bioinformatics_final_50":
            True,

        "consensus_level":
            level

    })


broad_df = pd.DataFrame(
    broad_rows
)


broad_df = (
    broad_df
    .sort_values(
        [
            "ML_folds_selected",
            "gene"
        ],
        ascending=[
            False,
            True
        ]
    )
)


broad_file = os.path.join(
    OUTPUT_DIR,
    "36_broad_consensus_candidates.csv"
)


broad_df.to_csv(
    broad_file,
    index=False
)


# ============================================================
# 18. ML RECURRENT GENES NOT IN FINAL 50
# ============================================================

ml_only = ml[
    ~ml["gene"].isin(
        final50_set
    )
].copy()


ml_only[
    "ML_recurrence_group"
] = (
    ml_only[
        "folds_selected"
    ]
    .apply(
        lambda x:
        "3_of_3"
        if x == 3
        else (
            "2_of_3"
            if x == 2
            else "1_of_3"
        )
    )
)


ml_only = (
    ml_only
    .sort_values(
        [
            "folds_selected",
            "gene"
        ],
        ascending=[
            False,
            True
        ]
    )
)


ml_only_file = os.path.join(
    OUTPUT_DIR,
    "36_ml_recurrent_genes_not_in_final50.csv"
)


ml_only.to_csv(
    ml_only_file,
    index=False
)


# ============================================================
# 19. FULL ML–BIOINFORMATICS INTEGRATION TABLE
# ============================================================

integration = ml.copy()


integration[
    "in_final_50"
] = (
    integration[
        "gene"
    ]
    .isin(
        final50_set
    )
)


def integration_status(
    row
):

    folds = row[
        "folds_selected"
    ]

    in_final = row[
        "in_final_50"
    ]


    if (
        in_final
        and folds == 3
    ):

        return (
            "Strong_consensus"
        )


    if (
        in_final
        and folds == 2
    ):

        return (
            "Moderate_consensus"
        )


    if (
        in_final
        and folds == 1
    ):

        return (
            "Weak_ML_support"
        )


    if in_final:

        return (
            "Bioinformatics_only"
        )


    if folds == 3:

        return (
            "ML_3of3_only"
        )


    if folds == 2:

        return (
            "ML_2of3_only"
        )


    return (
        "ML_1of3_only"
    )


integration[
    "consensus_status"
] = (
    integration
    .apply(
        integration_status,
        axis=1
    )
)


integration = (
    integration
    .sort_values(
        [
            "folds_selected",
            "in_final_50",
            "gene"
        ],
        ascending=[
            False,
            False,
            True
        ]
    )
)


integration_file = os.path.join(
    OUTPUT_DIR,
    "36_full_ml_bioinformatics_integration.csv"
)


integration.to_csv(
    integration_file,
    index=False
)


# ============================================================
# 20. FINAL 50 OVERLAP SUMMARY
# ============================================================

final50_summary = pd.DataFrame({

    "category": [

        "Final_50",

        "Final50_with_ML_3of3",

        "Final50_with_ML_2of3",

        "Final50_with_ML_1of3",

        "Final50_without_ML_recurrence"

    ],

    "gene_count": [

        len(final50_set),

        len(
            final50_set
            & ml_3of3
        ),

        len(
            final50_set
            & (
                ml_2of3
                - ml_3of3
            )
        ),

        len(
            final50_set
            & ml_1of3
        ),

        len(
            final50_set
            - ml_2of3
        )

    ]

})


final50_summary_file = os.path.join(
    OUTPUT_DIR,
    "36_final50_ml_support_summary.csv"
)


final50_summary.to_csv(
    final50_summary_file,
    index=False
)


# ============================================================
# 21. GLOBAL CONSENSUS SUMMARY
# ============================================================

summary = pd.DataFrame({

    "category": [

        "Final_50",

        "ML_3_of_3",

        "ML_2_or_more_of_3",

        "ML_1_of_3",

        "Final50_and_ML_3_of_3",

        "Final50_and_ML_exactly_2_of_3",

        "Final50_and_ML_1_of_3",

        "ML_3_of_3_only",

        "ML_2_or_more_of_3_only"

    ],

    "gene_count": [

        len(final50_set),

        len(ml_3of3),

        len(ml_2of3),

        len(ml_1of3),

        len(strong_consensus),

        len(moderate_consensus),

        len(
            final50_set
            & ml_1of3
        ),

        len(ml_3of3_only),

        len(ml_2of3_only)

    ]

})


summary_file = os.path.join(
    OUTPUT_DIR,
    "36_consensus_summary.csv"
)


summary.to_csv(
    summary_file,
    index=False
)


# ============================================================
# 22. TEXT REPORT
# ============================================================

report_file = os.path.join(
    OUTPUT_DIR,
    "36_consensus_report.txt"
)


with open(
    report_file,
    "w"
) as report:

    report.write(
        "STEP 36 - "
        "BIOINFORMATICS–ML CONSENSUS ANALYSIS\n"
    )

    report.write(
        "=" * 70
        + "\n\n"
    )


    report.write(
        "BIOINFORMATICS FEATURE SET\n"
    )

    report.write(
        f"Final biological candidates: "
        f"{len(final50_set)}\n\n"
    )


    report.write(
        "ML FEATURE RECURRENCE\n"
    )

    report.write(
        f"ML genes selected in 3/3 folds: "
        f"{len(ml_3of3)}\n"
    )

    report.write(
        f"ML genes selected in >=2/3 folds: "
        f"{len(ml_2of3)}\n"
    )

    report.write(
        f"ML genes selected in 1/3 folds: "
        f"{len(ml_1of3)}\n\n"
    )


    report.write(
        "CONSENSUS\n"
    )

    report.write(
        f"Final50 ∩ ML 3/3: "
        f"{len(strong_consensus)}\n"
    )

    report.write(
        f"Final50 ∩ ML exactly 2/3: "
        f"{len(moderate_consensus)}\n"
    )

    report.write(
        f"Final50 ∩ ML >=2/3: "
        f"{len(broad_consensus)}\n\n"
    )


    report.write(
        "STRONG CONSENSUS GENES\n"
    )

    report.write(
        "-" * 50
        + "\n"
    )


    if len(strong_consensus) == 0:

        report.write(
            "None\n"
        )

    else:

        for gene in sorted(
            strong_consensus
        ):

            report.write(
                gene + "\n"
            )


    report.write(
        "\nMODERATE CONSENSUS GENES\n"
    )

    report.write(
        "-" * 50
        + "\n"
    )


    if len(moderate_consensus) == 0:

        report.write(
            "None\n"
        )

    else:

        for gene in sorted(
            moderate_consensus
        ):

            report.write(
                gene + "\n"
            )


    report.write(
        "\nFINAL 50 NOT SELECTED BY ML >=2/3\n"
    )

    report.write(
        "-" * 50
        + "\n"
    )


    for gene in sorted(
        final50_not_ml_2of3
    ):

        report.write(
            gene + "\n"
        )


    report.write(
        "\nML 3/3 GENES NOT IN FINAL 50\n"
    )

    report.write(
        "-" * 50
        + "\n"
    )


    for gene in sorted(
        ml_3of3_only
    ):

        report.write(
            gene + "\n"
        )


# ============================================================
# 23. OUTPUT FILE LIST
# ============================================================

print(
    "\n===== OUTPUT FILES ====="
)


output_files = [

    strong_file,

    moderate_file,

    broad_file,

    final50_support_file,

    ml_only_file,

    integration_file,

    final50_summary_file,

    summary_file,

    report_file

]


for file_path in output_files:

    print(
        file_path
    )


# ============================================================
# 24. FINAL SUMMARY
# ============================================================

print(
    "\n" + "=" * 70
)

print(
    "STEP 36 COMPLETED"
)

print(
    "=" * 70
)


print(
    "\nFinal 50:",
    len(final50_set)
)


print(
    "ML 3/3:",
    len(ml_3of3)
)


print(
    "ML >=2/3:",
    len(ml_2of3)
)


print(
    "Strong consensus:",
    len(strong_consensus)
)


print(
    "Moderate consensus:",
    len(moderate_consensus)
)


print(
    "Broad consensus:",
    len(broad_consensus)
)


print(
    "\nResults saved to:"
)

print(
    OUTPUT_DIR
)