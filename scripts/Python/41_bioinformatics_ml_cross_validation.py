# ============================================================
# STEP 41
# BIOINFORMATICS–ML CROSS-VALIDATION
# ============================================================
#
# Purpose:
#
# Cross-check the final consensus candidates against:
#
#   1. Original bioinformatics prioritization
#   2. ML recurrence
#   3. Step 38 expression validation
#   4. Step 39 integrated evidence
#
# This step does NOT perform a new ML model.
# It integrates existing evidence.
#
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# 1. FILE PATHS
# ============================================================

BIOINFO_FILE = (
    "results/malignant/biological_validation/"
    "23_final_biological_candidate_annotation.csv"
)

STEP36_FULL = (
    "results/ai_ml/36_bioinformatics_ml_consensus/"
    "36_full_ml_bioinformatics_integration.csv"
)

STEP36_STRONG = (
    "results/ai_ml/36_bioinformatics_ml_consensus/"
    "36_strong_consensus_candidates.csv"
)

STEP36_MODERATE = (
    "results/ai_ml/36_bioinformatics_ml_consensus/"
    "36_moderate_consensus_candidates.csv"
)

STEP36_SUMMARY = (
    "results/ai_ml/36_bioinformatics_ml_consensus/"
    "36_consensus_summary.csv"
)

STEP38_FILE = (
    "results/ai_ml/38_consensus_expression_validation/"
    "38_consensus_expression_summary.csv"
)

STEP39_FILE = (
    "results/ai_ml/39_consensus_visualization/"
    "39_final_integrated_evidence_table.csv"
)


# ============================================================
# OUTPUT DIRECTORY
# ============================================================

OUTPUT_DIR = (
    "results/ai_ml/41_bioinformatics_ml_cross_validation"
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
# FINAL SIX
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
# START
# ============================================================

print("=" * 70)

print(
    "STEP 41 - BIOINFORMATICS–ML CROSS-VALIDATION"
)

print("=" * 70)


# ============================================================
# 2. CHECK FILES
# ============================================================

print(
    "\n===== INPUT FILE CHECK ====="
)


required_files = {

    "Bioinformatics annotation":
        BIOINFO_FILE,

    "Step 36 full integration":
        STEP36_FULL,

    "Step 36 strong consensus":
        STEP36_STRONG,

    "Step 36 moderate consensus":
        STEP36_MODERATE,

    "Step 36 summary":
        STEP36_SUMMARY,

    "Step 38 expression validation":
        STEP38_FILE,

    "Step 39 integrated evidence":
        STEP39_FILE

}


for name, path in required_files.items():

    if os.path.exists(path):

        print(
            f"{name}: FOUND"
        )

    else:

        raise FileNotFoundError(
            f"{name} not found:\n{path}"
        )


# ============================================================
# 3. LOAD BIOINFORMATICS DATA
# ============================================================

print(
    "\n===== LOADING BIOINFORMATICS DATA ====="
)


bioinfo = pd.read_csv(
    BIOINFO_FILE
)


bioinfo["gene"] = (
    bioinfo["gene"]
    .astype(str)
    .str.strip()
)


print(
    "Bioinformatics rows:",
    len(bioinfo)
)


print(
    "Bioinformatics columns:",
    len(bioinfo.columns)
)


# ============================================================
# 4. LOAD STEP 36 FULL INTEGRATION
# ============================================================

print(
    "\n===== LOADING STEP 36 FULL INTEGRATION ====="
)


step36 = pd.read_csv(
    STEP36_FULL
)


step36["gene"] = (
    step36["gene"]
    .astype(str)
    .str.strip()
)


print(
    "Step 36 integration rows:",
    len(step36)
)


print(
    "Columns:"
)

for column in step36.columns:

    print(
        " -",
        column
    )


# ============================================================
# 5. LOAD STRONG CONSENSUS
# ============================================================

strong = pd.read_csv(
    STEP36_STRONG
)


strong["gene"] = (
    strong["gene"]
    .astype(str)
    .str.strip()
)


strong_genes = set(
    strong["gene"]
)


print(
    "\nStrong consensus genes:",
    len(strong_genes)
)


print(
    sorted(
        strong_genes
    )
)


# ============================================================
# 6. LOAD MODERATE CONSENSUS
# ============================================================

moderate = pd.read_csv(
    STEP36_MODERATE
)


moderate["gene"] = (
    moderate["gene"]
    .astype(str)
    .str.strip()
)


moderate_genes = set(
    moderate["gene"]
)


print(
    "\nModerate consensus genes:",
    len(moderate_genes)
)


# ============================================================
# 7. LOAD STEP 36 SUMMARY
# ============================================================

summary36 = pd.read_csv(
    STEP36_SUMMARY
)


print(
    "\n===== STEP 36 SUMMARY ====="
)

print(
    summary36.to_string(
        index=False
    )
)


# ============================================================
# 8. LOAD STEP 38
# ============================================================

print(
    "\n===== LOADING STEP 38 ====="
)


step38 = pd.read_csv(
    STEP38_FILE
)


step38["gene"] = (
    step38["gene"]
    .astype(str)
    .str.strip()
)


step38_final = (
    step38[
        step38["gene"].isin(
            FINAL_SIX
        )
    ]
    .copy()
)


print(
    "Final six found in Step 38:",
    len(step38_final)
)


# ============================================================
# 9. LOAD STEP 39
# ============================================================

print(
    "\n===== LOADING STEP 39 ====="
)


step39 = pd.read_csv(
    STEP39_FILE
)


step39["gene"] = (
    step39["gene"]
    .astype(str)
    .str.strip()
)


step39_final = (
    step39[
        step39["gene"].isin(
            FINAL_SIX
        )
    ]
    .copy()
)


print(
    "Final six found in Step 39:",
    len(step39_final)
)


# ============================================================
# 10. BIOINFORMATICS EVIDENCE
# ============================================================

print(
    "\n===== EXTRACTING BIOINFORMATICS EVIDENCE ====="
)


bioinfo_final = (
    bioinfo[
        bioinfo["gene"].isin(
            FINAL_SIX
        )
    ]
    .copy()
)


bioinfo_final = (
    bioinfo_final
    .sort_values(
        "priority_rank"
        if "priority_rank" in bioinfo_final.columns
        else "gene"
    )
    .drop_duplicates(
        subset=["gene"]
    )
)


print(
    "Final six in original annotation:",
    len(bioinfo_final)
)


# ============================================================
# 11. PREPARE CROSS-VALIDATION TABLE
# ============================================================

cross = pd.DataFrame(
    {
        "gene": FINAL_SIX
    }
)


# ============================================================
# 12. BIOINFORMATICS PRESENCE
# ============================================================

cross[
    "bioinformatics_supported"
] = (
    cross["gene"]
    .isin(
        set(
            bioinfo["gene"]
        )
    )
)


# ============================================================
# 13. FINAL 50 PRESENCE
# ============================================================

cross[
    "bioinformatics_final_50"
] = False


if "in_final_50" in step36.columns:

    temp = (
        step36[
            [
                "gene",
                "in_final_50"
            ]
        ]
        .drop_duplicates(
            subset=["gene"]
        )
    )


    temp = temp[
        temp["gene"].isin(
            FINAL_SIX
        )
    ]


    cross = cross.merge(
        temp,
        on="gene",
        how="left"
    )


    cross[
        "bioinformatics_final_50"
    ] = (
        cross[
            "in_final_50"
        ]
        .fillna(False)
        .astype(bool)
    )


    cross = cross.drop(
        columns=[
            "in_final_50"
        ]
    )


# ============================================================
# 14. ML RECURRENCE
# ============================================================

if "folds_selected" in step36.columns:

    temp = (
        step36[
            [
                "gene",
                "folds_selected"
            ]
        ]
        .drop_duplicates(
            subset=["gene"]
        )
    )


    temp = temp[
        temp["gene"].isin(
            FINAL_SIX
        )
    ]


    cross = cross.merge(
        temp,
        on="gene",
        how="left"
    )


    cross = cross.rename(
        columns={
            "folds_selected":
            "ML_folds_selected"
        }
    )


else:

    cross[
        "ML_folds_selected"
    ] = np.nan


# ============================================================
# 15. ML CONSENSUS STATUS
# ============================================================

cross[
    "ML_consensus_status"
] = np.where(

    cross[
        "gene"
    ].isin(
        strong_genes
    ),

    "Strong_consensus",

    np.where(

        cross[
            "gene"
        ].isin(
            moderate_genes
        ),

        "Moderate_consensus",

        "Other"

    )

)


# ============================================================
# 16. ADD BIOINFORMATICS SCORES
# ============================================================

bioinfo_score_columns = [

    "priority_rank",

    "biomarker_score",

    "adjusted_biomarker_score",

    "avg_log2FC",

    "prevalence_difference",

    "mean_sample_prevalence",

    "consistency_category",

    "biological_class",

    "final_biological_category"

]


available_bioinfo_columns = [

    column
    for column in bioinfo_score_columns
    if column in bioinfo_final.columns

]


bioinfo_small = bioinfo_final[
    [
        "gene"
    ]
    +
    available_bioinfo_columns
].copy()


cross = cross.merge(
    bioinfo_small,
    on="gene",
    how="left"
)


# ============================================================
# 17. ADD STEP 38 EXPRESSION EVIDENCE
# ============================================================

expression_columns = [

    "mean_log2FC",

    "mean_malignant_detection_pct",

    "mean_diploid_detection_pct",

    "mean_detection_difference_pct",

    "samples_with_positive_log2FC",

    "samples_tested",

    "expression_consistency",

    "expression_validation_score"

]


# Avoid duplicate avg_log2FC conflict
expression_columns = [

    column
    for column in expression_columns
    if column in step38_final.columns

]


expression_small = step38_final[
    [
        "gene"
    ]
    +
    expression_columns
].copy()


# Rename if a bioinformatics column has the same name

if "mean_log2FC" in expression_small.columns:

    expression_small = (
        expression_small
        .rename(
            columns={
                "mean_log2FC":
                "expression_mean_log2FC"
            }
        )
    )


cross = cross.merge(
    expression_small,
    on="gene",
    how="left"
)


# ============================================================
# 18. ADD STEP 39 EVIDENCE
# ============================================================

step39_columns = [

    "final_evidence_rank",

    "final_evidence_category",

    "overall_exploratory_evidence_score",

    "external_support"

]


available_step39_columns = [

    column
    for column in step39_columns
    if column in step39_final.columns

]


if available_step39_columns:

    step39_small = step39_final[
        [
            "gene"
        ]
        +
        available_step39_columns
    ].copy()


    cross = cross.merge(
        step39_small,
        on="gene",
        how="left"
    )


# ============================================================
# 19. CREATE AGREEMENT FLAGS
# ============================================================

cross[
    "bioinformatics_and_ML_agreement"
] = (

    cross[
        "bioinformatics_final_50"
    ].astype(bool)

    &

    (
        pd.to_numeric(
            cross[
                "ML_folds_selected"
            ],
            errors="coerce"
        )
        >= 2
    )

)


cross[
    "strong_ML_recurrence"
] = (

    pd.to_numeric(
        cross[
            "ML_folds_selected"
        ],
        errors="coerce"
    )
    == 3

)


cross[
    "expression_consistent"
] = (

    cross[
        "expression_consistency"
    ]
    .fillna("")
    .astype(str)
    .str.contains(
        "Consistent",
        case=False,
        na=False
    )

)


# ============================================================
# 20. MULTI-LAYER EVIDENCE COUNT
# ============================================================

cross[
    "evidence_layers"
] = (

    cross[
        "bioinformatics_supported"
    ].astype(int)

    +

    cross[
        "bioinformatics_final_50"
    ].astype(int)

    +

    cross[
        "strong_ML_recurrence"
    ].astype(int)

    +

    cross[
        "expression_consistent"
    ].astype(int)

)


# ============================================================
# 21. CROSS-VALIDATION CLASSIFICATION
# ============================================================

def classify(row):

    layers = row[
        "evidence_layers"
    ]


    if layers >= 4:

        return (
            "Strong cross-method agreement"
        )


    if layers == 3:

        return (
            "Moderate cross-method agreement"
        )


    if layers == 2:

        return (
            "ML/bioinformatics supported"
        )


    return (
        "Requires further validation"
    )


cross[
    "cross_validation_class"
] = cross.apply(
    classify,
    axis=1
)


# ============================================================
# 22. SORT
# ============================================================

cross = (
    cross
    .sort_values(
        [
            "evidence_layers",
            "strong_ML_recurrence",
            "adjusted_biomarker_score"
            if "adjusted_biomarker_score"
            in cross.columns
            else "gene"
        ],
        ascending=[
            False,
            False,
            False
        ]
    )
    .reset_index(
        drop=True
    )
)


cross[
    "cross_validation_rank"
] = (
    cross.index + 1
)


# ============================================================
# 23. SAVE FULL TABLE
# ============================================================

full_output = os.path.join(
    OUTPUT_DIR,
    "41_bioinformatics_ml_cross_validation.csv"
)


cross.to_csv(
    full_output,
    index=False
)


print(
    "\nFull cross-validation table saved:"
)

print(
    full_output
)


# ============================================================
# 24. FINAL SIX EVIDENCE MATRIX
# ============================================================

matrix_columns = [

    "cross_validation_rank",

    "gene",

    "bioinformatics_final_50",

    "priority_rank",

    "biomarker_score",

    "adjusted_biomarker_score",

    "biological_class",

    "final_biological_category",

    "ML_folds_selected",

    "ML_consensus_status",

    "expression_mean_log2FC",

    "mean_malignant_detection_pct",

    "mean_diploid_detection_pct",

    "mean_detection_difference_pct",

    "samples_with_positive_log2FC",

    "samples_tested",

    "expression_consistency",

    "expression_validation_score",

    "final_evidence_category",

    "overall_exploratory_evidence_score",

    "evidence_layers",

    "cross_validation_class"

]


matrix_columns = [

    column
    for column in matrix_columns
    if column in cross.columns

]


final_matrix = cross[
    matrix_columns
].copy()


matrix_file = os.path.join(
    OUTPUT_DIR,
    "41_final_six_evidence_matrix.csv"
)


final_matrix.to_csv(
    matrix_file,
    index=False
)


# ============================================================
# 25. PRINT FINAL MATRIX
# ============================================================

print(
    "\n" + "=" * 70
)

print(
    "FINAL SIX CROSS-VALIDATION MATRIX"
)

print(
    "=" * 70
)


print(
    final_matrix.to_string(
        index=False
    )
)


# ============================================================
# 26. SUMMARY COUNTS
# ============================================================

print(
    "\n===== CROSS-VALIDATION SUMMARY ====="
)


print(
    "Final six:",
    len(cross)
)


print(
    "Bioinformatics-supported:",
    int(
        cross[
            "bioinformatics_supported"
        ].sum()
    )
)


print(
    "Present in final 50:",
    int(
        cross[
            "bioinformatics_final_50"
        ].sum()
    )
)


print(
    "ML 3/3:",
    int(
        cross[
            "strong_ML_recurrence"
        ].sum()
    )
)


print(
    "Expression consistent:",
    int(
        cross[
            "expression_consistent"
        ].sum()
    )
)


print(
    "\nCross-validation classes:"
)


print(
    cross[
        "cross_validation_class"
    ].value_counts()
    .to_string()
)


# ============================================================
# 27. WRITE SUMMARY REPORT
# ============================================================

summary_file = os.path.join(
    OUTPUT_DIR,
    "41_cross_validation_summary.txt"
)


with open(
    summary_file,
    "w"
) as report:

    report.write(
        "STEP 41 - BIOINFORMATICS–ML CROSS-VALIDATION\n"
    )

    report.write(
        "=" * 70
        + "\n\n"
    )


    report.write(
        "Purpose:\n"
    )

    report.write(
        "Cross-check the final six consensus candidates "
        "against the original bioinformatics prioritization, "
        "ML recurrence, and single-cell expression validation.\n\n"
    )


    report.write(
        "FINAL SIX GENES\n"
    )

    report.write(
        "-" * 70
        + "\n"
    )


    for _, row in cross.iterrows():

        report.write(
            f"\n{row['gene']}\n"
        )


        report.write(
            f"Bioinformatics supported: "
            f"{row['bioinformatics_supported']}\n"
        )


        report.write(
            f"Present in final 50: "
            f"{row['bioinformatics_final_50']}\n"
        )


        report.write(
            f"ML folds selected: "
            f"{row['ML_folds_selected']}\n"
        )


        report.write(
            f"ML consensus status: "
            f"{row['ML_consensus_status']}\n"
        )


        if (
            "adjusted_biomarker_score"
            in row.index
        ):

            report.write(
                f"Adjusted bioinformatics biomarker score: "
                f"{row['adjusted_biomarker_score']}\n"
            )


        if (
            "expression_mean_log2FC"
            in row.index
        ):

            report.write(
                f"Expression mean log2FC: "
                f"{row['expression_mean_log2FC']}\n"
            )


        report.write(
            f"Expression consistency: "
            f"{row['expression_consistency']}\n"
        )


        report.write(
            f"Evidence layers: "
            f"{row['evidence_layers']}\n"
        )


        report.write(
            f"Cross-validation classification: "
            f"{row['cross_validation_class']}\n"
        )


    report.write(
        "\n\nPROJECT INTERPRETATION\n"
    )

    report.write(
        "-" * 70
        + "\n"
    )


    report.write(
        "The cross-validation analysis evaluates whether "
        "the final consensus genes are supported by both "
        "the original bioinformatics prioritization and "
        "the machine-learning analysis.\n\n"
    )


    report.write(
        "The analysis does not establish clinical biomarker "
        "validity. The results represent computational and "
        "single-cell evidence supporting research candidates.\n\n"
    )


    report.write(
        "The strongest candidates are those for which "
        "bioinformatics prioritization, ML recurrence, and "
        "malignant-associated expression converge.\n"
    )


# ============================================================
# 28. FIGURE 1
# BIOINFORMATICS SCORE VS ML RECURRENCE
# ============================================================

if (
    "adjusted_biomarker_score"
    in cross.columns
):

    print(
        "\n===== GENERATING BIOINFORMATICS VS ML PLOT ====="
    )


    fig, ax = plt.subplots(
        figsize=(9, 6)
    )


    x = pd.to_numeric(
        cross[
            "ML_folds_selected"
        ],
        errors="coerce"
    )


    y = pd.to_numeric(
        cross[
            "adjusted_biomarker_score"
        ],
        errors="coerce"
    )


    ax.scatter(
        x,
        y
    )


    for _, row in cross.iterrows():

        if pd.notna(
            row[
                "ML_folds_selected"
            ]
        ) and pd.notna(
            row[
                "adjusted_biomarker_score"
            ]
        ):

            ax.annotate(
                row["gene"],
                (
                    row[
                        "ML_folds_selected"
                    ],
                    row[
                        "adjusted_biomarker_score"
                    ]
                ),
                xytext=(
                    5,
                    5
                ),
                textcoords="offset points"
            )


    ax.set_xlabel(
        "ML folds selected"
    )


    ax.set_ylabel(
        "Adjusted bioinformatics biomarker score"
    )


    ax.set_title(
        "Bioinformatics Prioritization vs ML Recurrence"
    )


    ax.set_xticks(
        [
            1,
            2,
            3
        ]
    )


    plt.tight_layout()


    plot1 = os.path.join(
        FIGURE_DIR,
        "41_bioinformatics_vs_ml_score.png"
    )


    plt.savefig(
        plot1,
        dpi=300,
        bbox_inches="tight"
    )


    plt.close()


    print(
        plot1
    )


# ============================================================
# 29. FIGURE 2
# EVIDENCE HEATMAP
# ============================================================

print(
    "\n===== GENERATING EVIDENCE HEATMAP ====="
)


heatmap_columns = [

    "bioinformatics_supported",

    "bioinformatics_final_50",

    "strong_ML_recurrence",

    "expression_consistent"

]


heatmap = (
    cross[
        [
            "gene"
        ]
        +
        heatmap_columns
    ]
    .set_index(
        "gene"
    )
    .astype(int)
)


fig, ax = plt.subplots(
    figsize=(9, 5)
)


im = ax.imshow(
    heatmap.values,
    aspect="auto"
)


ax.set_xticks(
    range(
        len(
            heatmap_columns
        )
    )
)


ax.set_xticklabels(
    [
        "Bioinformatics",
        "Final 50",
        "ML 3/3",
        "Expression 3/3"
    ],
    rotation=45,
    ha="right"
)


ax.set_yticks(
    range(
        len(
            heatmap.index
        )
    )
)


ax.set_yticklabels(
    heatmap.index
)


for i in range(
    heatmap.shape[0]
):

    for j in range(
        heatmap.shape[1]
    ):

        ax.text(
            j,
            i,
            str(
                heatmap.iloc[
                    i,
                    j
                ]
            ),
            ha="center",
            va="center"
        )


ax.set_title(
    "Cross-Method Evidence for Final Six Candidates"
)


plt.tight_layout()


plot2 = os.path.join(
    FIGURE_DIR,
    "41_final_six_evidence_heatmap.png"
)


plt.savefig(
    plot2,
    dpi=300,
    bbox_inches="tight"
)


plt.close()


print(
    plot2
)


# ============================================================
# 30. COMPLETION
# ============================================================

print(
    "\n" + "=" * 70
)

print(
    "STEP 41 COMPLETED"
)

print(
    "=" * 70
)


print(
    "\nOutputs:"
)

print(
    full_output
)

print(
    matrix_file
)

print(
    summary_file
)

print(
    plot1
    if "plot1" in locals()
    else "Plot 1 not generated"
)

print(
    plot2
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