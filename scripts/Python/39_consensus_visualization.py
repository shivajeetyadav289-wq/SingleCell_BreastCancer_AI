# ============================================================
# STEP 39
# CONSENSUS CANDIDATE VISUALIZATION & EVIDENCE INTEGRATION
# ============================================================
#
# Purpose:
#
# Integrate:
#   1. Bioinformatics evidence
#   2. ML recurrence
#   3. Single-cell expression validation
#   4. External validation where available
#
# Generate:
#   - Integrated evidence table
#   - Six-gene expression heatmap
#   - Malignant vs diploid detection plot
#   - Per-sample log2FC plot
#   - Final candidate classification
#   - Text report
#
# IMPORTANT:
# These are candidate biomarkers, not clinically validated
# biomarkers.
#
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# 1. CONFIGURATION
# ============================================================

STEP37_FILE = (
    "results/ai_ml/37_consensus_characterization/"
    "37_consensus_candidate_characterization.csv"
)

STEP38_SAMPLE_FILE = (
    "results/ai_ml/38_consensus_expression_validation/"
    "38_sample_level_expression_validation.csv"
)

STEP38_SUMMARY_FILE = (
    "results/ai_ml/38_consensus_expression_validation/"
    "38_consensus_expression_summary.csv"
)

STEP38_INTEGRATED_FILE = (
    "results/ai_ml/38_consensus_expression_validation/"
    "38_integrated_consensus_validation.csv"
)

OUTPUT_DIR = (
    "results/ai_ml/39_consensus_visualization"
)

FIGURE_DIR = os.path.join(
    OUTPUT_DIR,
    "figures"
)


# ============================================================
# 2. CREATE DIRECTORIES
# ============================================================

os.makedirs(
    OUTPUT_DIR,
    exist_ok=True
)

os.makedirs(
    FIGURE_DIR,
    exist_ok=True
)


# ============================================================
# 3. START
# ============================================================

print("=" * 70)
print("STEP 39 - CONSENSUS CANDIDATE VISUALIZATION")
print("=" * 70)


# ============================================================
# 4. CHECK INPUT FILES
# ============================================================

print(
    "\n===== INPUT FILE CHECK ====="
)


required_files = {

    "Step 37 characterization":
        STEP37_FILE,

    "Step 38 sample-level results":
        STEP38_SAMPLE_FILE,

    "Step 38 summary":
        STEP38_SUMMARY_FILE,

    "Step 38 integrated results":
        STEP38_INTEGRATED_FILE

}


for name, path in required_files.items():

    if os.path.exists(path):

        print(
            f"{name}: FOUND"
        )

    else:

        print(
            f"{name}: NOT FOUND"
        )


if not os.path.exists(
    STEP38_SUMMARY_FILE
):

    raise FileNotFoundError(
        "Step 38 summary file is required."
    )


# ============================================================
# 5. LOAD STEP 38 SUMMARY
# ============================================================

print(
    "\n===== LOADING STEP 38 SUMMARY ====="
)


summary = pd.read_csv(
    STEP38_SUMMARY_FILE
)


if "gene" not in summary.columns:

    raise ValueError(
        "Step 38 summary does not contain gene column."
    )


summary["gene"] = (
    summary["gene"]
    .astype(str)
    .str.strip()
)


print(
    "Consensus genes:",
    len(summary)
)


print(
    summary[
        "gene"
    ].tolist()
)


# ============================================================
# 6. LOAD STEP 38 SAMPLE-LEVEL RESULTS
# ============================================================

print(
    "\n===== LOADING SAMPLE-LEVEL RESULTS ====="
)


sample_results = pd.read_csv(
    STEP38_SAMPLE_FILE
)


sample_results["gene"] = (
    sample_results["gene"]
    .astype(str)
    .str.strip()
)


print(
    "Sample-level rows:",
    len(sample_results)
)


print(
    "Samples:",
    sorted(
        sample_results[
            "sample"
        ].unique()
    )
)


# ============================================================
# 7. LOAD STEP 37
# ============================================================

step37 = None


if os.path.exists(
    STEP37_FILE
):

    print(
        "\n===== LOADING STEP 37 ====="
    )


    step37 = pd.read_csv(
        STEP37_FILE
    )


    step37["gene"] = (
        step37["gene"]
        .astype(str)
        .str.strip()
    )


    print(
        "Step 37 rows:",
        len(step37)
    )


# ============================================================
# 8. CREATE INTEGRATED EVIDENCE TABLE
# ============================================================

print(
    "\n===== INTEGRATING EVIDENCE ====="
)


integrated = summary.copy()


# ------------------------------------------------------------
# Add Step 37 evidence
# ------------------------------------------------------------

if step37 is not None:

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


    available_columns = [

        column
        for column in step37_columns
        if column in step37.columns

    ]


    step37_small = (
        step37[
            available_columns
        ]
        .drop_duplicates(
            subset=["gene"]
        )
    )


    integrated = integrated.merge(
        step37_small,
        on="gene",
        how="left"
    )


# ============================================================
# 9. ADD EVIDENCE FLAGS
# ============================================================

# ------------------------------------------------------------
# ML recurrence
# ------------------------------------------------------------

if "ML_folds_selected" in integrated.columns:

    integrated[
        "ML_3_of_3"
    ] = (
        pd.to_numeric(
            integrated[
                "ML_folds_selected"
            ],
            errors="coerce"
        )
        == 3
    )

else:

    integrated[
        "ML_3_of_3"
    ] = np.nan


# ------------------------------------------------------------
# Expression consistency
# ------------------------------------------------------------

integrated[
    "expression_consistent_3_of_3"
] = (
    integrated[
        "samples_with_positive_log2FC"
    ]
    ==
    integrated[
        "samples_tested"
    ]
)


# ------------------------------------------------------------
# Strong expression support
# ------------------------------------------------------------

integrated[
    "strong_expression_support"
] = (

    (
        integrated[
            "mean_log2FC"
        ]
        >= 1
    )

    &

    (
        integrated[
            "samples_with_positive_log2FC"
        ]
        >= 2
    )

)


# ------------------------------------------------------------
# External support
# ------------------------------------------------------------

if "external_support" in integrated.columns:

    integrated[
        "external_validation_available"
    ] = (
        integrated[
            "external_support"
        ]
        .fillna("")
        .astype(str)
        .str.lower()
        .ne("")
    )

else:

    integrated[
        "external_validation_available"
    ] = False


# ============================================================
# 10. FINAL EVIDENCE CATEGORY
# ============================================================

def classify_candidate(row):

    ml_3 = (
        row.get(
            "ML_3_of_3",
            False
        )
        is True
    )


    expression_3 = (
        row.get(
            "expression_consistent_3_of_3",
            False
        )
        is True
    )


    strong_expression = (
        row.get(
            "strong_expression_support",
            False
        )
        is True
    )


    external = (
        row.get(
            "external_validation_available",
            False
        )
        is True
    )


    # --------------------------------------------------------
    # Strongest category
    # --------------------------------------------------------

    if (
        ml_3
        and expression_3
        and strong_expression
        and external
    ):

        return (
            "High-priority consensus candidate"
        )


    # --------------------------------------------------------
    # Strong consensus
    # --------------------------------------------------------

    if (
        ml_3
        and expression_3
        and strong_expression
    ):

        return (
            "Strong consensus candidate"
        )


    # --------------------------------------------------------
    # Moderate consensus
    # --------------------------------------------------------

    if (
        ml_3
        and expression_3
    ):

        return (
            "Moderate consensus candidate"
        )


    # --------------------------------------------------------
    # ML + expression support
    # --------------------------------------------------------

    if (
        ml_3
        and strong_expression
    ):

        return (
            "ML-supported expression candidate"
        )


    return (
        "Candidate requiring further validation"
    )


integrated[
    "final_evidence_category"
] = integrated.apply(
    classify_candidate,
    axis=1
)


# ============================================================
# 11. CREATE OVERALL EVIDENCE SCORE
# ============================================================
#
# This is an exploratory prioritization score.
#
# It is NOT a statistical significance measure.
#
# Components:
#
#   30% ML recurrence
#   40% expression validation
#   20% external support
#   10% biological category presence
#
# ============================================================

# ------------------------------------------------------------
# ML score
# ------------------------------------------------------------

if "ML_folds_selected" in integrated.columns:

    ml_folds = pd.to_numeric(
        integrated[
            "ML_folds_selected"
        ],
        errors="coerce"
    )

    ml_score = (
        ml_folds / 3.0
    )

else:

    ml_score = pd.Series(
        0.0,
        index=integrated.index
    )


# ------------------------------------------------------------
# Expression score
# ------------------------------------------------------------

expression_score = (
    pd.to_numeric(
        integrated[
            "expression_validation_score"
        ],
        errors="coerce"
    )
    .fillna(0)
)


# Normalize to 0-1 if necessary

if (
    expression_score.max()
    >
    1
):

    expression_score = (
        expression_score
        /
        expression_score.max()
    )


# ------------------------------------------------------------
# External score
# ------------------------------------------------------------

external_score = pd.Series(
    0.0,
    index=integrated.index
)


if "external_expression_prevalence" in integrated.columns:

    external_values = pd.to_numeric(
        integrated[
            "external_expression_prevalence"
        ],
        errors="coerce"
    )


    external_score = (
        external_values
        .fillna(0)
        .clip(
            lower=0,
            upper=1
        )
    )


# ------------------------------------------------------------
# Biological evidence
# ------------------------------------------------------------

biological_score = pd.Series(
    1.0,
    index=integrated.index
)


# ------------------------------------------------------------
# Overall score
# ------------------------------------------------------------

integrated[
    "overall_exploratory_evidence_score"
] = (

    0.30 * ml_score

    +

    0.40 * expression_score

    +

    0.20 * external_score

    +

    0.10 * biological_score

)


# ============================================================
# 12. FINAL RANK
# ============================================================

integrated = (
    integrated
    .sort_values(
        "overall_exploratory_evidence_score",
        ascending=False
    )
    .reset_index(
        drop=True
    )
)


integrated[
    "final_evidence_rank"
] = (
    integrated.index + 1
)


# ============================================================
# 13. SAVE INTEGRATED TABLE
# ============================================================

integrated_file = os.path.join(
    OUTPUT_DIR,
    "39_final_integrated_evidence_table.csv"
)


integrated.to_csv(
    integrated_file,
    index=False
)


print(
    "\nIntegrated evidence table saved:"
)

print(
    integrated_file
)


# ============================================================
# 14. DISPLAY FINAL TABLE
# ============================================================

display_columns = [

    "final_evidence_rank",

    "gene",

    "biological_category",

    "ML_folds_selected",

    "mean_log2FC",

    "samples_with_positive_log2FC",

    "samples_tested",

    "mean_malignant_detection_pct",

    "mean_diploid_detection_pct",

    "expression_consistency",

    "external_support",

    "final_evidence_category",

    "overall_exploratory_evidence_score"

]


display_columns = [

    column
    for column in display_columns
    if column in integrated.columns

]


print(
    "\n===== FINAL INTEGRATED EVIDENCE ====="
)


print(
    integrated[
        display_columns
    ].to_string(
        index=False
    )
)


# ============================================================
# 15. PLOT 1
# MALIGNANT VS DIPLOID DETECTION
# ============================================================

print(
    "\n===== GENERATING DETECTION PLOT ====="
)


plot_df = integrated.copy()


x = np.arange(
    len(plot_df)
)

width = 0.35


fig, ax = plt.subplots(
    figsize=(10, 6)
)


ax.bar(
    x - width / 2,
    plot_df[
        "mean_malignant_detection_pct"
    ],
    width,
    label="Malignant"
)


ax.bar(
    x + width / 2,
    plot_df[
        "mean_diploid_detection_pct"
    ],
    width,
    label="Diploid/non-malignant"
)


ax.set_xticks(
    x
)

ax.set_xticklabels(
    plot_df[
        "gene"
    ],
    rotation=45,
    ha="right"
)


ax.set_ylabel(
    "Mean detection (%)"
)


ax.set_title(
    "Consensus Gene Detection in Malignant vs Diploid Cells"
)


ax.legend()


plt.tight_layout()


detection_plot = os.path.join(
    FIGURE_DIR,
    "39_malignant_vs_diploid_detection.png"
)


plt.savefig(
    detection_plot,
    dpi=300,
    bbox_inches="tight"
)


plt.close()


print(
    detection_plot
)


# ============================================================
# 16. PLOT 2
# MEAN LOG2FC
# ============================================================

print(
    "\n===== GENERATING LOG2FC PLOT ====="
)


fig, ax = plt.subplots(
    figsize=(10, 6)
)


ax.bar(
    plot_df[
        "gene"
    ],
    plot_df[
        "mean_log2FC"
    ]
)


ax.axhline(
    0,
    linewidth=1
)


ax.set_ylabel(
    "Mean log2 fold-change"
)


ax.set_title(
    "Malignant vs Diploid Expression of Consensus Genes"
)


plt.xticks(
    rotation=45,
    ha="right"
)


plt.tight_layout()


log2fc_plot = os.path.join(
    FIGURE_DIR,
    "39_consensus_mean_log2FC.png"
)


plt.savefig(
    log2fc_plot,
    dpi=300,
    bbox_inches="tight"
)


plt.close()


print(
    log2fc_plot
)


# ============================================================
# 17. PLOT 3
# SAMPLE-LEVEL LOG2FC
# ============================================================

print(
    "\n===== GENERATING SAMPLE-LEVEL LOG2FC PLOT ====="
)


samples = (
    sample_results[
        "sample"
    ]
    .drop_duplicates()
    .tolist()
)


genes = (
    integrated[
        "gene"
    ]
    .tolist()
)


fig, ax = plt.subplots(
    figsize=(11, 7)
)


for gene in genes:

    gene_data = sample_results[
        sample_results[
            "gene"
        ] == gene
    ]


    values = []


    for sample in samples:

        row = gene_data[
            gene_data[
                "sample"
            ] == sample
        ]


        if len(row) == 0:

            values.append(
                np.nan
            )

        else:

            values.append(
                row[
                    "log2FC_malignant_vs_diploid"
                ].iloc[0]
            )


    ax.plot(
        samples,
        values,
        marker="o",
        label=gene
    )


ax.axhline(
    0,
    linewidth=1
)


ax.set_ylabel(
    "log2FC (malignant vs diploid)"
)


ax.set_xlabel(
    "Patient/sample"
)


ax.set_title(
    "Sample-Level Expression Consistency"
)


ax.legend(
    bbox_to_anchor=(
        1.02,
        1
    ),
    loc="upper left"
)


plt.tight_layout()


sample_plot = os.path.join(
    FIGURE_DIR,
    "39_sample_level_log2FC.png"
)


plt.savefig(
    sample_plot,
    dpi=300,
    bbox_inches="tight"
)


plt.close()


print(
    sample_plot
)


# ============================================================
# 18. PLOT 4
# EXPRESSION VALIDATION SCORE
# ============================================================

print(
    "\n===== GENERATING EVIDENCE SCORE PLOT ====="
)


fig, ax = plt.subplots(
    figsize=(10, 6)
)


ax.bar(
    plot_df[
        "gene"
    ],
    plot_df[
        "expression_validation_score"
    ]
)


ax.set_ylabel(
    "Expression validation score"
)


ax.set_ylim(
    0,
    1.05
)


ax.set_title(
    "Consensus Candidate Expression Validation"
)


plt.xticks(
    rotation=45,
    ha="right"
)


plt.tight_layout()


score_plot = os.path.join(
    FIGURE_DIR,
    "39_expression_validation_score.png"
)


plt.savefig(
    score_plot,
    dpi=300,
    bbox_inches="tight"
)


plt.close()


print(
    score_plot
)


# ============================================================
# 19. CREATE FINAL CANDIDATE PANEL
# ============================================================

final_panel_columns = [

    "final_evidence_rank",

    "gene",

    "biological_category",

    "ML_folds_selected",

    "mean_log2FC",

    "mean_malignant_detection_pct",

    "mean_diploid_detection_pct",

    "mean_detection_difference_pct",

    "samples_with_positive_log2FC",

    "samples_tested",

    "expression_consistency",

    "external_support",

    "final_evidence_category",

    "overall_exploratory_evidence_score"

]


final_panel_columns = [

    column
    for column in final_panel_columns
    if column in integrated.columns

]


final_panel = integrated[
    final_panel_columns
].copy()


final_panel_file = os.path.join(
    OUTPUT_DIR,
    "39_final_candidate_panel.csv"
)


final_panel.to_csv(
    final_panel_file,
    index=False
)


print(
    "\nFinal candidate panel saved:"
)

print(
    final_panel_file
)


# ============================================================
# 20. CREATE TEXT REPORT
# ============================================================

report_file = os.path.join(
    OUTPUT_DIR,
    "39_final_consensus_report.txt"
)


with open(
    report_file,
    "w"
) as report:

    report.write(
        "STEP 39 - FINAL CONSENSUS CANDIDATE ANALYSIS\n"
    )

    report.write(
        "=" * 70
        + "\n\n"
    )


    report.write(
        "Purpose:\n"
    )

    report.write(
        "Integrate bioinformatics prioritization, ML recurrence, "
        "single-cell expression validation, and external evidence "
        "for the six consensus candidate genes.\n\n"
    )


    report.write(
        "FINAL CANDIDATE PANEL\n"
    )

    report.write(
        "-" * 70
        + "\n"
    )


    for _, row in final_panel.iterrows():

        report.write(
            f"Rank {int(row['final_evidence_rank'])}: "
            f"{row['gene']}\n"
        )


        if (
            "biological_category"
            in row.index
        ):

            report.write(
                f"  Biological category: "
                f"{row['biological_category']}\n"
            )


        if (
            "ML_folds_selected"
            in row.index
        ):

            report.write(
                f"  ML recurrence: "
                f"{row['ML_folds_selected']}/3\n"
            )


        report.write(
            f"  Mean log2FC: "
            f"{row['mean_log2FC']:.3f}\n"
        )


        report.write(
            f"  Malignant detection: "
            f"{row['mean_malignant_detection_pct']:.2f}%\n"
        )


        report.write(
            f"  Diploid detection: "
            f"{row['mean_diploid_detection_pct']:.2f}%\n"
        )


        report.write(
            f"  Detection difference: "
            f"{row['mean_detection_difference_pct']:.2f}%\n"
        )


        report.write(
            f"  Expression consistency: "
            f"{row['expression_consistency']}\n"
        )


        report.write(
            f"  Final evidence category: "
            f"{row['final_evidence_category']}\n"
        )


        report.write(
            "\n"
        )


    report.write(
        "\nINTERPRETATION\n"
    )

    report.write(
        "-" * 70
        + "\n"
    )


    report.write(
        "The final six genes represent consensus candidates "
        "supported by multiple computational evidence layers. "
        "They should not be described as clinically validated "
        "biomarkers.\n\n"
    )


    report.write(
        "The strongest evidence comes from the combination of "
        "bioinformatics prioritization, recurrent ML feature "
        "selection, and consistent malignant-associated expression "
        "across the available informative samples.\n\n"
    )


    report.write(
        "LIMITATIONS\n"
    )

    report.write(
        "-" * 70
        + "\n"
    )


    report.write(
        "1. Only three informative patient samples were available.\n"
    )

    report.write(
        "2. LOPO performance was limited in previous ML analyses.\n"
    )

    report.write(
        "3. External validation was not available for every candidate.\n"
    )

    report.write(
        "4. The integrated evidence score is exploratory.\n"
    )

    report.write(
        "5. Independent datasets and experimental validation are "
        "required before clinical biomarker claims.\n"
    )


# ============================================================
# 21. FINAL OUTPUT
# ============================================================

print(
    "\n" + "=" * 70
)

print(
    "STEP 39 COMPLETED"
)

print(
    "=" * 70
)


print(
    "\n===== FINAL CANDIDATE PANEL ====="
)


print(
    final_panel.to_string(
        index=False
    )
)


print(
    "\n===== OUTPUT FILES ====="
)


print(
    "Integrated evidence:"
)

print(
    integrated_file
)


print(
    "\nFinal candidate panel:"
)

print(
    final_panel_file
)


print(
    "\nReport:"
)

print(
    report_file
)


print(
    "\nFigures:"
)

print(
    detection_plot
)

print(
    log2fc_plot
)

print(
    sample_plot
)

print(
    score_plot
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