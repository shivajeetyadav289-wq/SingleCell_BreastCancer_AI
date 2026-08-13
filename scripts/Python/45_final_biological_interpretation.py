# ============================================================
# STEP 45
# FINAL BIOLOGICAL INTERPRETATION AND CANDIDATE PRIORITIZATION
# ============================================================
#
# Purpose:
#
# Convert the integrated computational evidence from Step 44
# into a biologically interpretable final candidate panel.
#
# Evidence considered:
#
#   1. Bioinformatics support
#   2. ML recurrence
#   3. Discovery expression validation
#   4. External validation using GSE176078
#   5. Tumor association
#
# IMPORTANT:
#
# This is an exploratory biological interpretation step.
# It does NOT establish clinical biomarker validity.
#
# No new ML model is trained.
# No new feature selection is performed.
# No new external dataset is downloaded.
#
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# 1. INPUT
# ============================================================

INPUT_FILE = (
    "results/ai_ml/44_final_multilayer_evidence/"
    "44_final_candidate_evidence_matrix.csv"
)


# ============================================================
# 2. OUTPUT DIRECTORY
# ============================================================

OUTPUT_DIR = (
    "results/ai_ml/45_final_biological_interpretation"
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
# 3. START
# ============================================================

print("=" * 70)
print("STEP 45 - FINAL BIOLOGICAL INTERPRETATION")
print("=" * 70)


# ============================================================
# 4. INPUT CHECK
# ============================================================

print("\n===== INPUT FILE CHECK =====")

if not os.path.exists(INPUT_FILE):

    raise FileNotFoundError(
        f"\nStep 44 output not found:\n{INPUT_FILE}\n\n"
        "Check the Step 44 output directory before running Step 45."
    )

print(
    "Step 44 evidence matrix: FOUND"
)


# ============================================================
# 5. LOAD DATA
# ============================================================

df = pd.read_csv(
    INPUT_FILE
)

print(
    "\nRows:",
    len(df)
)

print(
    "Columns:",
    len(df.columns)
)

print(
    "\nColumns:"
)

for column in df.columns:

    print(
        " -",
        column
    )


# ============================================================
# 6. REQUIRED COLUMN CHECK
# ============================================================

required_columns = [
    "gene",
    "final_multilayer_score",
    "supporting_evidence_layers"
]

missing = [
    column
    for column in required_columns
    if column not in df.columns
]

if missing:

    raise ValueError(
        "\nMissing required columns:\n"
        + "\n".join(missing)
    )


# ============================================================
# 7. CLEAN GENE NAMES
# ============================================================

df["gene"] = (
    df["gene"]
    .astype(str)
    .str.strip()
)


# ============================================================
# 8. NUMERIC CONVERSION
# ============================================================

numeric_columns = [
    "final_multilayer_score",
    "supporting_evidence_layers",
    "ML_folds_selected",
    "mean_log2FC",
    "malignant_detection_pct",
    "diploid_detection_pct",
    "malignant_detection_difference_pct",
    "malignant_vs_diploid_log2FC",
    "sample_consistency",
    "external_validation_score"
]

for column in numeric_columns:

    if column in df.columns:

        df[column] = pd.to_numeric(
            df[column],
            errors="coerce"
        )


# ============================================================
# 9. DISPLAY INPUT
# ============================================================

print(
    "\n===== STEP 44 EVIDENCE MATRIX ====="
)

display_columns = [
    "gene",
    "final_multilayer_score",
    "supporting_evidence_layers",
    "ML_folds_selected",
    "malignant_detection_difference_pct",
    "sample_consistency"
]

display_columns = [
    column
    for column in display_columns
    if column in df.columns
]

print(
    df[
        display_columns
    ].to_string(
        index=False
    )
)


# ============================================================
# 10. BIOLOGICAL INTERPRETATION
# ============================================================

def interpret_gene(row):

    gene = row["gene"]

    score = row[
        "final_multilayer_score"
    ]

    layers = row[
        "supporting_evidence_layers"
    ]

    tumor_class = str(
        row.get(
            "tumor_specificity_class",
            ""
        )
    )

    external = str(
        row.get(
            "external_support_category",
            ""
        )
    )

    expression_consistency = str(
        row.get(
            "expression_consistency",
            ""
        )
    )

    ml_folds = row.get(
        "ML_folds_selected",
        np.nan
    )

    tumor_difference = row.get(
        "malignant_detection_difference_pct",
        np.nan
    )

    sample_consistency = row.get(
        "sample_consistency",
        np.nan
    )


    # --------------------------------------------------------
    # HIGHEST PRIORITY
    # --------------------------------------------------------

    if (
        score >= 80
        and
        layers >= 4
        and
        (
            pd.isna(ml_folds)
            or
            ml_folds >= 3
        )
        and
        (
            pd.isna(sample_consistency)
            or
            sample_consistency >= 0.67
        )
    ):

        return (
            "Tier_1_High_priority_candidate"
        )


    # --------------------------------------------------------
    # STRONG
    # --------------------------------------------------------

    if (
        score >= 60
        and
        layers >= 3
    ):

        return (
            "Tier_2_Strong_candidate"
        )


    # --------------------------------------------------------
    # MODERATE
    # --------------------------------------------------------

    if (
        score >= 40
        and
        layers >= 2
    ):

        return (
            "Tier_3_Candidate_requiring_validation"
        )


    return (
        "Tier_3_Candidate_requiring_validation"
    )


df[
    "final_biological_tier"
] = df.apply(
    interpret_gene,
    axis=1
)


# ============================================================
# 11. BIOLOGICAL ROLE INTERPRETATION
# ============================================================

biological_roles = {

    "PRSS23":
        "Epithelial/tumor-associated candidate; "
        "supported by recurrent ML selection and "
        "consistent malignant enrichment.",

    "FCGRT":
        "Epithelial/tumor-associated candidate; "
        "shows consistent malignant-cell enrichment "
        "across informative samples.",

    "ETFB":
        "Mitochondrial/metabolic candidate; "
        "shows recurrent ML support and consistent "
        "malignant enrichment.",

    "MRPS7":
        "Mitochondrial/ribosomal candidate; "
        "shows strong malignant enrichment and "
        "consistent sample-level expression.",

    "EPN3":
        "Epithelial-associated candidate; "
        "shows strong malignant-vs-diploid expression "
        "difference but requires additional validation "
        "because of weaker consistency.",

    "COX8A":
        "Mitochondrial/OXPHOS-associated candidate; "
        "highly expressed in malignant cells but also "
        "highly prevalent in diploid cells, indicating "
        "lower tumor specificity."
}


df[
    "biological_interpretation"
] = df["gene"].map(
    biological_roles
)


# ============================================================
# 12. FINAL RESEARCH INTERPRETATION
# ============================================================

def research_interpretation(row):

    gene = row["gene"]

    tumor_class = str(
        row.get(
            "tumor_specificity_class",
            ""
        )
    )

    score = row[
        "final_multilayer_score"
    ]


    if gene in [
        "PRSS23",
        "FCGRT",
        "ETFB",
        "MRPS7"
    ]:

        return (
            "Strong computationally supported "
            "research candidate"
        )


    if gene == "EPN3":

        return (
            "Promising candidate requiring "
            "additional validation"
        )


    if gene == "COX8A":

        return (
            "Strongly supported expression candidate "
            "but lower tumor specificity"
        )


    if score >= 60:

        return (
            "Strong computational research candidate"
        )


    return (
        "Requires further validation"
    )


df[
    "final_research_interpretation"
] = df.apply(
    research_interpretation,
    axis=1
)


# ============================================================
# 13. FINAL PRIORITY RANK
# ============================================================

df = df.sort_values(
    [
        "final_multilayer_score",
        "supporting_evidence_layers"
    ],
    ascending=[
        False,
        False
    ]
).reset_index(
    drop=True
)


df[
    "final_research_rank"
] = (
    df.index + 1
)


# ============================================================
# 14. EVIDENCE CONVERGENCE DESCRIPTION
# ============================================================

def evidence_description(row):

    evidence = []

    if bool(
        row.get(
            "bioinformatics_support",
            False
        )
    ):

        evidence.append(
            "bioinformatics"
        )


    if bool(
        row.get(
            "ML_3_of_3",
            False
        )
    ):

        evidence.append(
            "ML_recurrence"
        )


    if bool(
        row.get(
            "discovery_expression_support",
            False
        )
    ):

        evidence.append(
            "discovery_expression"
        )


    if bool(
        row.get(
            "external_support_flag",
            False
        )
    ):

        evidence.append(
            "external_validation"
        )


    if bool(
        row.get(
            "tumor_association_flag",
            False
        )
    ):

        evidence.append(
            "tumor_association"
        )


    if not evidence:

        return "Limited evidence"


    return (
        " + ".join(
            evidence
        )
    )


df[
    "evidence_convergence"
] = df.apply(
    evidence_description,
    axis=1
)


# ============================================================
# 15. SAVE COMPLETE FINAL PANEL
# ============================================================

complete_output = os.path.join(
    OUTPUT_DIR,
    "45_final_biological_candidate_panel.csv"
)

df.to_csv(
    complete_output,
    index=False
)


# ============================================================
# 16. CREATE FINAL REPORT TABLE
# ============================================================

report_columns = [
    "final_research_rank",
    "gene",
    "final_multilayer_score",
    "supporting_evidence_layers",
    "ML_folds_selected",
    "mean_log2FC",
    "malignant_detection_pct",
    "diploid_detection_pct",
    "malignant_detection_difference_pct",
    "sample_consistency",
    "tumor_specificity_class",
    "external_support_category",
    "final_biological_tier",
    "biological_interpretation",
    "final_research_interpretation",
    "evidence_convergence"
]

report_columns = [
    column
    for column in report_columns
    if column in df.columns
]


final_panel = df[
    report_columns
].copy()


panel_output = os.path.join(
    OUTPUT_DIR,
    "45_final_research_candidate_panel.csv"
)

final_panel.to_csv(
    panel_output,
    index=False
)


# ============================================================
# 17. PRINT FINAL PANEL
# ============================================================

print(
    "\n" + "=" * 70
)

print(
    "STEP 45 - FINAL RESEARCH CANDIDATE PANEL"
)

print(
    "=" * 70
)

print(
    final_panel.to_string(
        index=False
    )
)


# ============================================================
# 18. TIER SUMMARY
# ============================================================

print(
    "\n===== FINAL TIER DISTRIBUTION ====="
)

print(
    final_panel[
        "final_biological_tier"
    ]
    .value_counts()
    .to_string()
)


# ============================================================
# 19. FINAL CANDIDATE SUMMARY
# ============================================================

print(
    "\n===== FINAL BIOLOGICAL INTERPRETATION ====="
)

for _, row in final_panel.iterrows():

    print(
        f"\n{int(row['final_research_rank'])}. "
        f"{row['gene']}"
    )

    print(
        "Score:",
        row[
            "final_multilayer_score"
        ]
    )

    print(
        "Evidence layers:",
        row[
            "supporting_evidence_layers"
        ]
    )

    print(
        "Tier:",
        row[
            "final_biological_tier"
        ]
    )

    print(
        "Interpretation:",
        row[
            "final_research_interpretation"
        ]
    )

    print(
        "Evidence convergence:",
        row[
            "evidence_convergence"
        ]
    )


# ============================================================
# 20. FINAL REPORT
# ============================================================

report_file = os.path.join(
    OUTPUT_DIR,
    "45_final_biological_interpretation_report.txt"
)


with open(
    report_file,
    "w"
) as f:

    f.write(
        "STEP 45 - FINAL BIOLOGICAL INTERPRETATION\n"
    )

    f.write(
        "=" * 70 + "\n\n"
    )

    f.write(
        "Purpose:\n"
    )

    f.write(
        "Interpret the final six candidate genes using "
        "integrated bioinformatics, ML, discovery expression, "
        "external validation, and tumor-association evidence.\n\n"
    )

    f.write(
        "FINAL CANDIDATE PANEL\n"
    )

    f.write(
        "-" * 70 + "\n"
    )


    for _, row in final_panel.iterrows():

        f.write(
            f"\nRank {int(row['final_research_rank'])}: "
            f"{row['gene']}\n"
        )

        f.write(
            f"Multi-layer score: "
            f"{row['final_multilayer_score']}/100\n"
        )

        f.write(
            f"Evidence layers: "
            f"{row['supporting_evidence_layers']}/5\n"
        )

        f.write(
            f"Final tier: "
            f"{row['final_biological_tier']}\n"
        )

        f.write(
            f"Biological interpretation: "
            f"{row['biological_interpretation']}\n"
        )

        f.write(
            f"Research interpretation: "
            f"{row['final_research_interpretation']}\n"
        )

        f.write(
            f"Evidence convergence: "
            f"{row['evidence_convergence']}\n"
        )


    f.write(
        "\n\nKEY BIOLOGICAL CONCLUSIONS\n"
    )

    f.write(
        "-" * 70 + "\n\n"
    )

    f.write(
        "PRSS23, FCGRT, ETFB and MRPS7 represent the "
        "strongest tumor-associated research candidates "
        "based on the evidence generated in this project.\n\n"
    )

    f.write(
        "EPN3 remains a promising candidate but requires "
        "additional validation because its consistency "
        "was weaker than the strongest candidates.\n\n"
    )

    f.write(
        "COX8A shows strong expression and computational "
        "support but relatively lower tumor specificity "
        "because it is also highly detected in diploid cells.\n\n"
    )

    f.write(
        "IMPORTANT LIMITATION:\n"
    )

    f.write(
        "These findings represent computational and "
        "single-cell research evidence. They do not "
        "establish clinical diagnostic, prognostic, "
        "therapeutic, or causal validity.\n"
    )

    f.write(
        "Experimental validation is required before "
        "clinical interpretation.\n"
    )


# ============================================================
# 21. FIGURE 1
# FINAL SCORES
# ============================================================

fig, ax = plt.subplots(
    figsize=(9, 6)
)


ax.bar(
    final_panel["gene"],
    final_panel[
        "final_multilayer_score"
    ]
)


ax.set_ylabel(
    "Exploratory evidence score / 100"
)

ax.set_xlabel(
    "Gene"
)

ax.set_title(
    "Final Candidate Evidence Scores"
)

ax.set_ylim(
    0,
    105
)


plt.tight_layout()


figure1 = os.path.join(
    FIGURE_DIR,
    "45_final_candidate_scores.png"
)


plt.savefig(
    figure1,
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# ============================================================
# 22. FIGURE 2
# MALIGNANT DETECTION DIFFERENCE
# ============================================================

if (
    "malignant_detection_difference_pct"
    in final_panel.columns
):

    fig, ax = plt.subplots(
        figsize=(9, 6)
    )


    ax.bar(
        final_panel["gene"],
        final_panel[
            "malignant_detection_difference_pct"
        ]
    )


    ax.axhline(
        0,
        linewidth=1
    )


    ax.set_ylabel(
        "Malignant detection difference (percentage points)"
    )

    ax.set_xlabel(
        "Gene"
    )

    ax.set_title(
        "Tumor-Associated Detection Enrichment"
    )


    plt.tight_layout()


    figure2 = os.path.join(
        FIGURE_DIR,
        "45_tumor_association_detection_difference.png"
    )


    plt.savefig(
        figure2,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()


# ============================================================
# 23. FIGURE 3
# EVIDENCE LAYER COUNT
# ============================================================

fig, ax = plt.subplots(
    figsize=(9, 6)
)


ax.bar(
    final_panel["gene"],
    final_panel[
        "supporting_evidence_layers"
    ]
)


ax.set_ylabel(
    "Supporting evidence layers"
)

ax.set_xlabel(
    "Gene"
)

ax.set_title(
    "Evidence Convergence Across Final Candidates"
)

ax.set_ylim(
    0,
    5.5
)


plt.tight_layout()


figure3 = os.path.join(
    FIGURE_DIR,
    "45_evidence_layer_convergence.png"
)


plt.savefig(
    figure3,
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# ============================================================
# 24. COMPLETION
# ============================================================

print(
    "\n" + "=" * 70
)

print(
    "STEP 45 COMPLETED"
)

print(
    "=" * 70
)

print(
    "\nFinal complete panel:"
)

print(
    complete_output
)

print(
    "\nFinal research candidate panel:"
)

print(
    panel_output
)

print(
    "\nFinal report:"
)

print(
    report_file
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