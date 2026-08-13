# ============================================================
# STEP 44
# FINAL MULTI-LAYER EVIDENCE INTEGRATION
# ============================================================
#
# Purpose:
#
# Integrate all major evidence generated during the project:
#
# 1. Original bioinformatics prioritization
# 2. ML recurrence / consensus
# 3. Discovery-dataset expression validation
# 4. GSE176078 external validation
# 5. Tumor-specificity validation
#
# This is an evidence-integration step.
# No new ML model is trained.
#
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# 1. INPUT FILES
# ============================================================

STEP39_FILE = (
    "results/ai_ml/39_consensus_visualization/"
    "39_final_candidate_panel.csv"
)

STEP41_FILE = (
    "results/ai_ml/41_bioinformatics_ml_cross_validation/"
    "41_final_six_evidence_matrix.csv"
)

STEP42_FILE = (
    "results/ai_ml/42_external_validation/"
    "42_six_gene_external_validation.csv"
)

STEP43_FILE = (
    "results/ai_ml/43_celltype_tumor_specificity/"
    "43_celltype_tumor_specificity.csv"
)


# ============================================================
# 2. OUTPUT DIRECTORY
# ============================================================

OUTPUT_DIR = (
    "results/ai_ml/44_final_multilayer_evidence"
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
# 3. FINAL SIX GENES
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
print("STEP 44 - FINAL MULTI-LAYER EVIDENCE INTEGRATION")
print("=" * 70)


# ============================================================
# 5. FILE CHECK
# ============================================================

print("\n===== INPUT FILE CHECK =====")

required_files = {
    "Step 39": STEP39_FILE,
    "Step 41": STEP41_FILE,
    "Step 42": STEP42_FILE,
    "Step 43": STEP43_FILE
}

for name, path in required_files.items():

    if os.path.exists(path):

        print(
            f"{name}: FOUND"
        )

    else:

        raise FileNotFoundError(
            f"\nMissing {name} file:\n{path}"
        )


# ============================================================
# 6. LOAD STEP 39
# ============================================================

print("\n===== LOADING STEP 39 =====")

step39 = pd.read_csv(
    STEP39_FILE
)

step39["gene"] = (
    step39["gene"]
    .astype(str)
    .str.strip()
)

step39 = step39[
    step39["gene"].isin(FINAL_SIX)
].copy()

print(
    "Step 39 genes:",
    len(step39)
)


# ============================================================
# 7. LOAD STEP 41
# ============================================================

print("\n===== LOADING STEP 41 =====")

step41 = pd.read_csv(
    STEP41_FILE
)

step41["gene"] = (
    step41["gene"]
    .astype(str)
    .str.strip()
)

step41 = step41[
    step41["gene"].isin(FINAL_SIX)
].copy()

print(
    "Step 41 genes:",
    len(step41)
)


# ============================================================
# 8. LOAD STEP 42
# ============================================================

print("\n===== LOADING STEP 42 =====")

step42 = pd.read_csv(
    STEP42_FILE
)

step42["gene"] = (
    step42["gene"]
    .astype(str)
    .str.strip()
)

step42 = step42[
    step42["gene"].isin(FINAL_SIX)
].copy()

print(
    "Step 42 genes:",
    len(step42)
)


# ============================================================
# 9. LOAD STEP 43
# ============================================================

print("\n===== LOADING STEP 43 =====")

step43 = pd.read_csv(
    STEP43_FILE
)

step43["gene"] = (
    step43["gene"]
    .astype(str)
    .str.strip()
)

step43 = step43[
    step43["gene"].isin(FINAL_SIX)
].copy()

print(
    "Step 43 genes:",
    len(step43)
)


# ============================================================
# 10. CREATE MASTER TABLE
# ============================================================

print(
    "\n===== BUILDING MASTER EVIDENCE TABLE ====="
)

master = pd.DataFrame({
    "gene": FINAL_SIX
})


# ============================================================
# 11. STEP 39 EVIDENCE
# ============================================================

step39_columns = [
    "final_evidence_rank",
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

available = [
    c
    for c in step39_columns
    if c in step39.columns
]

if available:

    temp = step39[
        ["gene"] + available
    ].copy()

    master = master.merge(
        temp,
        on="gene",
        how="left"
    )


# ============================================================
# 12. STEP 41 EVIDENCE
# ============================================================

step41_columns = [
    "bioinformatics_final_50",
    "priority_rank",
    "biomarker_score",
    "adjusted_biomarker_score",
    "biological_class",
    "final_biological_category",
    "ML_folds_selected",
    "ML_consensus_status",
    "expression_mean_log2FC",
    "expression_consistency",
    "expression_validation_score",
    "cross_validation_class"
]

available = [
    c
    for c in step41_columns
    if c in step41.columns
]

temp = step41[
    ["gene"] + available
].copy()

# Prevent duplicate ML columns from creating _x/_y
for column in temp.columns:

    if column != "gene" and column in master.columns:

        temp = temp.rename(
            columns={
                column:
                f"{column}_step41"
            }
        )

master = master.merge(
    temp,
    on="gene",
    how="left"
)


# ============================================================
# 13. STEP 42 EXTERNAL VALIDATION
# ============================================================

step42_columns = [
    "external_gene_present",
    "external_validated_gene",
    "external_gene_not_found",
    "external_support_category",
    "external_validation_score",
    "external_expression_prevalence",
    "external_validation_status",
    "evidence_layers",
    "final_external_validation_class"
]

available = [
    c
    for c in step42_columns
    if c in step42.columns
]

temp = step42[
    ["gene"] + available
].copy()

master = master.merge(
    temp,
    on="gene",
    how="left"
)


# ============================================================
# 14. STEP 43 TUMOR SPECIFICITY
# ============================================================

step43_columns = [
    "malignant_detection_pct",
    "diploid_detection_pct",
    "malignant_detection_difference_pct",
    "malignant_vs_diploid_log2FC",
    "epithelial_detection_pct",
    "non_epithelial_detection_pct",
    "epithelial_detection_difference_pct",
    "malignant_epithelial_detection_pct",
    "samples_tested",
    "samples_positive_log2FC",
    "sample_consistency",
    "mean_sample_log2FC",
    "tumor_specificity_class",
    "final_specificity_interpretation"
]

available = [
    c
    for c in step43_columns
    if c in step43.columns
]

temp = step43[
    ["gene"] + available
].copy()

master = master.merge(
    temp,
    on="gene",
    how="left"
)


# ============================================================
# 15. CLEAN DUPLICATE COLUMNS
# ============================================================

print(
    "\n===== CLEANING MASTER TABLE ====="
)

# Remove accidental pandas duplicate suffixes where possible
duplicate_suffixes = [
    "_x",
    "_y"
]

for base in [
    "ML_folds_selected",
    "expression_consistency",
    "samples_tested"
]:

    x = f"{base}_x"
    y = f"{base}_y"

    if x in master.columns and y in master.columns:

        master[base] = master[x].combine_first(
            master[y]
        )

        master = master.drop(
            columns=[
                x,
                y
            ]
        )

    elif x in master.columns:

        master = master.rename(
            columns={
                x: base
            }
        )

    elif y in master.columns:

        master = master.rename(
            columns={
                y: base
            }
        )


# ============================================================
# 16. DEFINE BOOLEAN EVIDENCE FLAGS
# ============================================================

print(
    "\n===== DEFINING EVIDENCE FLAGS ====="
)


master[
    "bioinformatics_support"
] = (
    master[
        "bioinformatics_final_50"
    ]
    .fillna(False)
    .astype(bool)
)


master[
    "ML_3_of_3"
] = (
    pd.to_numeric(
        master[
            "ML_folds_selected"
        ],
        errors="coerce"
    )
    == 3
)


master[
    "discovery_expression_support"
] = (
    master[
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


master[
    "external_support_flag"
] = (
    master[
        "external_validation_score"
    ]
    .fillna(0)
    >= 2
)


master[
    "tumor_association_flag"
] = (
    master[
        "tumor_specificity_class"
    ]
    .isin(
        [
            "Tumor_associated",
            "Strong_tumor_epithelial_specificity"
        ]
    )
)


# ============================================================
# 17. EVIDENCE LAYER COUNT
# ============================================================

master[
    "supporting_evidence_layers"
] = (

    master[
        "bioinformatics_support"
    ].astype(int)

    +

    master[
        "ML_3_of_3"
    ].astype(int)

    +

    master[
        "discovery_expression_support"
    ].astype(int)

    +

    master[
        "external_support_flag"
    ].astype(int)

    +

    master[
        "tumor_association_flag"
    ].astype(int)

)


# ============================================================
# 18. CREATE FINAL EVIDENCE SCORE
# ============================================================
#
# Maximum = 100
#
# Components:
#
# Bioinformatics support       20
# ML recurrence                20
# Discovery expression        20
# External validation         20
# Tumor association           20
#
# This is an exploratory project score,
# NOT a clinical biomarker score.
#
# ============================================================

master[
    "bioinformatics_score_component"
] = np.where(
    master[
        "bioinformatics_support"
    ],
    20,
    0
)


master[
    "ML_score_component"
] = np.where(
    master[
        "ML_3_of_3"
    ],
    20,
    0
)


master[
    "expression_score_component"
] = np.where(
    master[
        "discovery_expression_support"
    ],
    20,
    0
)


master[
    "external_score_component"
] = np.where(
    master[
        "external_support_flag"
    ],
    20,
    0
)


master[
    "tumor_specificity_score_component"
] = np.where(
    master[
        "tumor_association_flag"
    ],
    20,
    0
)


master[
    "final_multilayer_score"
] = (

    master[
        "bioinformatics_score_component"
    ]

    +

    master[
        "ML_score_component"
    ]

    +

    master[
        "expression_score_component"
    ]

    +

    master[
        "external_score_component"
    ]

    +

    master[
        "tumor_specificity_score_component"
    ]

)


# ============================================================
# 19. FINAL CANDIDATE CLASS
# ============================================================

def classify_candidate(row):

    score = row[
        "final_multilayer_score"
    ]

    tumor_class = row[
        "tumor_specificity_class"
    ]

    external = row[
        "external_support_flag"
    ]

    expression = row[
        "discovery_expression_support"
    ]

    ml = row[
        "ML_3_of_3"
    ]


    # Highest-confidence computational candidate
    if (
        score >= 80
        and
        ml
        and
        expression
        and
        external
    ):

        return (
            "High-priority integrated candidate"
        )


    # Strong but missing one evidence layer
    if score >= 60:

        return (
            "Strong integrated candidate"
        )


    # Intermediate
    if score >= 40:

        return (
            "Moderate integrated candidate"
        )


    return (
        "Candidate requiring further validation"
    )


master[
    "final_candidate_class"
] = master.apply(
    classify_candidate,
    axis=1
)


# ============================================================
# 20. EPN3 SPECIAL FLAG
# ============================================================

master[
    "additional_validation_flag"
] = ""

master.loc[
    master["gene"] == "EPN3",
    "additional_validation_flag"
] = (
    "Weaker external/sample consistency; "
    "requires further validation"
)


# ============================================================
# 21. COX8A SPECIAL INTERPRETATION
# ============================================================

master.loc[
    master["gene"] == "COX8A",
    "additional_validation_flag"
] = (
    "Highly expressed but relatively broad expression; "
    "lower tumor specificity"
)


# ============================================================
# 22. SORT FINAL TABLE
# ============================================================

master = master.sort_values(
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


master[
    "final_multilayer_rank"
] = (
    master.index + 1
)


# ============================================================
# 23. SAVE COMPLETE TABLE
# ============================================================

complete_file = os.path.join(
    OUTPUT_DIR,
    "44_complete_multilayer_evidence_table.csv"
)

master.to_csv(
    complete_file,
    index=False
)


# ============================================================
# 24. CREATE FINAL COMPACT TABLE
# ============================================================

compact_columns = [
    "final_multilayer_rank",
    "gene",

    "bioinformatics_support",
    "ML_3_of_3",
    "discovery_expression_support",
    "external_support_flag",
    "tumor_association_flag",

    "supporting_evidence_layers",
    "final_multilayer_score",

    "mean_log2FC",
    "malignant_detection_pct",
    "diploid_detection_pct",
    "malignant_detection_difference_pct",

    "external_support_category",
    "tumor_specificity_class",

    "final_candidate_class",
    "additional_validation_flag"
]

compact_columns = [
    c
    for c in compact_columns
    if c in master.columns
]


compact = master[
    compact_columns
].copy()


compact_file = os.path.join(
    OUTPUT_DIR,
    "44_final_candidate_evidence_matrix.csv"
)

compact.to_csv(
    compact_file,
    index=False
)


# ============================================================
# 25. PRINT FINAL TABLE
# ============================================================

print(
    "\n" + "=" * 70
)

print(
    "FINAL MULTI-LAYER EVIDENCE MATRIX"
)

print(
    "=" * 70
)

print(
    compact.to_string(
        index=False
    )
)


# ============================================================
# 26. SUMMARY
# ============================================================

print(
    "\n===== FINAL PROJECT SUMMARY ====="
)


print(
    "Genes:",
    len(master)
)


print(
    "\nEvidence layers:"
)

print(
    "Bioinformatics supported:",
    int(
        master[
            "bioinformatics_support"
        ].sum()
    )
)

print(
    "ML 3/3:",
    int(
        master[
            "ML_3_of_3"
        ].sum()
    )
)

print(
    "Discovery expression supported:",
    int(
        master[
            "discovery_expression_support"
        ].sum()
    )
)

print(
    "External validation supported:",
    int(
        master[
            "external_support_flag"
        ].sum()
    )
)

print(
    "Tumor-associated:",
    int(
        master[
            "tumor_association_flag"
        ].sum()
    )
)


# ============================================================
# 27. FINAL CLASS DISTRIBUTION
# ============================================================

print(
    "\n===== FINAL CANDIDATE CLASSIFICATION ====="
)

print(
    master[
        "final_candidate_class"
    ]
    .value_counts()
    .to_string()
)


# ============================================================
# 28. SAVE SUMMARY REPORT
# ============================================================

summary_file = os.path.join(
    OUTPUT_DIR,
    "44_final_multilayer_evidence_report.txt"
)


with open(
    summary_file,
    "w"
) as f:

    f.write(
        "STEP 44 - FINAL MULTI-LAYER EVIDENCE INTEGRATION\n"
    )

    f.write(
        "=" * 70 + "\n\n"
    )

    f.write(
        "Purpose:\n"
    )

    f.write(
        "Integrate bioinformatics, machine-learning, "
        "discovery expression, external validation, "
        "and tumor-specificity evidence for the final "
        "six candidate genes.\n\n"
    )

    f.write(
        "IMPORTANT:\n"
    )

    f.write(
        "The final score is an exploratory computational "
        "evidence score and is NOT a clinical biomarker "
        "score or probability of clinical validity.\n\n"
    )

    f.write(
        "FINAL CANDIDATES\n"
    )

    f.write(
        "-" * 70 + "\n"
    )


    for _, row in master.iterrows():

        f.write(
            f"\n{row['final_multilayer_rank']}. "
            f"{row['gene']}\n"
        )

        f.write(
            f"Final multi-layer score: "
            f"{row['final_multilayer_score']}/100\n"
        )

        f.write(
            f"Evidence layers: "
            f"{row['supporting_evidence_layers']}/5\n"
        )

        f.write(
            f"Bioinformatics: "
            f"{row['bioinformatics_support']}\n"
        )

        f.write(
            f"ML 3/3: "
            f"{row['ML_3_of_3']}\n"
        )

        f.write(
            f"Discovery expression: "
            f"{row['discovery_expression_support']}\n"
        )

        f.write(
            f"External validation: "
            f"{row['external_support_flag']}\n"
        )

        f.write(
            f"Tumor association: "
            f"{row['tumor_association_flag']}\n"
        )

        f.write(
            f"Final classification: "
            f"{row['final_candidate_class']}\n"
        )

        if row[
            "additional_validation_flag"
        ]:

            f.write(
                f"Additional note: "
                f"{row['additional_validation_flag']}\n"
            )


    f.write(
        "\n\nLIMITATIONS\n"
    )

    f.write(
        "-" * 70 + "\n"
    )

    f.write(
        "These results represent computational and "
        "single-cell evidence. They do not establish "
        "clinical diagnostic, prognostic, therapeutic, "
        "or causal validity.\n"
    )

    f.write(
        "Experimental validation would be required "
        "before clinical interpretation.\n"
    )


# ============================================================
# 29. FIGURE 1
# FINAL SCORE
# ============================================================

print(
    "\n===== GENERATING FIGURES ====="
)


fig, ax = plt.subplots(
    figsize=(9, 6)
)


ax.bar(
    master[
        "gene"
    ],
    master[
        "final_multilayer_score"
    ]
)


ax.set_ylabel(
    "Exploratory multi-layer evidence score / 100"
)

ax.set_xlabel(
    "Gene"
)

ax.set_title(
    "Final Multi-Layer Evidence Score"
)

ax.set_ylim(
    0,
    105
)

plt.tight_layout()


figure1 = os.path.join(
    FIGURE_DIR,
    "44_final_multilayer_scores.png"
)


plt.savefig(
    figure1,
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# ============================================================
# 30. FIGURE 2
# EVIDENCE LAYER HEATMAP
# ============================================================

heatmap_columns = [
    "bioinformatics_support",
    "ML_3_of_3",
    "discovery_expression_support",
    "external_support_flag",
    "tumor_association_flag"
]


heatmap = (
    master[
        [
            "gene"
        ] + heatmap_columns
    ]
    .set_index(
        "gene"
    )
    .astype(int)
)


fig, ax = plt.subplots(
    figsize=(10, 5)
)


image = ax.imshow(
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
        "ML 3/3",
        "Discovery expression",
        "External validation",
        "Tumor association"
    ],
    rotation=35,
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
    "Evidence Layer Support Across Final Candidates"
)


plt.tight_layout()


figure2 = os.path.join(
    FIGURE_DIR,
    "44_evidence_layer_heatmap.png"
)


plt.savefig(
    figure2,
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# ============================================================
# 31. FIGURE 3
# MALIGNANT DETECTION
# ============================================================

if (
    "malignant_detection_pct"
    in master.columns
):

    fig, ax = plt.subplots(
        figsize=(9, 6)
    )


    ax.bar(
        master[
            "gene"
        ],
        master[
            "malignant_detection_pct"
        ]
    )


    ax.set_ylabel(
        "Malignant-cell detection (%)"
    )

    ax.set_xlabel(
        "Gene"
    )

    ax.set_title(
        "Final Candidates: Malignant-Cell Detection"
    )


    ax.set_ylim(
        0,
        105
    )


    plt.tight_layout()


    figure3 = os.path.join(
        FIGURE_DIR,
        "44_malignant_detection.png"
    )


    plt.savefig(
        figure3,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

else:

    figure3 = None


# ============================================================
# 32. COMPLETION
# ============================================================

print(
    "\n" + "=" * 70
)

print(
    "STEP 44 COMPLETED"
)

print(
    "=" * 70
)

print(
    "\nMain complete table:"
)

print(
    complete_file
)

print(
    "\nFinal compact evidence matrix:"
)

print(
    compact_file
)

print(
    "\nReport:"
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

if figure3:

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