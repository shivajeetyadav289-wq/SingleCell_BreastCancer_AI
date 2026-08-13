# ============================================================
# STEP 40
# BIOLOGICAL INTERPRETATION & MULTI-LAYER EVIDENCE INTEGRATION
# ============================================================
#
# Purpose:
#
# Integrate the original bioinformatics biological annotation
# with:
#
#   Step 36 - Bioinformatics + ML consensus
#   Step 38 - Expression validation
#   Step 39 - Final integrated evidence
#
# The goal is to produce a biologically interpretable final
# six-gene candidate panel.
#
# IMPORTANT:
#
# These genes are candidate biomarkers.
# They are NOT clinically validated biomarkers.
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

STEP36_FILE = (
    "results/ai_ml/36_bioinformatics_ml_consensus/"
    "36_consensus_candidate_summary.csv"
)

STEP37_FILE = (
    "results/ai_ml/37_consensus_characterization/"
    "37_consensus_candidate_characterization.csv"
)

STEP38_FILE = (
    "results/ai_ml/38_consensus_expression_validation/"
    "38_consensus_expression_summary.csv"
)

STEP39_FILE = (
    "results/ai_ml/39_consensus_visualization/"
    "39_final_integrated_evidence_table.csv"
)


OUTPUT_DIR = (
    "results/ai_ml/40_biological_interpretation"
)

FIGURE_DIR = os.path.join(
    OUTPUT_DIR,
    "figures"
)


# ============================================================
# 2. FINAL SIX GENES
# ============================================================

FINAL_GENES = [

    "COX8A",
    "MRPS7",
    "PRSS23",
    "ETFB",
    "FCGRT",
    "EPN3"

]


# ============================================================
# 3. CREATE OUTPUT DIRECTORIES
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
# 4. START
# ============================================================

print("=" * 70)

print(
    "STEP 40 - BIOLOGICAL INTERPRETATION"
)

print("=" * 70)


# ============================================================
# 5. INPUT CHECK
# ============================================================

print(
    "\n===== INPUT FILE CHECK ====="
)


required_files = {

    "Bioinformatics annotation":
        BIOINFO_FILE,

    "Step 38 expression validation":
        STEP38_FILE,

    "Step 39 integrated evidence":
        STEP39_FILE

}


optional_files = {

    "Step 36 consensus":
        STEP36_FILE,

    "Step 37 characterization":
        STEP37_FILE

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


for name, path in optional_files.items():

    if os.path.exists(path):

        print(
            f"{name}: FOUND"
        )

    else:

        print(
            f"{name}: NOT FOUND"
        )


# ============================================================
# 6. LOAD BIOINFORMATICS ANNOTATION
# ============================================================

print(
    "\n===== LOADING BIOINFORMATICS ANNOTATION ====="
)


bioinfo = pd.read_csv(
    BIOINFO_FILE
)


print(
    "Rows:",
    len(bioinfo)
)

print(
    "Columns:",
    len(bioinfo.columns)
)


if "gene" not in bioinfo.columns:

    raise ValueError(
        "Bioinformatics file does not contain 'gene'."
    )


bioinfo["gene"] = (
    bioinfo["gene"]
    .astype(str)
    .str.strip()
)


# ============================================================
# 7. CHECK FINAL GENES
# ============================================================

bioinfo_genes = set(
    bioinfo["gene"]
)


missing_bioinfo = [

    gene
    for gene in FINAL_GENES
    if gene not in bioinfo_genes

]


print(
    "\n===== FINAL GENE CHECK ====="
)

print(
    "Final genes:",
    len(FINAL_GENES)
)

print(
    "Present in bioinformatics annotation:",
    len(FINAL_GENES) - len(missing_bioinfo)
)


if missing_bioinfo:

    print(
        "\nMissing:"
    )

    for gene in missing_bioinfo:

        print(
            " -",
            gene
        )


# ============================================================
# 8. SELECT BIOINFORMATICS COLUMNS
# ============================================================

bioinfo_columns = [

    "gene",

    "p_val",

    "avg_log2FC",

    "pct.1",

    "pct.2",

    "p_val_adj",

    "prevalence_difference",

    "absolute_log2FC",

    "direction",

    "samples_tested",

    "samples_prevalence_25pct",

    "samples_prevalence_50pct",

    "minimum_sample_prevalence",

    "mean_sample_prevalence",

    "maximum_sample_prevalence",

    "consistency_category",

    "biological_review_flag",

    "fc_score",

    "prevalence_score",

    "specificity_score",

    "significance_score",

    "consistency_score",

    "biomarker_score",

    "adjusted_biomarker_score",

    "priority_rank",

    "biological_class",

    "final_biological_category"

]


bioinfo_columns = [

    column
    for column in bioinfo_columns
    if column in bioinfo.columns

]


bioinfo_final = (
    bioinfo[
        bioinfo["gene"].isin(
            FINAL_GENES
        )
    ][
        bioinfo_columns
    ]
    .drop_duplicates(
        subset=["gene"]
    )
    .copy()
)


# ============================================================
# 9. LOAD STEP 38
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
            FINAL_GENES
        )
    ]
    .copy()
)


print(
    "Step 38 genes:",
    len(
        step38_final
    )
)


# ============================================================
# 10. LOAD STEP 39
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
            FINAL_GENES
        )
    ]
    .copy()
)


print(
    "Step 39 genes:",
    len(
        step39_final
    )
)


# ============================================================
# 11. LOAD STEP 37 IF AVAILABLE
# ============================================================

step37_final = None


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


    step37_final = (
        step37[
            step37["gene"].isin(
                FINAL_GENES
            )
        ]
        .drop_duplicates(
            subset=["gene"]
        )
        .copy()
    )


    print(
        "Step 37 genes:",
        len(
            step37_final
        )
    )


# ============================================================
# 12. MERGE BIOINFORMATICS + STEP 38 + STEP 39
# ============================================================

print(
    "\n===== INTEGRATING EVIDENCE ====="
)


integrated = bioinfo_final.copy()


# ------------------------------------------------------------
# Add Step 38
# ------------------------------------------------------------

step38_columns = [

    "gene",

    "mean_malignant_expression",

    "mean_diploid_expression",

    "mean_log2FC",

    "mean_malignant_detection_pct",

    "mean_diploid_detection_pct",

    "mean_detection_difference_pct",

    "samples_with_positive_log2FC",

    "samples_with_log2FC_ge_1",

    "samples_with_positive_detection_difference",

    "sample_consistency_score",

    "expression_consistency",

    "expression_validation_score"

]


step38_columns = [

    column
    for column in step38_columns
    if column in step38_final.columns

]


integrated = integrated.merge(

    step38_final[
        step38_columns
    ],

    on="gene",

    how="left",

    suffixes=(
        "",
        "_step38"
    )

)


# ------------------------------------------------------------
# Add Step 39
# ------------------------------------------------------------

step39_columns = [

    "gene",

    "ML_folds_selected",

    "final_evidence_category",

    "overall_exploratory_evidence_score",

    "external_support"

]


step39_columns = [

    column
    for column in step39_columns
    if column in step39_final.columns

]


integrated = integrated.merge(

    step39_final[
        step39_columns
    ],

    on="gene",

    how="left",

    suffixes=(
        "",
        "_step39"
    )

)


# ============================================================
# 13. ADD STEP 37 BIOLOGICAL INFORMATION
# ============================================================

if step37_final is not None:

    step37_columns = [

        "gene",

        "biological_category",

        "ML_folds_selected",

        "ML_recurrence",

        "external_expression_prevalence",

        "external_mean_expression",

        "external_support"

    ]


    step37_columns = [

        column
        for column in step37_columns
        if column in step37_final.columns

    ]


    if len(
        step37_columns
    ) > 1:

        integrated = integrated.merge(

            step37_final[
                step37_columns
            ],

            on="gene",

            how="left",

            suffixes=(
                "",
                "_step37"
            )

        )


# ============================================================
# 14. REMOVE DUPLICATE COLUMN CONFLICTS
# ============================================================

# Some previous steps may contain the same concept under
# slightly different names. Keep the Step 40 representation
# simple and explicit.


# ============================================================
# 15. CREATE BIOLOGICAL INTERPRETATION
# ============================================================

def biological_interpretation(row):

    category = str(
        row.get(
            "final_biological_category",
            ""
        )
    )


    biological_class = str(
        row.get(
            "biological_class",
            ""
        )
    )


    gene = str(
        row.get(
            "gene",
            ""
        )
    )


    # --------------------------------------------------------
    # Use categories generated by the original pipeline.
    # Do NOT invent categories.
    # --------------------------------------------------------

    if (
        "OXPHOS"
        in category
        or
        "Mitochondrial"
        in category
        or
        "OXPHOS"
        in biological_class
        or
        "Mitochondrial"
        in biological_class
    ):

        return (
            "Mitochondrial/OXPHOS-associated "
            "candidate with multi-layer support."
        )


    if (
        "Epithelial"
        in category
        or
        "Tumor"
        in category
        or
        "Epithelial"
        in biological_class
    ):

        return (
            "Epithelial/tumor-associated "
            "candidate with multi-layer support."
        )


    if (
        "Proliferation"
        in category
        or
        "Proliferation"
        in biological_class
    ):

        return (
            "Proliferation-associated "
            "candidate with multi-layer support."
        )


    if (
        "Stress"
        in category
        or
        "Stress"
        in biological_class
    ):

        return (
            "Stress-response-associated "
            "candidate with multi-layer support."
        )


    return (
        "Bioinformatics-prioritized candidate "
        "with ML and expression evidence."
    )


integrated[
    "biological_interpretation"
] = integrated.apply(
    biological_interpretation,
    axis=1
)


# ============================================================
# 16. CREATE EVIDENCE FLAGS
# ============================================================

integrated[
    "bioinformatics_supported"
] = True


integrated[
    "ML_3_of_3"
] = (

    pd.to_numeric(
        integrated.get(
            "ML_folds_selected",
            np.nan
        ),
        errors="coerce"
    )
    >= 3

)


integrated[
    "expression_3_of_3"
] = (

    pd.to_numeric(
        integrated[
            "samples_with_positive_log2FC"
        ],
        errors="coerce"
    )

    >=

    pd.to_numeric(
        integrated[
            "samples_tested"
        ],
        errors="coerce"
    )

)


integrated[
    "strong_malignant_expression"
] = (

    pd.to_numeric(
        integrated[
            "mean_log2FC"
        ],
        errors="coerce"
    )
    >= 1

)


# ============================================================
# 17. EVIDENCE LAYER COUNT
# ============================================================

integrated[
    "evidence_layer_count"
] = (

    integrated[
        "bioinformatics_supported"
    ].astype(int)

    +

    integrated[
        "ML_3_of_3"
    ].astype(int)

    +

    integrated[
        "expression_3_of_3"
    ].astype(int)

    +

    integrated[
        "strong_malignant_expression"
    ].astype(int)

)


# ============================================================
# 18. FINAL BIOLOGICAL PRIORITY
# ============================================================

def final_priority(row):

    layers = row[
        "evidence_layer_count"
    ]


    gene = row[
        "gene"
    ]


    if layers >= 4:

        return (
            "Strong multi-layer candidate"
        )


    if layers == 3:

        return (
            "Supported multi-layer candidate"
        )


    if layers == 2:

        return (
            "Moderate candidate"
        )


    return (
        "Candidate requiring validation"
    )


integrated[
    "final_biological_priority"
] = integrated.apply(
    final_priority,
    axis=1
)


# ============================================================
# 19. SORT FINAL RESULTS
# ============================================================

integrated = (
    integrated
    .sort_values(
        [
            "evidence_layer_count",
            "adjusted_biomarker_score",
            "mean_log2FC"
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


integrated[
    "biological_interpretation_rank"
] = (
    integrated.index + 1
)


# ============================================================
# 20. SAVE FINAL TABLE
# ============================================================

output_file = os.path.join(
    OUTPUT_DIR,
    "40_final_gene_biological_interpretation.csv"
)


integrated.to_csv(
    output_file,
    index=False
)


print(
    "\nFinal interpretation table saved:"
)

print(
    output_file
)


# ============================================================
# 21. DISPLAY CORE RESULT
# ============================================================

display_columns = [

    "biological_interpretation_rank",

    "gene",

    "biological_class",

    "final_biological_category",

    "adjusted_biomarker_score",

    "priority_rank",

    "ML_folds_selected",

    "mean_log2FC",

    "mean_malignant_detection_pct",

    "mean_diploid_detection_pct",

    "samples_with_positive_log2FC",

    "samples_tested",

    "expression_consistency",

    "evidence_layer_count",

    "final_biological_priority",

    "biological_interpretation"

]


display_columns = [

    column
    for column in display_columns
    if column in integrated.columns

]


print(
    "\n===== FINAL BIOLOGICAL INTERPRETATION ====="
)


print(
    integrated[
        display_columns
    ].to_string(
        index=False
    )
)


# ============================================================
# 22. BIOLOGICAL CATEGORY SUMMARY
# ============================================================

print(
    "\n===== BIOLOGICAL CATEGORY SUMMARY ====="
)


category_column = None


if (
    "final_biological_category"
    in integrated.columns
):

    category_column = (
        "final_biological_category"
    )

elif (
    "biological_class"
    in integrated.columns
):

    category_column = (
        "biological_class"
    )


if category_column is not None:

    category_summary = (

        integrated[
            category_column
        ]
        .fillna(
            "Unknown"
        )
        .value_counts()
        .reset_index()

    )


    category_summary.columns = [

        "biological_category",

        "gene_count"

    ]


    print(
        category_summary.to_string(
            index=False
        )
    )


    category_file = os.path.join(
        OUTPUT_DIR,
        "40_pathway_category_summary.csv"
    )


    category_summary.to_csv(
        category_file,
        index=False
    )


else:

    category_summary = pd.DataFrame()


# ============================================================
# 23. FIGURE 1
# EVIDENCE LAYER HEATMAP
# ============================================================

print(
    "\n===== GENERATING EVIDENCE HEATMAP ====="
)


heatmap_columns = [

    "bioinformatics_supported",

    "ML_3_of_3",

    "expression_3_of_3",

    "strong_malignant_expression"

]


heatmap_data = (
    integrated[
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
    heatmap_data.values,
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
        "Expression 3/3",
        "Strong expression"
    ],
    rotation=45,
    ha="right"
)


ax.set_yticks(
    range(
        len(
            heatmap_data.index
        )
    )
)


ax.set_yticklabels(
    heatmap_data.index
)


for i in range(
    heatmap_data.shape[0]
):

    for j in range(
        heatmap_data.shape[1]
    ):

        ax.text(
            j,
            i,
            str(
                heatmap_data.iloc[
                    i,
                    j
                ]
            ),
            ha="center",
            va="center"
        )


ax.set_title(
    "Evidence Layers Supporting Final Consensus Candidates"
)


plt.tight_layout()


heatmap_file = os.path.join(
    FIGURE_DIR,
    "40_candidate_evidence_heatmap.png"
)


plt.savefig(
    heatmap_file,
    dpi=300,
    bbox_inches="tight"
)


plt.close()


print(
    heatmap_file
)


# ============================================================
# 24. FIGURE 2
# BIOLOGICAL CATEGORY DISTRIBUTION
# ============================================================

if not category_summary.empty:

    print(
        "\n===== GENERATING CATEGORY PLOT ====="
    )


    fig, ax = plt.subplots(
        figsize=(9, 5)
    )


    ax.bar(
        category_summary[
            "biological_category"
        ],
        category_summary[
            "gene_count"
        ]
    )


    ax.set_ylabel(
        "Number of final candidates"
    )


    ax.set_xlabel(
        "Biological category"
    )


    ax.set_title(
        "Biological Categories of Final Consensus Candidates"
    )


    plt.xticks(
        rotation=45,
        ha="right"
    )


    plt.tight_layout()


    category_plot = os.path.join(
        FIGURE_DIR,
        "40_biological_category_distribution.png"
    )


    plt.savefig(
        category_plot,
        dpi=300,
        bbox_inches="tight"
    )


    plt.close()


    print(
        category_plot
    )


# ============================================================
# 25. CREATE BIOLOGICAL REPORT
# ============================================================

report_file = os.path.join(
    OUTPUT_DIR,
    "40_final_candidate_biological_report.txt"
)


with open(
    report_file,
    "w"
) as report:

    report.write(
        "STEP 40 - FINAL BIOLOGICAL INTERPRETATION\n"
    )

    report.write(
        "=" * 70
        + "\n\n"
    )


    report.write(
        "FINAL CONSENSUS GENES\n"
    )

    report.write(
        "-" * 70
        + "\n"
    )


    for _, row in integrated.iterrows():

        report.write(
            f"\nRank {int(row['biological_interpretation_rank'])}: "
            f"{row['gene']}\n"
        )


        if (
            "biological_class"
            in row.index
        ):

            report.write(
                f"Biological class: "
                f"{row['biological_class']}\n"
            )


        if (
            "final_biological_category"
            in row.index
        ):

            report.write(
                f"Final biological category: "
                f"{row['final_biological_category']}\n"
            )


        if (
            "adjusted_biomarker_score"
            in row.index
        ):

            report.write(
                f"Bioinformatics adjusted biomarker score: "
                f"{row['adjusted_biomarker_score']:.4f}\n"
            )


        if (
            "ML_folds_selected"
            in row.index
        ):

            report.write(
                f"ML recurrence: "
                f"{row['ML_folds_selected']}/3\n"
            )


        if (
            "mean_log2FC"
            in row.index
        ):

            report.write(
                f"Mean malignant-vs-diploid log2FC: "
                f"{row['mean_log2FC']:.4f}\n"
            )


        if (
            "mean_malignant_detection_pct"
            in row.index
        ):

            report.write(
                f"Malignant detection: "
                f"{row['mean_malignant_detection_pct']:.2f}%\n"
            )


        if (
            "mean_diploid_detection_pct"
            in row.index
        ):

            report.write(
                f"Diploid detection: "
                f"{row['mean_diploid_detection_pct']:.2f}%\n"
            )


        if (
            "samples_with_positive_log2FC"
            in row.index
        ):

            report.write(
                f"Positive expression samples: "
                f"{int(row['samples_with_positive_log2FC'])}/"
                f"{int(row['samples_tested'])}\n"
            )


        report.write(
            f"Evidence layers: "
            f"{int(row['evidence_layer_count'])}\n"
        )


        report.write(
            f"Final biological priority: "
            f"{row['final_biological_priority']}\n"
        )


        report.write(
            f"Interpretation: "
            f"{row['biological_interpretation']}\n"
        )


    report.write(
        "\n\nPROJECT-LEVEL INTERPRETATION\n"
    )

    report.write(
        "-" * 70
        + "\n"
    )


    report.write(
        "The final six genes represent a refined consensus "
        "candidate panel derived from the original "
        "bioinformatics candidate space and subsequently "
        "supported by machine-learning recurrence and "
        "single-cell expression validation.\n\n"
    )


    report.write(
        "The candidates should be considered research "
        "biomarker candidates rather than clinically validated "
        "biomarkers.\n\n"
    )


    report.write(
        "The strongest evidence is the convergence of "
        "independent computational layers rather than any "
        "single score.\n\n"
    )


    report.write(
        "LIMITATIONS\n"
    )

    report.write(
        "-" * 70
        + "\n"
    )


    report.write(
        "1. Only three informative patient samples were "
        "available for malignant-vs-diploid validation.\n"
    )

    report.write(
        "2. LOPO machine-learning performance showed limited "
        "generalization.\n"
    )

    report.write(
        "3. The integrated evidence score is exploratory.\n"
    )

    report.write(
        "4. Independent external datasets and experimental "
        "validation are required.\n"
    )


# ============================================================
# 26. FINAL OUTPUT
# ============================================================

print(
    "\n" + "=" * 70
)

print(
    "STEP 40 COMPLETED"
)

print(
    "=" * 70
)


print(
    "\nFinal interpretation:"
)

print(
    output_file
)


print(
    "\nCategory summary:"
)

if category_column is not None:

    print(
        category_file
    )


print(
    "\nBiological report:"
)

print(
    report_file
)


print(
    "\nFigures:"
)

print(
    heatmap_file
)

if not category_summary.empty:

    print(
        category_plot
    )


print(
    "\n" + "=" * 70
)