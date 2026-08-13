# ============================================================
# STEP 42
# ML CONSENSUS EXTERNAL VALIDATION
# ============================================================
#
# Purpose:
# Integrate the six ML consensus genes with the already
# completed GSE176078 external validation results.
#
# IMPORTANT:
# - GSE176078 is NOT downloaded again.
# - No new external expression analysis is performed.
# - Existing Step 26/27 results are reused.
#
# Final ML consensus genes:
# COX8A
# MRPS7
# PRSS23
# ETFB
# FCGRT
# EPN3
#
# ============================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# 1. PATHS
# ============================================================

STEP39_FILE = (
    "results/ai_ml/39_consensus_visualization/"
    "39_final_candidate_panel.csv"
)

STEP41_FILE = (
    "results/ai_ml/41_bioinformatics_ml_cross_validation/"
    "41_final_six_evidence_matrix.csv"
)

EXTERNAL_EXPRESSION_FILE = (
    "results/malignant/external_validation/"
    "27_external_expression_validation.csv"
)

EXTERNAL_SAMPLE_FILE = (
    "results/malignant/external_validation/"
    "27_external_sample_consistency.csv"
)

EXTERNAL_CANDIDATE_FILE = (
    "results/malignant/external_validation/"
    "26_external_validation_candidates.csv"
)

EXTERNAL_VALIDATED_LIST = (
    "results/malignant/external_validation/"
    "27_external_validated_gene_list.txt"
)

EXTERNAL_NOT_FOUND = (
    "results/malignant/external_validation/"
    "27_candidates_not_found_externally.txt"
)


# ============================================================
# 2. OUTPUT DIRECTORY
# ============================================================

OUTPUT_DIR = (
    "results/ai_ml/42_external_validation"
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
# 3. SIX ML CONSENSUS GENES
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
print("STEP 42 - ML CONSENSUS EXTERNAL VALIDATION")
print("=" * 70)


# ============================================================
# 5. FILE CHECK
# ============================================================

print("\n===== INPUT FILE CHECK =====")

required_files = {
    "Step 39 final candidate panel": STEP39_FILE,
    "Step 41 final six evidence matrix": STEP41_FILE,
    "External expression validation": EXTERNAL_EXPRESSION_FILE,
    "External sample consistency": EXTERNAL_SAMPLE_FILE,
    "External candidate validation": EXTERNAL_CANDIDATE_FILE,
    "External validated gene list": EXTERNAL_VALIDATED_LIST,
    "External not-found list": EXTERNAL_NOT_FOUND,
}

for name, path in required_files.items():

    if os.path.exists(path):
        print(f"{name}: FOUND")
    else:
        raise FileNotFoundError(
            f"\nMissing input file:\n{name}\n{path}"
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

step39_final = step39[
    step39["gene"].isin(FINAL_SIX)
].copy()

print(
    "Step 39 genes found:",
    len(step39_final)
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

step41_final = step41[
    step41["gene"].isin(FINAL_SIX)
].copy()

print(
    "Step 41 genes found:",
    len(step41_final)
)


# ============================================================
# 8. LOAD EXTERNAL EXPRESSION RESULTS
# ============================================================

print("\n===== LOADING GSE176078 EXPRESSION VALIDATION =====")

external_expression = pd.read_csv(
    EXTERNAL_EXPRESSION_FILE
)

external_expression["gene"] = (
    external_expression["gene"]
    .astype(str)
    .str.strip()
)

external_expression_final = external_expression[
    external_expression["gene"].isin(FINAL_SIX)
].copy()

print(
    "Six genes present in external expression table:",
    len(external_expression_final)
)


# ============================================================
# 9. LOAD EXTERNAL SAMPLE CONSISTENCY
# ============================================================

print("\n===== LOADING EXTERNAL SAMPLE CONSISTENCY =====")

external_sample = pd.read_csv(
    EXTERNAL_SAMPLE_FILE
)

external_sample["gene"] = (
    external_sample["gene"]
    .astype(str)
    .str.strip()
)

external_sample_final = external_sample[
    external_sample["gene"].isin(FINAL_SIX)
].copy()


print(
    "Genes with sample-level external data:",
    len(external_sample_final)
)


# ============================================================
# 10. LOAD EXTERNAL CANDIDATE VALIDATION
# ============================================================

print("\n===== LOADING EXTERNAL CANDIDATE VALIDATION =====")

external_candidates = pd.read_csv(
    EXTERNAL_CANDIDATE_FILE
)

external_candidates["gene"] = (
    external_candidates["gene"]
    .astype(str)
    .str.strip()
)

external_candidates_final = external_candidates[
    external_candidates["gene"].isin(FINAL_SIX)
].copy()


print(
    "Genes found in external candidate table:",
    len(external_candidates_final)
)


# ============================================================
# 11. LOAD EXTERNAL VALIDATED GENE LIST
# ============================================================

print("\n===== LOADING EXTERNAL VALIDATED GENE LIST =====")

with open(
    EXTERNAL_VALIDATED_LIST,
    "r"
) as f:

    validated_genes = {
        line.strip()
        for line in f
        if line.strip()
    }


print(
    "Genes in external validated list:",
    len(validated_genes)
)


# ============================================================
# 12. LOAD NOT FOUND LIST
# ============================================================

with open(
    EXTERNAL_NOT_FOUND,
    "r"
) as f:

    not_found_genes = {
        line.strip()
        for line in f
        if line.strip()
    }


print(
    "Genes in external not-found list:",
    len(not_found_genes)
)


# ============================================================
# 13. CREATE MASTER SIX-GENE TABLE
# ============================================================

print("\n===== BUILDING SIX-GENE EXTERNAL VALIDATION TABLE =====")

master = pd.DataFrame({
    "gene": FINAL_SIX
})


# ============================================================
# 14. STEP 41 EVIDENCE
# ============================================================

step41_columns = [
    "ML_folds_selected",
    "ML_consensus_status",
    "bioinformatics_final_50",
    "expression_mean_log2FC",
    "mean_malignant_detection_pct",
    "mean_diploid_detection_pct",
    "mean_detection_difference_pct",
    "expression_consistency",
    "expression_validation_score",
    "cross_validation_class"
]

available_step41_columns = [
    c for c in step41_columns
    if c in step41_final.columns
]


if available_step41_columns:

    step41_small = step41_final[
        ["gene"] + available_step41_columns
    ].copy()

    master = master.merge(
        step41_small,
        on="gene",
        how="left"
    )


# ============================================================
# 15. STEP 39 EVIDENCE
# ============================================================

step39_columns = [
    "final_evidence_rank",
    "final_evidence_category",
    "overall_exploratory_evidence_score",
    "external_support"
]

available_step39_columns = [
    c for c in step39_columns
    if c in step39_final.columns
]


if available_step39_columns:

    step39_small = step39_final[
        ["gene"] + available_step39_columns
    ].copy()

    master = master.merge(
        step39_small,
        on="gene",
        how="left"
    )


# ============================================================
# 16. EXTERNAL CANDIDATE PRESENCE
# ============================================================

master[
    "external_gene_present"
] = master["gene"].isin(
    set(
        external_candidates["gene"]
    )
)


# ============================================================
# 17. EXTERNAL VALIDATED LIST
# ============================================================

master[
    "external_validated_gene"
] = master["gene"].isin(
    validated_genes
)


# ============================================================
# 18. EXTERNAL NOT FOUND
# ============================================================

master[
    "external_gene_not_found"
] = master["gene"].isin(
    not_found_genes
)


# ============================================================
# 19. MERGE EXTERNAL EXPRESSION DATA
# ============================================================

external_expression_columns = [
    "external_cells_expressing",
    "external_expression_prevalence",
    "external_mean_expression",
    "samples_tested",
    "samples_detected_25pct",
    "samples_detected_50pct",
    "minimum_sample_prevalence",
    "mean_sample_prevalence",
    "maximum_sample_prevalence",
    "external_validation_status"
]

available_external_expression_columns = [
    c for c in external_expression_columns
    if c in external_expression_final.columns
]


if available_external_expression_columns:

    ext_small = external_expression_final[
        ["gene"] + available_external_expression_columns
    ].copy()

    master = master.merge(
        ext_small,
        on="gene",
        how="left"
    )


# ============================================================
# 20. MERGE SAMPLE CONSISTENCY
# ============================================================

sample_columns = [
    "samples_tested",
    "samples_detected_25pct",
    "samples_detected_50pct",
    "minimum_sample_prevalence",
    "mean_sample_prevalence",
    "maximum_sample_prevalence"
]

available_sample_columns = [
    c for c in sample_columns
    if c in external_sample_final.columns
]


if available_sample_columns:

    sample_small = external_sample_final[
        ["gene"] + available_sample_columns
    ].copy()

    sample_small = sample_small.rename(
        columns={
            c: f"external_{c}"
            for c in available_sample_columns
        }
    )

    master = master.merge(
        sample_small,
        on="gene",
        how="left"
    )


# ============================================================
# 21. DEFINE EXTERNAL SUPPORT CATEGORY
# ============================================================

def external_category(row):

    status = str(
        row.get(
            "external_validation_status",
            ""
        )
    ).lower()

    validated = bool(
        row.get(
            "external_validated_gene",
            False
        )
    )

    present = bool(
        row.get(
            "external_gene_present",
            False
        )
    )

    not_found = bool(
        row.get(
            "external_gene_not_found",
            False
        )
    )


    # Strong external expression support
    if (
        "strong_external_expression_support"
        in status
    ):
        return "Strong_external_expression_support"


    # Existing external support
    if validated and present:

        return "External_expression_supported"


    # Low prevalence
    if (
        "low_prevalence"
        in status
    ):
        return "Low_external_expression_support"


    # Present but not enough evidence
    if present:

        return "External_gene_present_only"


    # Explicitly not found
    if not_found:

        return "Not_found_externally"


    return "No_external_result"


master[
    "external_support_category"
] = master.apply(
    external_category,
    axis=1
)


# ============================================================
# 22. EXTERNAL VALIDATION SCORE
# ============================================================

def external_score(row):

    category = row[
        "external_support_category"
    ]

    if category == (
        "Strong_external_expression_support"
    ):
        return 2

    if category == (
        "External_expression_supported"
    ):
        return 2

    if category == (
        "Low_external_expression_support"
    ):
        return 1

    if category == (
        "External_gene_present_only"
    ):
        return 1

    return 0


master[
    "external_validation_score"
] = master.apply(
    external_score,
    axis=1
)


# ============================================================
# 23. MULTI-LAYER VALIDATION SCORE
# ============================================================

master[
    "ML_supported"
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
    "discovery_expression_supported"
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
    "bioinformatics_supported"
] = (
    master[
        "bioinformatics_final_50"
    ]
    .fillna(False)
    .astype(bool)
)


master[
    "external_supported"
] = (
    master[
        "external_validation_score"
    ] >= 2
)


master[
    "evidence_layers"
] = (

    master[
        "bioinformatics_supported"
    ].astype(int)

    +

    master[
        "ML_supported"
    ].astype(int)

    +

    master[
        "discovery_expression_supported"
    ].astype(int)

    +

    master[
        "external_supported"
    ].astype(int)

)


# ============================================================
# 24. FINAL EXTERNAL VALIDATION CLASS
# ============================================================

def final_external_class(row):

    layers = row[
        "evidence_layers"
    ]

    external = row[
        "external_support_category"
    ]

    if (
        layers >= 4
        and external
        in [
            "Strong_external_expression_support",
            "External_expression_supported"
        ]
    ):
        return (
            "Strong multi-layer externally supported candidate"
        )


    if (
        layers >= 3
        and external
        in [
            "Strong_external_expression_support",
            "External_expression_supported"
        ]
    ):
        return (
            "Externally supported consensus candidate"
        )


    if external == (
        "Low_external_expression_support"
    ):
        return (
            "ML consensus - weak external expression"
        )


    if external == (
        "External_gene_present_only"
    ):
        return (
            "ML consensus - external presence only"
        )


    if external == (
        "Not_found_externally"
    ):
        return (
            "ML consensus - not externally detected"
        )


    return (
        "Requires further external validation"
    )


master[
    "final_external_validation_class"
] = master.apply(
    final_external_class,
    axis=1
)


# ============================================================
# 25. SORT BY EVIDENCE
# ============================================================

master = master.sort_values(
    [
        "evidence_layers",
        "external_validation_score"
    ],
    ascending=[
        False,
        False
    ]
).reset_index(
    drop=True
)


master[
    "external_validation_rank"
] = (
    master.index + 1
)


# ============================================================
# 26. SAVE MASTER TABLE
# ============================================================

master_file = os.path.join(
    OUTPUT_DIR,
    "42_six_gene_external_validation.csv"
)

master.to_csv(
    master_file,
    index=False
)


# ============================================================
# 27. PRINT FINAL RESULTS
# ============================================================

print("\n" + "=" * 70)

print(
    "STEP 42 - SIX-GENE EXTERNAL VALIDATION RESULTS"
)

print("=" * 70)

display_columns = [
    "external_validation_rank",
    "gene",
    "ML_supported",
    "bioinformatics_supported",
    "discovery_expression_supported",
    "external_gene_present",
    "external_validated_gene",
    "external_support_category",
    "external_validation_score",
    "evidence_layers",
    "final_external_validation_class"
]

display_columns = [
    c for c in display_columns
    if c in master.columns
]

print(
    master[
        display_columns
    ].to_string(
        index=False
    )
)


# ============================================================
# 28. SUMMARY
# ============================================================

strong_external = master[
    "external_supported"
].sum()

weak_external = (
    master[
        "external_support_category"
    ]
    .eq(
        "Low_external_expression_support"
    )
    .sum()
)

not_supported = (
    master[
        "external_support_category"
    ]
    .isin(
        [
            "Not_found_externally",
            "No_external_result"
        ]
    )
    .sum()
)


print(
    "\n===== STEP 42 SUMMARY ====="
)

print(
    "Final ML consensus genes:",
    len(master)
)

print(
    "Externally supported:",
    int(strong_external)
)

print(
    "Low external expression support:",
    int(weak_external)
)

print(
    "Not externally supported:",
    int(not_supported)
)


print(
    "\nExternal validation categories:"
)

print(
    master[
        "external_support_category"
    ]
    .value_counts()
    .to_string()
)


# ============================================================
# 29. SAVE SUMMARY TEXT
# ============================================================

summary_file = os.path.join(
    OUTPUT_DIR,
    "42_external_validation_summary.txt"
)


with open(
    summary_file,
    "w"
) as f:

    f.write(
        "STEP 42 - ML CONSENSUS EXTERNAL VALIDATION\n"
    )

    f.write(
        "=" * 70 + "\n\n"
    )

    f.write(
        "External dataset: GSE176078\n"
    )

    f.write(
        "Validation type: Existing external expression "
        "validation integration\n\n"
    )

    f.write(
        "Final ML consensus genes:\n"
    )

    for gene in FINAL_SIX:

        f.write(
            f"- {gene}\n"
        )

    f.write(
        "\nRESULTS\n"
    )

    f.write(
        "-" * 70 + "\n"
    )

    for _, row in master.iterrows():

        f.write(
            f"\nGene: {row['gene']}\n"
        )

        f.write(
            f"ML supported: "
            f"{row['ML_supported']}\n"
        )

        f.write(
            f"Bioinformatics supported: "
            f"{row['bioinformatics_supported']}\n"
        )

        f.write(
            f"Discovery expression supported: "
            f"{row['discovery_expression_supported']}\n"
        )

        f.write(
            f"External gene present: "
            f"{row['external_gene_present']}\n"
        )

        f.write(
            f"External support category: "
            f"{row['external_support_category']}\n"
        )

        f.write(
            f"External validation score: "
            f"{row['external_validation_score']}\n"
        )

        f.write(
            f"Evidence layers: "
            f"{row['evidence_layers']}\n"
        )

        f.write(
            f"Final classification: "
            f"{row['final_external_validation_class']}\n"
        )


    f.write(
        "\n\nINTERPRETATION\n"
    )

    f.write(
        "-" * 70 + "\n"
    )

    f.write(
        "This analysis integrates the ML consensus panel "
        "with the previously completed GSE176078 external "
        "validation analysis.\n\n"
    )

    f.write(
        "External expression support indicates reproducible "
        "expression in an independent breast-cancer dataset. "
        "It does not establish clinical biomarker validity, "
        "causality, diagnostic performance, or prognostic utility.\n"
    )


# ============================================================
# 30. FIGURE 1
# EXTERNAL VALIDATION STATUS
# ============================================================

print(
    "\n===== GENERATING FIGURES ====="
)


category_counts = (
    master[
        "external_support_category"
    ]
    .value_counts()
)


fig, ax = plt.subplots(
    figsize=(9, 6)
)


ax.bar(
    category_counts.index,
    category_counts.values
)


ax.set_ylabel(
    "Number of genes"
)

ax.set_xlabel(
    "External validation category"
)

ax.set_title(
    "External Validation of ML Consensus Genes"
)

plt.xticks(
    rotation=35,
    ha="right"
)

plt.tight_layout()


figure1 = os.path.join(
    FIGURE_DIR,
    "42_external_validation_status.png"
)


plt.savefig(
    figure1,
    dpi=300,
    bbox_inches="tight"
)

plt.close()


# ============================================================
# 31. FIGURE 2
# DISCOVERY VS EXTERNAL PREVALENCE
# ============================================================

if (
    "expression_mean_log2FC"
    in master.columns
    and
    "external_expression_prevalence"
    in master.columns
):

    fig, ax = plt.subplots(
        figsize=(9, 6)
    )

    genes = master[
        "gene"
    ]

    external_prev = (
        pd.to_numeric(
            master[
                "external_expression_prevalence"
            ],
            errors="coerce"
        )
        * 100
    )

    discovery_prev = pd.to_numeric(
        master[
            "mean_malignant_detection_pct"
        ],
        errors="coerce"
    )

    x = np.arange(
        len(genes)
    )

    width = 0.36

    ax.bar(
        x - width / 2,
        discovery_prev,
        width,
        label="Discovery malignant detection"
    )

    ax.bar(
        x + width / 2,
        external_prev,
        width,
        label="GSE176078 expression prevalence"
    )

    ax.set_xticks(
        x
    )

    ax.set_xticklabels(
        genes
    )

    ax.set_ylabel(
        "Expression / detection (%)"
    )

    ax.set_xlabel(
        "Gene"
    )

    ax.set_title(
        "Discovery vs External Expression Prevalence"
    )

    ax.legend()

    plt.tight_layout()

    figure2 = os.path.join(
        FIGURE_DIR,
        "42_discovery_vs_external_expression.png"
    )

    plt.savefig(
        figure2,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

else:

    figure2 = None


# ============================================================
# 32. COMPLETION
# ============================================================

print(
    "\n" + "=" * 70
)

print(
    "STEP 42 COMPLETED"
)

print(
    "=" * 70
)

print(
    "\nMain output:"
)

print(
    master_file
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

if figure2:

    print(
        figure2
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