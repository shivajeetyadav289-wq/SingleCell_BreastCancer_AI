# ============================================================
# STEP 31
# CONSTRUCT ML FEATURE SETS
# Single-Cell Breast Cancer AI/ML Project
# ============================================================

import os
import pandas as pd


# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

FILTERED_DIR = (
    "results/ai_ml/30_feature_selection"
)

BIOINFO_FILE = (
    "results/malignant/integrated_prioritization/"
    "24_integrated_malignant_biomarker_ranking.csv"
)

FINAL_PANEL_FILE = (
    "results/malignant/final_candidates/"
    "25_final_research_candidate_panel.csv"
)

OUTPUT_DIR = (
    "results/ai_ml/31_feature_sets"
)

os.makedirs(OUTPUT_DIR, exist_ok=True)


print("=" * 60)
print("STEP 31 - FEATURE SET CONSTRUCTION")
print("=" * 60)


# ------------------------------------------------------------
# Load expression-filtered genes
# ------------------------------------------------------------

filtered_gene_file = os.path.join(
    FILTERED_DIR,
    "30_filtered_genes.tsv"
)

filtered_genes = pd.read_csv(
    filtered_gene_file,
    sep="\t",
    header=None,
    names=["gene"]
)

filtered_genes["gene"] = (
    filtered_genes["gene"]
    .astype(str)
    .str.strip()
)

filtered_gene_set = set(
    filtered_genes["gene"]
)

print("\n===== EXPRESSION FILTERED GENES =====")

print(
    "Genes available:",
    len(filtered_genes)
)


# ------------------------------------------------------------
# Load 986 bioinformatics candidates
# ------------------------------------------------------------

print("\n===== BIOINFORMATICS CANDIDATES =====")

bioinfo = pd.read_csv(
    BIOINFO_FILE
)

print(
    "Rows in source file:",
    len(bioinfo)
)

print(
    "Columns:",
    len(bioinfo.columns)
)


# Verify gene column

if "gene" not in bioinfo.columns:

    raise ValueError(
        "Could not find 'gene' column in "
        "the 986 candidate file."
    )


bioinfo_genes = (
    bioinfo["gene"]
    .astype(str)
    .str.strip()
)

bioinfo_genes = (
    bioinfo_genes[
        bioinfo_genes.notna()
    ]
    .drop_duplicates()
)


print(
    "Unique bioinformatics candidates:",
    len(bioinfo_genes)
)


# ------------------------------------------------------------
# Intersect 986 candidates with expression features
# ------------------------------------------------------------

bioinfo_present = [
    g for g in bioinfo_genes
    if g in filtered_gene_set
]

bioinfo_absent = [
    g for g in bioinfo_genes
    if g not in filtered_gene_set
]


print(
    "Present after expression filtering:",
    len(bioinfo_present)
)

print(
    "Absent after expression filtering:",
    len(bioinfo_absent)
)


# ------------------------------------------------------------
# Load 50 final candidates
# ------------------------------------------------------------

print("\n===== FINAL 50 CANDIDATES =====")

final_panel = pd.read_csv(
    FINAL_PANEL_FILE
)

print(
    "Rows in source file:",
    len(final_panel)
)

if "gene" not in final_panel.columns:

    raise ValueError(
        "Could not find 'gene' column in "
        "the final 50 candidate file."
    )


final_genes = (
    final_panel["gene"]
    .astype(str)
    .str.strip()
)

final_genes = (
    final_genes[
        final_genes.notna()
    ]
    .drop_duplicates()
)


print(
    "Unique final candidates:",
    len(final_genes)
)


# ------------------------------------------------------------
# Intersect 50 candidates with expression features
# ------------------------------------------------------------

final_present = [
    g for g in final_genes
    if g in filtered_gene_set
]

final_absent = [
    g for g in final_genes
    if g not in filtered_gene_set
]


print(
    "Present after expression filtering:",
    len(final_present)
)

print(
    "Absent after expression filtering:",
    len(final_absent)
)


# ------------------------------------------------------------
# Verify final panel is contained in 986 candidates
# ------------------------------------------------------------

bioinfo_set = set(bioinfo_genes)

final_in_bioinfo = [
    g for g in final_genes
    if g in bioinfo_set
]

final_not_in_bioinfo = [
    g for g in final_genes
    if g not in bioinfo_set
]


print("\n===== CANDIDATE HIERARCHY CHECK =====")

print(
    "Final 50 also present in 986:",
    len(final_in_bioinfo)
)

print(
    "Final 50 not present in 986:",
    len(final_not_in_bioinfo)
)


# ------------------------------------------------------------
# Save feature sets
# ------------------------------------------------------------

# A. All expression-filtered genes

pd.DataFrame({
    "gene": filtered_genes["gene"]
}).to_csv(
    os.path.join(
        OUTPUT_DIR,
        "31_feature_set_A_all_filtered_genes.tsv"
    ),
    sep="\t",
    index=False
)


# B. 986 bioinformatics candidates

pd.DataFrame({
    "gene": bioinfo_present
}).to_csv(
    os.path.join(
        OUTPUT_DIR,
        "31_feature_set_B_bioinformatics_candidates.tsv"
    ),
    sep="\t",
    index=False
)


# C. 50 final research candidates

pd.DataFrame({
    "gene": final_present
}).to_csv(
    os.path.join(
        OUTPUT_DIR,
        "31_feature_set_C_final_50_candidates.tsv"
    ),
    sep="\t",
    index=False
)


# ------------------------------------------------------------
# Save absent candidate lists
# ------------------------------------------------------------

pd.DataFrame({
    "gene": bioinfo_absent
}).to_csv(
    os.path.join(
        OUTPUT_DIR,
        "31_bioinformatics_candidates_absent.tsv"
    ),
    sep="\t",
    index=False
)


pd.DataFrame({
    "gene": final_absent
}).to_csv(
    os.path.join(
        OUTPUT_DIR,
        "31_final_50_candidates_absent.tsv"
    ),
    sep="\t",
    index=False
)


# ------------------------------------------------------------
# Save hierarchy information
# ------------------------------------------------------------

pd.DataFrame({
    "gene": final_in_bioinfo
}).to_csv(
    os.path.join(
        OUTPUT_DIR,
        "31_final_50_in_986_candidates.tsv"
    ),
    sep="\t",
    index=False
)


pd.DataFrame({
    "gene": final_not_in_bioinfo
}).to_csv(
    os.path.join(
        OUTPUT_DIR,
        "31_final_50_not_in_986_candidates.tsv"
    ),
    sep="\t",
    index=False
)


# ------------------------------------------------------------
# Feature-set summary
# ------------------------------------------------------------

summary = pd.DataFrame({

    "feature_set": [
        "A_all_expression_filtered",
        "B_bioinformatics_986",
        "C_final_50"
    ],

    "source_genes": [
        len(filtered_genes),
        len(bioinfo_genes),
        len(final_genes)
    ],

    "genes_present_in_expression": [
        len(filtered_genes),
        len(bioinfo_present),
        len(final_present)
    ],

    "genes_absent_from_expression": [
        0,
        len(bioinfo_absent),
        len(final_absent)
    ]
})


summary.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "31_feature_set_summary.csv"
    ),
    index=False
)


# ------------------------------------------------------------
# Final
# ------------------------------------------------------------

print("\n" + "=" * 60)
print("STEP 31 COMPLETED")
print("=" * 60)

print("\nFeature sets:")

print(
    "A - All filtered genes:",
    len(filtered_genes)
)

print(
    "B - Bioinformatics candidates:",
    len(bioinfo_present)
)

print(
    "C - Final 50 candidates:",
    len(final_present)
)

print("\nResults saved to:")

print(OUTPUT_DIR)