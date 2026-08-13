# ============================================================
# STEP 35
# DATA-DRIVEN FEATURE SELECTION
#
# Starting feature space:
#   13,732 expression-filtered genes
#
# Outer validation:
#   Leave-one-sample-out using BC11, BC12, BC17
#
# Inner validation:
#   StratifiedKFold on training cells
#
# Feature selection:
#   SelectKBest + ANOVA F-test
#
# Model:
#   Logistic Regression
#
# IMPORTANT:
# The outer test sample is NEVER used for feature selection.
# ============================================================

import os
import warnings
import numpy as np
import pandas as pd

from scipy.io import mmread

from sklearn.feature_selection import (
    SelectKBest,
    f_classif
)

from sklearn.preprocessing import StandardScaler

from sklearn.linear_model import LogisticRegression

from sklearn.model_selection import StratifiedKFold

from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    confusion_matrix
)


# ============================================================
# PATHS
# ============================================================

DATA_DIR = (
    "results/ai_ml/28_ml_dataset"
)

FILTER_DIR = (
    "results/ai_ml/30_feature_selection"
)

OUTPUT_DIR = (
    "results/ai_ml/35_data_driven_features"
)

os.makedirs(
    OUTPUT_DIR,
    exist_ok=True
)


print("=" * 60)
print("STEP 35 - DATA-DRIVEN FEATURE SELECTION")
print("=" * 60)


# ============================================================
# LOAD EXPRESSION MATRIX
# ============================================================

print("\nLoading expression matrix...")

X_all = mmread(
    os.path.join(
        DATA_DIR,
        "28_expression_matrix.mtx.gz"
    )
).tocsr()


# ============================================================
# LOAD GENES
# ============================================================

genes = pd.read_csv(
    os.path.join(
        DATA_DIR,
        "28_genes.tsv"
    ),
    sep="\t",
    header=None,
    names=["gene"]
)

gene_names = (
    genes["gene"]
    .astype(str)
    .str.strip()
    .values
)


# ============================================================
# LOAD METADATA
# ============================================================

meta = pd.read_csv(
    os.path.join(
        DATA_DIR,
        "28_ml_metadata.csv"
    ),
    index_col=0
)


# ============================================================
# DIMENSION CHECKS
# ============================================================

assert (
    X_all.shape[0]
    == len(gene_names)
)

assert (
    X_all.shape[1]
    == len(meta)
)


print(
    "\nOriginal expression matrix:",
    X_all.shape
)

print(
    "Metadata:",
    meta.shape
)


# ============================================================
# LOAD STEP 30 FILTERED GENES
# ============================================================

filtered_gene_file = os.path.join(
    FILTER_DIR,
    "30_filtered_genes.tsv"
)


if not os.path.exists(
    filtered_gene_file
):

    raise FileNotFoundError(
        f"File not found: "
        f"{filtered_gene_file}"
    )


# No header assumption
filtered = pd.read_csv(
    filtered_gene_file,
    sep="\t",
    header=None
)


filtered_genes = (
    filtered.iloc[:, 0]
    .astype(str)
    .str.strip()
    .tolist()
)


# Remove possible accidental header
filtered_genes = [
    g
    for g in filtered_genes
    if g.lower() != "gene"
]


# Remove blanks
filtered_genes = [
    g
    for g in filtered_genes
    if g != ""
]


# Remove duplicates
filtered_genes = list(
    dict.fromkeys(
        filtered_genes
    )
)


print(
    "\nExpression-filtered genes:",
    len(filtered_genes)
)


if len(filtered_genes) != 13732:

    raise ValueError(
        "Expected 13,732 filtered genes, "
        f"found {len(filtered_genes)}"
    )


# ============================================================
# MAP GENES
# ============================================================

gene_to_index = {
    gene: i
    for i, gene
    in enumerate(gene_names)
}


missing_genes = [
    gene
    for gene in filtered_genes
    if gene not in gene_to_index
]


if missing_genes:

    raise ValueError(
        "Filtered genes missing from "
        "expression matrix: "
        + ", ".join(
            missing_genes[:20]
        )
    )


filtered_indices = [
    gene_to_index[gene]
    for gene in filtered_genes
]


# ============================================================
# CREATE 13,732-GENE MATRIX
# ============================================================

X = X_all[
    filtered_indices,
    :
].T.tocsr()


print(
    "Filtered ML matrix:",
    X.shape
)


# ============================================================
# TARGET AND SAMPLE INFORMATION
# ============================================================

y = (
    meta["ML_label"]
    .astype(int)
    .values
)


groups = (
    meta["patient_id"]
    .astype(str)
    .values
)


assert (
    X.shape[0]
    == len(y)
)

assert (
    X.shape[0]
    == len(groups)
)


# ============================================================
# CHECK CLASS DISTRIBUTION
# ============================================================

print(
    "\n===== CLASS DISTRIBUTION ====="
)

print(
    meta["ML_class"]
    .value_counts()
)


print(
    "\n===== SAMPLE DISTRIBUTION ====="
)

sample_table = pd.crosstab(
    groups,
    meta["ML_class"]
)

print(
    sample_table
)


# ============================================================
# OUTER TEST SAMPLES
# ============================================================

outer_samples = [
    "BC11",
    "BC12",
    "BC17"
]


# ============================================================
# FEATURE NUMBERS
# ============================================================

K_values = [
    50,
    100,
    250,
    500
]


# Fixed regularization for this experiment.
#
# We are primarily evaluating the effect of
# data-driven feature selection.
#
C_value = 0.1


print(
    "\nCandidate feature sizes:",
    K_values
)

print(
    "Logistic Regression C:",
    C_value
)


# ============================================================
# RESULT CONTAINERS
# ============================================================

outer_results = []

selected_features_all = []

inner_results_all = []


# ============================================================
# OUTER LOPO
# ============================================================

for test_sample in outer_samples:

    print("\n" + "-" * 60)

    print(
        "OUTER TEST SAMPLE:",
        test_sample
    )


    # --------------------------------------------------------
    # Outer split
    # --------------------------------------------------------

    train_mask = (
        groups != test_sample
    )

    test_mask = (
        groups == test_sample
    )


    X_train = X[
        train_mask
    ]

    X_test = X[
        test_mask
    ]

    y_train = y[
        train_mask
    ]

    y_test = y[
        test_mask
    ]


    train_groups = groups[
        train_mask
    ]


    print(
        "Training cells:",
        len(y_train)
    )

    print(
        "Test cells:",
        len(y_test)
    )

    print(
        "Training samples:",
        sorted(
            np.unique(
                train_groups
            )
        )
    )

    print(
        "Test malignant:",
        int(
            np.sum(
                y_test == 1
            )
        )
    )

    print(
        "Test diploid:",
        int(
            np.sum(
                y_test == 0
            )
        )
    )


    # ========================================================
    # INNER STRATIFIED CV
    # ========================================================

    inner_cv = StratifiedKFold(
        n_splits=5,
        shuffle=True,
        random_state=42
    )


    inner_scores = []


    for K in K_values:

        fold_scores = []


        print(
            f"\nInner CV: K={K}"
        )


        for (
            inner_train_idx,
            inner_val_idx
        ) in inner_cv.split(
            X_train,
            y_train
        ):


            X_inner_train = X_train[
                inner_train_idx
            ]

            X_inner_val = X_train[
                inner_val_idx
            ]


            y_inner_train = y_train[
                inner_train_idx
            ]

            y_inner_val = y_train[
                inner_val_idx
            ]


            # ------------------------------------------------
            # FEATURE SELECTION
            # INNER TRAIN ONLY
            # ------------------------------------------------

            selector = SelectKBest(
                score_func=f_classif,
                k=K
            )


            with warnings.catch_warnings():

                warnings.simplefilter(
                    "ignore"
                )

                X_inner_train_selected = (
                    selector.fit_transform(
                        X_inner_train,
                        y_inner_train
                    )
                )


            X_inner_val_selected = (
                selector.transform(
                    X_inner_val
                )
            )


            # ------------------------------------------------
            # SCALING
            # INNER TRAIN ONLY
            # ------------------------------------------------

            scaler = StandardScaler(
                with_mean=False
            )


            X_inner_train_scaled = (
                scaler.fit_transform(
                    X_inner_train_selected
                )
            )


            X_inner_val_scaled = (
                scaler.transform(
                    X_inner_val_selected
                )
            )


            # ------------------------------------------------
            # LOGISTIC REGRESSION
            # ------------------------------------------------

            model = LogisticRegression(

                C=C_value,

                solver="liblinear",

                class_weight="balanced",

                max_iter=2000,

                random_state=42

            )


            model.fit(
                X_inner_train_scaled,
                y_inner_train
            )


            probabilities = (
                model.predict_proba(
                    X_inner_val_scaled
                )[:, 1]
            )


            auc = roc_auc_score(
                y_inner_val,
                probabilities
            )


            fold_scores.append(
                auc
            )


        mean_auc = np.mean(
            fold_scores
        )


        print(
            f"Mean inner ROC-AUC: "
            f"{mean_auc:.4f}"
        )


        inner_results_all.append({

            "test_sample":
                test_sample,

            "K":
                K,

            "inner_auc_mean":
                mean_auc,

            "inner_auc_std":
                np.std(
                    fold_scores
                )

        })


        inner_scores.append({
            "K": K,
            "mean_auc": mean_auc
        })


    # ========================================================
    # SELECT BEST K
    # ========================================================

    inner_df = pd.DataFrame(
        inner_scores
    )


    best_row = (
        inner_df
        .sort_values(
            "mean_auc",
            ascending=False
        )
        .iloc[0]
    )


    best_K = int(
        best_row["K"]
    )


    best_inner_auc = float(
        best_row["mean_auc"]
    )


    print(
        "\nSelected K:",
        best_K
    )

    print(
        "Best inner ROC-AUC:",
        f"{best_inner_auc:.4f}"
    )


    # ========================================================
    # FINAL FEATURE SELECTION
    # OUTER TRAINING DATA ONLY
    # ========================================================

    selector = SelectKBest(
        score_func=f_classif,
        k=best_K
    )


    with warnings.catch_warnings():

        warnings.simplefilter(
            "ignore"
        )

        X_train_selected = (
            selector.fit_transform(
                X_train,
                y_train
            )
        )


    X_test_selected = (
        selector.transform(
            X_test
        )
    )


    # ========================================================
    # SELECTED GENE NAMES
    # ========================================================

    selected_mask = (
        selector.get_support()
    )


    selected_genes = [

        filtered_genes[i]

        for i, selected

        in enumerate(
            selected_mask
        )

        if selected

    ]


    print(
        "Selected genes:",
        len(selected_genes)
    )


    # Save selected genes
    for rank, gene in enumerate(
        selected_genes,
        start=1
    ):

        selected_features_all.append({

            "test_sample":
                test_sample,

            "rank":
                rank,

            "gene":
                gene,

            "K":
                best_K

        })


    # ========================================================
    # SCALE
    # OUTER TRAINING DATA ONLY
    # ========================================================

    scaler = StandardScaler(
        with_mean=False
    )


    X_train_scaled = (
        scaler.fit_transform(
            X_train_selected
        )
    )


    X_test_scaled = (
        scaler.transform(
            X_test_selected
        )
    )


    # ========================================================
    # FINAL OUTER MODEL
    # ========================================================

    model = LogisticRegression(

        C=C_value,

        solver="liblinear",

        class_weight="balanced",

        max_iter=2000,

        random_state=42

    )


    model.fit(
        X_train_scaled,
        y_train
    )


    # ========================================================
    # OUTER TEST PREDICTION
    # ========================================================

    y_probability = (
        model.predict_proba(
            X_test_scaled
        )[:, 1]
    )


    y_pred = (
        y_probability >= 0.5
    ).astype(int)


    # ========================================================
    # METRICS
    # ========================================================

    roc_auc = roc_auc_score(
        y_test,
        y_probability
    )


    pr_auc = average_precision_score(
        y_test,
        y_probability
    )


    accuracy = accuracy_score(
        y_test,
        y_pred
    )


    precision = precision_score(
        y_test,
        y_pred,
        zero_division=0
    )


    recall = recall_score(
        y_test,
        y_pred,
        zero_division=0
    )


    f1 = f1_score(
        y_test,
        y_pred,
        zero_division=0
    )


    tn, fp, fn, tp = (
        confusion_matrix(
            y_test,
            y_pred,
            labels=[
                0,
                1
            ]
        ).ravel()
    )


    specificity = (
        tn / (tn + fp)
        if (tn + fp) > 0
        else np.nan
    )


    # ========================================================
    # PRINT PERFORMANCE
    # ========================================================

    print(
        "\nOUTER TEST PERFORMANCE"
    )

    print(
        f"ROC-AUC: {roc_auc:.4f}"
    )

    print(
        f"PR-AUC: {pr_auc:.4f}"
    )

    print(
        f"Accuracy: {accuracy:.4f}"
    )

    print(
        f"Precision: {precision:.4f}"
    )

    print(
        f"Recall: {recall:.4f}"
    )

    print(
        f"Specificity: {specificity:.4f}"
    )

    print(
        f"F1: {f1:.4f}"
    )


    # ========================================================
    # SAVE OUTER RESULTS
    # ========================================================

    outer_results.append({

        "test_sample":
            test_sample,

        "selected_K":
            best_K,

        "inner_best_auc":
            best_inner_auc,

        "roc_auc":
            roc_auc,

        "pr_auc":
            pr_auc,

        "accuracy":
            accuracy,

        "precision":
            precision,

        "recall":
            recall,

        "specificity":
            specificity,

        "f1":
            f1,

        "n_test_cells":
            len(y_test),

        "malignant_test_cells":
            int(
                np.sum(
                    y_test == 1
                )
            ),

        "diploid_test_cells":
            int(
                np.sum(
                    y_test == 0
                )
            )

    })


# ============================================================
# SAVE OUTER RESULTS
# ============================================================

results_df = pd.DataFrame(
    outer_results
)


results_df.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "35_data_driven_LOPO_results.csv"
    ),
    index=False
)


# ============================================================
# SAVE INNER CV RESULTS
# ============================================================

inner_results_df = pd.DataFrame(
    inner_results_all
)


inner_results_df.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "35_inner_feature_selection_results.csv"
    ),
    index=False
)


# ============================================================
# SAVE SELECTED GENES
# ============================================================

selected_df = pd.DataFrame(
    selected_features_all
)


selected_df.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "35_selected_genes_by_outer_fold.csv"
    ),
    index=False
)


# ============================================================
# GENE RECURRENCE
# ============================================================

if len(selected_df) > 0:

    recurrence = (
        selected_df
        .groupby("gene")
        .agg(
            folds_selected=(
                "test_sample",
                "nunique"
            )
        )
        .reset_index()
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


    recurrence.to_csv(
        os.path.join(
            OUTPUT_DIR,
            "35_gene_selection_recurrence.csv"
        ),
        index=False
    )


# ============================================================
# SUMMARY
# ============================================================

summary = pd.DataFrame({

    "metric": [

        "roc_auc",
        "pr_auc",
        "accuracy",
        "precision",
        "recall",
        "specificity",
        "f1"

    ],

    "mean": [

        results_df["roc_auc"].mean(),

        results_df["pr_auc"].mean(),

        results_df["accuracy"].mean(),

        results_df["precision"].mean(),

        results_df["recall"].mean(),

        results_df["specificity"].mean(),

        results_df["f1"].mean()

    ],

    "std": [

        results_df["roc_auc"].std(),

        results_df["pr_auc"].std(),

        results_df["accuracy"].std(),

        results_df["precision"].std(),

        results_df["recall"].std(),

        results_df["specificity"].std(),

        results_df["f1"].std()

    ]

})


summary.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "35_data_driven_summary.csv"
    ),
    index=False
)


# ============================================================
# FINAL
# ============================================================

print("\n" + "=" * 60)
print("STEP 35 COMPLETED")
print("=" * 60)

print(
    "\nOuter LOPO results:"
)

print(
    results_df.to_string(
        index=False
    )
)


print(
    "\nMean performance:"
)

print(
    summary.to_string(
        index=False
    )
)


print(
    "\nResults saved to:"
)

print(
    OUTPUT_DIR
)