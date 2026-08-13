# ============================================================
# STEP 34
# LOGISTIC REGRESSION - 986 BIOINFORMATICS CANDIDATES
# Single-Cell Breast Cancer AI/ML Project
# ============================================================

import os
import numpy as np
import pandas as pd

from scipy.io import mmread

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    confusion_matrix
)


# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

DATA_DIR = "results/ai_ml/28_ml_dataset"

FEATURE_DIR = "results/ai_ml/31_feature_sets"

OUTPUT_DIR = (
    "results/ai_ml/34_logistic_regression_986"
)

os.makedirs(OUTPUT_DIR, exist_ok=True)


print("=" * 60)
print("STEP 34 - LOGISTIC REGRESSION")
print("986 BIOINFORMATICS CANDIDATES")
print("=" * 60)


# ------------------------------------------------------------
# Load expression matrix
# ------------------------------------------------------------

print("\nLoading expression matrix...")

X_all = mmread(
    os.path.join(
        DATA_DIR,
        "28_expression_matrix.mtx.gz"
    )
).tocsr()


# ------------------------------------------------------------
# Load genes
# ------------------------------------------------------------

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
    .values
)


# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

meta = pd.read_csv(
    os.path.join(
        DATA_DIR,
        "28_ml_metadata.csv"
    ),
    index_col=0
)


# ------------------------------------------------------------
# Dimension checks
# ------------------------------------------------------------

assert X_all.shape[0] == len(gene_names)
assert X_all.shape[1] == len(meta)


print(
    "\nExpression matrix:",
    X_all.shape
)

print(
    "Metadata:",
    meta.shape
)


# ------------------------------------------------------------
# Load 986 candidate genes
# ------------------------------------------------------------

feature_file = os.path.join(
    FEATURE_DIR,
    "31_feature_set_B_bioinformatics_candidates.tsv"
)

features = pd.read_csv(
    feature_file,
    sep="\t"
)

feature_genes = (
    features["gene"]
    .astype(str)
    .str.strip()
    .tolist()
)


print(
    "\nBioinformatics feature count:",
    len(feature_genes)
)


# ------------------------------------------------------------
# Verify candidate count
# ------------------------------------------------------------

if len(feature_genes) != 986:

    raise ValueError(
        f"Expected 986 candidates, "
        f"found {len(feature_genes)}"
    )


# ------------------------------------------------------------
# Map gene names to matrix indices
# ------------------------------------------------------------

gene_to_index = {
    gene: i
    for i, gene in enumerate(gene_names)
}


missing_genes = [
    gene
    for gene in feature_genes
    if gene not in gene_to_index
]


if len(missing_genes) > 0:

    raise ValueError(
        "Missing candidate genes: "
        + ", ".join(missing_genes)
    )


feature_indices = [
    gene_to_index[gene]
    for gene in feature_genes
]


# ------------------------------------------------------------
# Extract feature matrix
# ------------------------------------------------------------

X = X_all[
    feature_indices,
    :
].T.tocsr()


print(
    "ML feature matrix:",
    X.shape
)


# ------------------------------------------------------------
# Target
# ------------------------------------------------------------

y = meta[
    "ML_label"
].astype(int).values


groups = meta[
    "patient_id"
].astype(str).values


assert X.shape[0] == len(y)
assert X.shape[0] == len(groups)


# ------------------------------------------------------------
# LOPO evaluation
# ------------------------------------------------------------

samples = [
    "BC11",
    "BC12",
    "BC17"
]


results = []

coefficient_results = []


for test_sample in samples:

    print("\n" + "-" * 60)

    print(
        "TEST SAMPLE:",
        test_sample
    )


    train_mask = groups != test_sample

    test_mask = groups == test_sample


    X_train = X[train_mask]

    X_test = X[test_mask]

    y_train = y[train_mask]

    y_test = y[test_mask]


    print(
        "Training cells:",
        X_train.shape[0]
    )

    print(
        "Test cells:",
        X_test.shape[0]
    )

    print(
        "Training malignant:",
        np.sum(y_train == 1)
    )

    print(
        "Training diploid:",
        np.sum(y_train == 0)
    )

    print(
        "Test malignant:",
        np.sum(y_test == 1)
    )

    print(
        "Test diploid:",
        np.sum(y_test == 0)
    )


    # --------------------------------------------------------
    # Standardization
    # --------------------------------------------------------

    scaler = StandardScaler(
        with_mean=False
    )


    X_train_scaled = scaler.fit_transform(
        X_train
    )

    X_test_scaled = scaler.transform(
        X_test
    )


    # --------------------------------------------------------
    # Logistic Regression
    # --------------------------------------------------------

    model = LogisticRegression(
        C=1.0,
        solver="liblinear",
        class_weight="balanced",
        max_iter=2000,
        random_state=42
    )


    model.fit(
        X_train_scaled,
        y_train
    )


    # --------------------------------------------------------
    # Predictions
    # --------------------------------------------------------

    y_probability = model.predict_proba(
        X_test_scaled
    )[:, 1]
    
    print(
    "Test probability range:",
    f"{y_probability.min():.6f}",
    "to",
    f"{y_probability.max():.6f}"
    )

    print(
        "Mean probability:",
        f"{y_probability.mean():.6f}"
    )

    print(
        "Predicted malignant cells:",
        int(np.sum(y_probability >= 0.5))
    )

    print(
        "Predicted diploid cells:",
        int(np.sum(y_probability < 0.5))
    )

    print(
        "Mean probability - true malignant:",
        f"{y_probability[y_test == 1].mean():.6f}"
    )

    print(
        "Mean probability - true diploid:",
        f"{y_probability[y_test == 0].mean():.6f}"
    )


    y_pred = (
        y_probability >= 0.5
    ).astype(int)


    # --------------------------------------------------------
    # Metrics
    # --------------------------------------------------------

    auc = roc_auc_score(
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


    tn, fp, fn, tp = confusion_matrix(
        y_test,
        y_pred,
        labels=[0, 1]
    ).ravel()


    specificity = (
        tn / (tn + fp)
        if (tn + fp) > 0
        else np.nan
    )


    print(
        f"ROC-AUC: {auc:.4f}"
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


    results.append({

        "test_sample":
            test_sample,

        "n_test_cells":
            len(y_test),

        "malignant_test_cells":
            int(np.sum(y_test == 1)),

        "diploid_test_cells":
            int(np.sum(y_test == 0)),

        "roc_auc":
            auc,

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
            f1
    })


    # --------------------------------------------------------
    # Coefficients
    # --------------------------------------------------------

    coefficients = model.coef_[0]

    fold_coefficients = pd.DataFrame({

        "gene":
            feature_genes,

        "coefficient":
            coefficients,

        "test_sample":
            test_sample

    })


    coefficient_results.append(
        fold_coefficients
    )


# ------------------------------------------------------------
# Save fold results
# ------------------------------------------------------------

results_df = pd.DataFrame(
    results
)


results_df.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "34_logistic_regression_986_LOPO_results.csv"
    ),
    index=False
)


# ------------------------------------------------------------
# Save coefficients
# ------------------------------------------------------------

coef_df = pd.concat(
    coefficient_results,
    ignore_index=True
)


coef_df.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "34_logistic_regression_986_coefficients.csv"
    ),
    index=False
)


# ------------------------------------------------------------
# Summary statistics
# ------------------------------------------------------------

mean_metrics = pd.DataFrame({

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


mean_metrics.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "34_logistic_regression_986_summary.csv"
    ),
    index=False
)


# ------------------------------------------------------------
# Final
# ------------------------------------------------------------

print("\n" + "=" * 60)
print("STEP 34 COMPLETED")
print("=" * 60)

print("\nMean LOPO performance:")

print(
    mean_metrics.to_string(
        index=False
    )
)

print("\nResults saved to:")

print(OUTPUT_DIR)