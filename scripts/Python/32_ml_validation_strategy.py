# ============================================================
# STEP 32
# ML VALIDATION STRATEGY
# Single-Cell Breast Cancer AI/ML Project
# ============================================================

import os
import pandas as pd


DATA_DIR = "results/ai_ml/28_ml_dataset"

OUTPUT_DIR = "results/ai_ml/32_validation"

os.makedirs(OUTPUT_DIR, exist_ok=True)


print("=" * 60)
print("STEP 32 - ML VALIDATION STRATEGY")
print("=" * 60)


# ------------------------------------------------------------
# Load metadata
# ------------------------------------------------------------

metadata_file = os.path.join(
    DATA_DIR,
    "28_ml_metadata.csv"
)

meta = pd.read_csv(
    metadata_file,
    index_col=0
)


# ------------------------------------------------------------
# Overall dataset
# ------------------------------------------------------------

print("\n===== OVERALL DATASET =====")

print(
    "Cells:",
    len(meta)
)

print(
    "Samples:",
    meta["patient_id"].nunique()
)

print("\nClass distribution:")

print(
    meta["ML_class"]
    .value_counts()
)


# ------------------------------------------------------------
# Identify samples containing malignant cells
# ------------------------------------------------------------

sample_class = pd.crosstab(
    meta["patient_id"],
    meta["ML_class"]
)

print("\n===== SAMPLE × CLASS =====")

print(sample_class)


malignant_samples = sample_class.index[
    sample_class.get(
        "malignant",
        pd.Series(0, index=sample_class.index)
    ) > 0
].tolist()


print("\nSamples containing malignant cells:")

for sample in malignant_samples:
    print("-", sample)


print(
    "\nNumber of informative samples:",
    len(malignant_samples)
)


# ------------------------------------------------------------
# Create LOPO strategy
# ------------------------------------------------------------

print("\n===== LEAVE-ONE-SAMPLE-OUT STRATEGY =====")

lopo_rows = []

for test_sample in malignant_samples:

    train_samples = [
        s for s in malignant_samples
        if s != test_sample
    ]

    lopo_rows.append({

        "test_sample": test_sample,

        "training_samples":
            ";".join(train_samples),

        "n_training_samples":
            len(train_samples)

    })

    print(
        f"Test: {test_sample} | "
        f"Train: {', '.join(train_samples)}"
    )


lopo = pd.DataFrame(lopo_rows)


# ------------------------------------------------------------
# Save strategy
# ------------------------------------------------------------

sample_class.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "32_sample_class_distribution.csv"
    )
)

lopo.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "32_leave_one_sample_out_strategy.csv"
    ),
    index=False
)


# ------------------------------------------------------------
# Validation report
# ------------------------------------------------------------

report = pd.DataFrame({

    "metric": [
        "total_cells",
        "total_samples",
        "malignant_cells",
        "diploid_cells",
        "informative_samples"
    ],

    "value": [
        len(meta),
        meta["patient_id"].nunique(),
        int(
            (meta["ML_label"] == 1).sum()
        ),
        int(
            (meta["ML_label"] == 0).sum()
        ),
        len(malignant_samples)
    ]

})


report.to_csv(
    os.path.join(
        OUTPUT_DIR,
        "32_validation_strategy_summary.csv"
    ),
    index=False
)


# ------------------------------------------------------------
# Final
# ------------------------------------------------------------

print("\n" + "=" * 60)
print("STEP 32 COMPLETED")
print("=" * 60)

print(
    "\nPrimary informative samples:",
    ", ".join(malignant_samples)
)

print(
    "\nResults saved to:"
)

print(OUTPUT_DIR)