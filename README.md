# Single-Cell RNA-seq and Machine Learning for Breast Cancer Biomarker Discovery and Validation

<p align="center">

**A sample-aware computational workflow integrating single-cell RNA-seq bioinformatics, machine learning, external validation, and tumor-specificity analysis for breast cancer biomarker candidate prioritization**

</p>

---

## Overview

Breast cancer is a highly heterogeneous disease in which malignant and non-malignant cell populations can coexist within the same tumor microenvironment.

Single-cell RNA sequencing (scRNA-seq) provides an opportunity to study this heterogeneity at cellular resolution and identify genes associated with malignant cell populations.

In this project, I developed a **sample-aware single-cell RNA-seq and machine-learning workflow** to prioritize breast cancer biomarker candidates.

The workflow integrates:

- Single-cell RNA-seq analysis
- Cell-level biological characterization
- CNV-supported malignant-cell identification
- Bioinformatics candidate prioritization
- Expression-based feature filtering
- Machine-learning feature selection
- Logistic regression
- Leave-One-Sample-Out (LOPO) validation
- ML feature recurrence and consensus
- Discovery-dataset expression validation
- External validation using GSE176078
- Tumor/cell-type specificity analysis
- Multi-layer evidence integration

The analysis progressively reduced the candidate search space:

```text
36,626 genes
      ↓
13,732 expression-filtered genes
      ↓
986 bioinformatics candidates
      ↓
50 final bioinformatics candidates
      ↓
6 consensus candidates
