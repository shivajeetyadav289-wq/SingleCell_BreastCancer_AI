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

The six final consensus candidates were:

COX8A, EPN3, ETFB, FCGRT, MRPS7 and PRSS23

Among these, PRSS23, FCGRT, ETFB and MRPS7 showed the strongest tumor-associated expression evidence in the discovery analysis

Important: These genes are computationally prioritized research candidates and should not be interpreted as clinically validated biomarkers.

Research Question

Can single-cell RNA-seq bioinformatics combined with machine learning identify breast cancer gene candidates that show reproducible association with malignant cells across independent samples?

The project was designed to evaluate candidates using multiple evidence layers rather than relying on a single statistical or machine-learning result

Bioinformatics evidence
          +
Machine-learning recurrence
          +
Discovery expression validation
          +
External validation
          +
Tumor specificity
          ↓
Final candidate prioritization

Objectives
Primary Objective

To develop a sample-aware computational framework for breast cancer biomarker candidate prioritization using single-cell RNA-seq and machine learning.

Specific Objectives
Characterize breast cancer single-cell RNA-seq data.
Identify malignant-cell candidates using CNV-supported evidence.
Prioritize biologically relevant genes using bioinformatics analysis.
Perform expression-based feature filtering.
Construct machine-learning feature sets.
Apply logistic regression for malignant-cell classification.
Evaluate cross-sample generalization using Leave-One-Sample-Out validation.
Identify recurrent ML-selected genes.
Integrate machine-learning and bioinformatics evidence.
Validate consensus candidates at the expression level.
Perform external validation using GSE176078.
Evaluate tumor/cell-type specificity.
Integrate multiple evidence layers into a final candidate panel.
Interpret the candidates while explicitly documenting limitations.

Dataset
Discovery Dataset

The primary discovery dataset used in this project was:

GSE228499

The processed single-cell object contained:
| Property    |    Value |
| ----------- | -------: |
| Cells       |   24,575 |
| Genes       |   36,626 |
| Samples     |        9 |
| Assay       |      RNA |
| Seurat      |    5.5.1 |
| Assay class | `Assay5` |

BC03
BC05
BC06
BC08
BC11
BC12
BC14
BC15
BC17
Malignant Cell Identification

Malignant-cell candidates were identified using integrated biological annotation and CNV-supported evidence.

The original processed object contained:
24,575 cells
36,626 genes

The machine-learning analysis subsequently focused on cells with defined ML labels.

The final ML dataset contained:

| ML class              |     Cells |
| --------------------- | --------: |
| Malignant             |     1,017 |
| Diploid/non-malignant |     2,362 |
| **Total**             | **3,379** |

Informative Samples

Only three samples contained malignant cells:

| Sample | Malignant | Diploid/non-malignant | Total |
| ------ | --------: | --------------------: | ----: |
| BC11   |       278 |                   183 |   461 |
| BC12   |       380 |                    32 |   412 |
| BC17   |       359 |                     2 |   361 |


Therefore, the primary sample-aware machine-learning validation was performed across:

BC11
BC12
BC17

This limited number of informative samples is an important limitation of the study.

Workflow

Major workflow
GSE228499
    ↓
Single-cell processing
    ↓
Cell annotation
    ↓
CNV-supported malignant-cell identification
    ↓
Bioinformatics candidate prioritization
    ↓
Expression filtering
    ↓
13,732 genes
    ↓
986 bioinformatics candidates
    ↓
50 final candidates
    ↓
Machine Learning
    ↓
Logistic Regression
    ↓
LOPO validation
    ↓
ML feature recurrence
    ↓
6 consensus candidates
    ↓
Expression validation
    ↓
GSE176078 external validation
    ↓
Tumor-specificity analysis
    ↓
Final evidence integration

Bioinformatics Candidate Prioritization

Before machine learning, genes were prioritized through the bioinformatics pipeline.

The analysis incorporated evidence including:

malignant-cell characterization
differential expression
prevalence
sample-level consistency
biological annotation
epithelial/tumor-associated evidence
mitochondrial/OXPHOS characterization
external validation
biological review

This produced:

986 bioinformatics candidates

which were subsequently narrowed to:

50 final bioinformatics candidates

Feature Filtering

Expression-based feature filtering was performed before machine learning.

Filtering results:
| Stage                     |  Genes |
| ------------------------- | -----: |
| Original genes            | 36,626 |
| Expression-filtered genes | 13,732 |
| Removed                   | 22,894 |

iltering criteria:

Detection threshold: 1%
Minimum detected cells: 34

Final expression feature space:

13,732 genes × 3,379 cells

Machine Learning Strategy

The primary machine-learning model was:

Logistic Regression

Three feature levels were evaluated:
| Feature set                     | Number of genes |
| ------------------------------- | --------------: |
| All expression-filtered genes   |          13,732 |
| Bioinformatics candidates       |             986 |
| Final bioinformatics candidates |              50 |

The purpose of the ML analysis was not simply to maximize predictive accuracy.

The main objective was to determine whether genes repeatedly emerged as useful features across independent sample-level validation folds.

Leave-One-Sample-Out Validation

Randomly splitting individual cells can result in overly optimistic estimates because cells from the same patient/sample are biologically related.

Therefore, the project used Leave-One-Sample-Out (LOPO) validation.

The strategy was:
Test BC11
    ↓
Train using remaining samples

Test BC12
    ↓
Train using remaining samples

Test BC17
    ↓
Train using remaining samples
This approach asks:

Can the model generalize to a sample that was not used during training?

This is more appropriate for evaluating sample-level generalization than a simple random cell-level split.

50-Gene Logistic Regression Results

The initial 50-gene model produced the following mean LOPO performance:
| Metric      |  Mean |
| ----------- | ----: |
| ROC-AUC     | 0.824 |
| PR-AUC      | 0.902 |
| Accuracy    | 0.325 |
| Precision   | 0.913 |
| Recall      | 0.274 |
| Specificity | 0.909 |
| F1          | 0.368 |

The results showed useful ranking performance but relatively low recall and F1.

Therefore, the model was not treated as a clinically useful classifier.

Instead, the ML analysis was used as an additional feature-prioritization and recurrence layer.

986-Gene Machine Learning Analysis

The broader 986-gene feature set was also evaluated.

Mean LOPO performance:
| Metric      |  Mean |
| ----------- | ----: |
| ROC-AUC     | 0.241 |
| PR-AUC      | 0.765 |
| Accuracy    | 0.160 |
| Precision   | 0.000 |
| Recall      | 0.000 |
| Specificity | 1.000 |
| F1          | 0.000 |

Data-Driven Feature Selection

A second feature-selection analysis evaluated:

50 genes
100 genes
250 genes
500 genes

The inner validation results were very high.

However, outer-sample performance was substantially weaker.

Mean outer LOPO performance:

| Metric      |  Mean |
| ----------- | ----: |
| ROC-AUC     | 0.437 |
| PR-AUC      | 0.882 |
| Accuracy    | 0.167 |
| Recall      | 0.009 |
| Specificity | 1.000 |
| F1          | 0.017 |

The selected feature size was:

500 genes

This difference between internal validation and outer-sample performance is an important result.

It demonstrates why sample-aware validation is essential when working with single-cell data from multiple patients/samples.

ML Consensus Analysis

Genes were evaluated for recurrence across the three informative LOPO folds.

Results:
| Category                   | Number of genes |
| -------------------------- | --------------: |
| Final bioinformatics panel |              50 |
| ML-selected in all 3 folds |              80 |
| ML-selected in ≥2 folds    |             469 |
| Strong consensus           |               6 |
| Moderate consensus         |              33 |
| Broad consensus            |              39 |

The six strong consensus genes were:

COX8A
EPN3
ETFB
FCGRT
MRPS7
PRSS23

These six genes were present in the final bioinformatics candidate panel and were repeatedly selected during the ML analysis.

Final Six Consensus Candidates:
| Gene       | Interpretation                                           |
| ---------- | -------------------------------------------------------- |
| **PRSS23** | Strong tumor-associated research candidate               |
| **FCGRT**  | Strong tumor-associated research candidate               |
| **ETFB**   | Strong mitochondrial/metabolic research candidate        |
| **MRPS7**  | Strong mitochondrial/ribosomal research candidate        |
| **EPN3**   | Promising candidate requiring further validation         |
| **COX8A**  | Strong expression candidate with lower tumor specificity |

Discovery Expression Validation

The six consensus genes were evaluated across BC11, BC12 and BC17.
| Gene   | Mean log2FC | Malignant detection | Diploid detection | Difference |
| ------ | ----------: | ------------------: | ----------------: | ---------: |
| PRSS23 |       2.986 |              73.94% |            37.85% |   36.09 pp |
| FCGRT  |       2.119 |              70.80% |            29.17% |   41.63 pp |
| ETFB   |       2.372 |              65.09% |            21.46% |   43.63 pp |
| MRPS7  |       2.642 |              74.53% |            39.20% |   35.33 pp |
| COX8A  |       0.716 |              98.13% |            80.69% |   17.44 pp |
| EPN3   |       0.321 |              45.23% |             7.41% |   37.82 pp |

Expression Validation Scores:

| Gene   | Expression validation score |
| ------ | --------------------------: |
| MRPS7  |                       0.961 |
| PRSS23 |                       0.954 |
| ETFB   |                       0.885 |
| FCGRT  |                       0.824 |
| COX8A  |                       0.663 |
| EPN3   |                       0.333 |

MRPS7 and PRSS23 showed the strongest expression-validation scores.

EPN3 showed weaker sample consistency.

Tumor-Specificity Analysis

The six consensus candidates were further evaluated for malignant/tumor association.

PRSS23
Malignant detection: 73.94%
Diploid detection: 37.85%
Detection difference: 36.09 percentage points
Interpretation: Tumor-associated candidate
FCGRT
Malignant detection: 70.80%
Diploid detection: 29.17%
Detection difference: 41.63 percentage points
Interpretation: Tumor-associated candidate
ETFB
Malignant detection: 65.09%
Diploid detection: 21.46%
Detection difference: 43.63 percentage points
Interpretation: Tumor-associated candidate
MRPS7
Malignant detection: 74.53%
Diploid detection: 39.20%
Detection difference: 35.33 percentage points
Interpretation: Tumor-associated candidate
EPN3
Malignant detection: 45.23%
Diploid detection: 7.41%
Detection difference: 37.82 percentage points
Interpretation: Weak tumor association / requires further validation
COX8A
Malignant detection: 98.13%
Diploid detection: 80.69%
Detection difference: 17.44 percentage points
Interpretation: Weak tumor association

COX8A showed very high expression in malignant cells but was also highly detected in diploid cells, indicating lower tumor specificity.

External Validation

An independent dataset:

GSE176078

was used as an external validation dataset.

External validation was treated as an independent evidence layer rather than as proof of clinical biomarker validity.

The final evidence framework distinguishes between:

candidates with external support
candidates without sufficient external evidence
candidates requiring further validation

External validation does not establish clinical biomarker validity.

Final Evidence Integration

The final candidate assessment integrated multiple evidence layers:

Bioinformatics support
        +
ML recurrence
        +
Discovery expression
        +
External validation
        +
Tumor specificity
        ↓
Integrated candidate evidence

The integrated score generated by the project is an:

Exploratory computational evidence score

It should not be interpreted as:

probability of clinical validity
diagnostic probability
prognostic probability
treatment-response probability
statistical confidence

Final Candidate Interpretation
PRSS23

PRSS23 showed strong convergence across the computational workflow.

Evidence included:

3/3 ML recurrence
strong malignant-cell enrichment
high malignant-cell detection
consistent expression across informative samples
tumor-associated expression pattern

Interpretation: Strong tumor-associated research candidate

FCGRT

FCGRT showed:

3/3 ML recurrence
high malignant-cell detection
strong malignant-vs-diploid detection difference
consistent expression
tumor-associated expression

Interpretation: Strong tumor-associated research candidate

ETFB

ETFB showed:

3/3 ML recurrence
strong malignant enrichment
largest malignant-vs-diploid detection difference among the six
consistent expression
mitochondrial/metabolic biological context

Interpretation: Strong mitochondrial/metabolic research candidate

MRPS7

MRPS7 showed:

3/3 ML recurrence
high malignant-cell detection
strong expression difference
strong cross-sample expression validation
mitochondrial/ribosomal biological context

Interpretation: Strong mitochondrial/ribosomal research candidate

EPN3

EPN3 showed:

recurrent ML support
strong malignant-vs-diploid detection difference
low diploid detection

However, its sample consistency was weaker.

Interpretation: Promising candidate requiring additional validation

COX8A

COX8A showed:

3/3 ML recurrence
very high malignant-cell detection
consistent expression
strong computational evidence

However, it was also highly detected in diploid cells.

Interpretation: Strong expression candidate with lower tumor specificity

Final Candidate Summary

| Gene       | ML recurrence | Tumor association | Expression evidence | Final interpretation                            |
| ---------- | ------------- | ----------------- | ------------------- | ----------------------------------------------- |
| **PRSS23** | 3/3           | Strong            | Strong              | High-priority research candidate                |
| **FCGRT**  | 3/3           | Strong            | Strong              | High-priority research candidate                |
| **ETFB**   | 3/3           | Strong            | Strong              | High-priority metabolic candidate               |
| **MRPS7**  | 3/3           | Strong            | Strong              | High-priority mitochondrial/ribosomal candidate |
| **EPN3**   | 3/3           | Weak/moderate     | Less consistent     | Requires further validation                     |
| **COX8A**  | 3/3           | Weak              | Strong but broad    | Lower tumor specificity                         |

Key Findings
1. Progressive candidate reduction

The workflow reduced the candidate space from:

36,626 genes
        ↓
13,732 expression-filtered genes
        ↓
986 bioinformatics candidates
        ↓
50 final candidates
        ↓
6 consensus candidates

2. Sample-aware validation revealed generalization challenges

The initial 50-gene model achieved:

Mean ROC-AUC = 0.824
Mean PR-AUC  = 0.902

However, data-driven feature selection produced substantially weaker outer-sample performance.

This demonstrates that:

Strong internal validation does not necessarily translate into strong generalization to an unseen biological sample.

3. ML recurrence provided an additional evidence layer

Rather than selecting genes from a single model fit, recurrent feature selection across independent sample folds was used.

This identified six strong consensus genes.

4. Tumor-associated evidence differed among candidates

The strongest tumor-associated candidates were:

PRSS23, FCGRT, ETFB and MRPS7

EPN3 showed promising but less consistent tumor association.

COX8A showed broad expression and therefore lower tumor specificity.

Limitations
Limited number of informative samples

Only three samples contained malignant cells:

BC11
BC12
BC17

Therefore, the effective sample-level validation size is small.

Class imbalance

Some test samples contained very few diploid/non-malignant cells.

For example:

BC17

Malignant: 359
Diploid:     2

Therefore, sample-specific classification metrics should be interpreted cautiously.

Internal vs outer-sample performance

Some inner-validation results were extremely high, while outer LOPO performance was considerably lower.

This suggests potential sample-specific signal and reinforces the importance of sample-aware validation.

Computational malignant-cell identification

CNV-supported malignant-cell identification is computational evidence and is not equivalent to experimental pathological confirmation.

Association does not establish causality

A gene enriched in malignant cells does not automatically mean that the gene causes cancer development or progression.

Clinical validation has not been established

The final genes should be considered:

Computationally prioritized breast cancer biomarker research candidates

They are not clinically validated diagnostic, prognostic, or therapeutic biomarkers.

External validation limitations

GSE176078 provides an independent evidence layer, but larger independent cohorts are required to establish robust generalization.

Experimental validation

Experimental validation is required before any clinical interpretation.

Software and Tools
Programming
Python
R
Bash/Linux

Single-cell Analysis
Seurat
Seurat v5
Assay5

Machine Learning
scikit-learn
Logistic Regression
Feature Selection
Leave-One-Sample-Out Validation

Data Analysis
pandas
NumPy
SciPy
Visualization
Matplotlib
Seaborn

Bioinformatics
Single-cell RNA-seq analysis
CNV-supported malignant-cell analysis
Differential expression
Candidate prioritization
Biological annotation
External validation

Future Work

Future extensions could include:

Validation in larger independent single-cell datasets.
Validation using bulk breast cancer cohorts.
Survival/prognostic analysis.
Protein-level validation.
Spatial transcriptomics validation.
Functional pathway analysis.
qPCR validation.
Immunohistochemistry validation.
Functional experiments.
Evaluation of multi-gene biomarker panels.

Final Conclusion

This study developed a sample-aware single-cell RNA-seq and machine-learning workflow for breast cancer biomarker candidate prioritization.

By integrating bioinformatics candidate prioritization, CNV-supported malignant-cell identification, expression filtering, logistic regression, LOPO validation, ML feature recurrence, discovery expression validation, external validation using GSE176078, and tumor-specificity analysis, the workflow identified six consensus candidate genes:

COX8A, EPN3, ETFB, FCGRT, MRPS7 and PRSS23.

Among these candidates, PRSS23, FCGRT, ETFB and MRPS7 showed the strongest tumor-associated expression evidence in the discovery dataset.

EPN3 showed a strong malignant-vs-diploid detection difference but weaker sample consistency.

COX8A showed very high malignant-cell detection but relatively low tumor specificity because it was also highly detected in diploid cells.

An important methodological finding was that strong internal validation performance did not consistently translate into strong performance on completely held-out samples. This highlights the importance of sample-aware validation when applying machine learning to single-cell datasets.

Overall, the project demonstrates a strategy for combining bioinformatics, machine learning, sample-aware validation, expression analysis, external validation, and biological interpretation to prioritize candidate genes.

The final genes should be considered:

Computationally prioritized breast cancer biomarker research candidates, not clinically validated biomarkers.

Larger independent cohorts and experimental validation will be required to establish their diagnostic, prognostic, therapeutic, or causal relevance

Author

Shivajeet Yadav

Bioinformatics | Computational Biology | NGS | Single-Cell Genomics | Machine Learning.
