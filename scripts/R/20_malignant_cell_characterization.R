# ============================================================
# STEP 20: MALIGNANT CELL CHARACTERIZATION
# ============================================================
#
# Project:
# AI-Assisted Discovery of Breast Cancer Biomarkers
# Using scRNA-seq and Machine Learning
#
# Dataset:
# GSE228499
#
# Main comparison:
# CNV-supported malignant epithelial cells
#             VS
# CopyKAT-diploid epithelial cells
#
# Additional comparison:
# Luminal malignant epithelial
#             VS
# Basal-like malignant epithelial
#
# IMPORTANT:
# These are research annotations and are NOT
# clinical diagnostic classifications.
# ============================================================


# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

suppressPackageStartupMessages({

  library(Seurat)
  library(dplyr)
  library(ggplot2)

})


# ------------------------------------------------------------
# 2. Project directories
# ------------------------------------------------------------

project_dir <- path.expand(
  "~/Projects/SingleCell_BreastCancer_AI"
)

input_file <- file.path(
  project_dir,
  "data",
  "processed",
  "GSE228499_final_malignant_annotated.rds"
)

result_dir <- file.path(
  project_dir,
  "results",
  "malignant",
  "characterization"
)

figure_dir <- file.path(
  project_dir,
  "figures",
  "malignant",
  "characterization"
)

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. Start message
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("STEP 20: MALIGNANT CELL CHARACTERIZATION\n")
cat("====================================================\n")


# ------------------------------------------------------------
# 4. Check input file
# ------------------------------------------------------------

if (!file.exists(input_file)) {

  stop(
    "\nInput Seurat object not found:\n",
    input_file,
    "\n"
  )

}


# ------------------------------------------------------------
# 5. Load Seurat object
# ------------------------------------------------------------

cat("\n")
cat("Loading Seurat object...\n")

obj <- readRDS(
  input_file
)

cat(
  "Genes:",
  nrow(obj),
  "\n"
)

cat(
  "Cells:",
  ncol(obj),
  "\n"
)

cat(
  "Default assay:",
  DefaultAssay(obj),
  "\n"
)


# ------------------------------------------------------------
# 6. Validate required metadata
# ------------------------------------------------------------

required_metadata <- c(

  "sample_id",
  "cell_type",
  "malignant_candidate",
  "malignant_status",
  "copykat_prediction",
  "epithelial_evidence"

)

missing_metadata <- setdiff(

  required_metadata,

  colnames(
    obj@meta.data
  )

)

if (length(missing_metadata) > 0) {

  stop(

    "\nMissing required metadata columns:\n",

    paste(
      missing_metadata,
      collapse = "\n"
    ),

    "\n"

  )

}


# ------------------------------------------------------------
# 7. Identify malignant candidates
# ------------------------------------------------------------

malignant_cells <- colnames(obj)[
  obj$malignant_candidate
]

cat("\n")

cat(
  "CNV-supported malignant candidates:",
  length(malignant_cells),
  "\n"
)

if (length(malignant_cells) == 0) {

  stop(
    "No malignant candidates were identified."
  )

}


# ------------------------------------------------------------
# 8. Create malignant object
# ------------------------------------------------------------

malignant_obj <- subset(
  obj,
  cells = malignant_cells
)

cat(
  "Malignant object:",
  nrow(malignant_obj),
  "genes x",
  ncol(malignant_obj),
  "cells\n"
)


# ============================================================
# PART A
# MALIGNANT CELL COMPOSITION
# ============================================================


# ------------------------------------------------------------
# 9. Cell-type composition
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("MALIGNANT CELL-TYPE COMPOSITION\n")
cat("====================================================\n")

celltype_summary <- malignant_obj@meta.data %>%

  count(
    cell_type,
    name = "cells"
  ) %>%

  mutate(
    percentage =
      100 *
      cells /
      sum(cells)
  ) %>%

  arrange(
    desc(cells)
  )


print(
  celltype_summary
)


write.csv(

  celltype_summary,

  file.path(
    result_dir,
    "20_malignant_celltype_composition.csv"
  ),

  row.names = FALSE

)


# ------------------------------------------------------------
# 10. Malignant cells by sample
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("MALIGNANT CELLS BY SAMPLE\n")
cat("====================================================\n")

sample_summary <- malignant_obj@meta.data %>%

  count(
    sample_id,
    name = "malignant_cells"
  ) %>%

  mutate(
    percentage =
      100 *
      malignant_cells /
      sum(malignant_cells)
  ) %>%

  arrange(
    desc(malignant_cells)
  )


print(
  sample_summary
)


write.csv(

  sample_summary,

  file.path(
    result_dir,
    "20_malignant_cells_by_sample.csv"
  ),

  row.names = FALSE

)


# ============================================================
# PART B
# UMAP VISUALIZATION
# ============================================================


# ------------------------------------------------------------
# 11. UMAP visualizations
# ------------------------------------------------------------

if ("umap" %in% Reductions(obj)) {

  cat("\n")
  cat("Creating UMAP visualizations...\n")


  # ----------------------------------------------------------
  # Overall malignant status
  # ----------------------------------------------------------

  pdf(

    file.path(
      figure_dir,
      "20_malignant_status_UMAP.pdf"
    ),

    width = 10,
    height = 8

  )

  print(

    DimPlot(

      obj,

      reduction = "umap",

      group.by = "malignant_status",

      raster = TRUE

    ) +

      ggtitle(
        "CNV-Supported Malignant Cell Candidates"
      )

  )

  dev.off()


  # ----------------------------------------------------------
  # Malignant cell type
  # ----------------------------------------------------------

  if ("umap" %in% Reductions(malignant_obj)) {

    pdf(

      file.path(
        figure_dir,
        "20_malignant_cells_UMAP_celltype.pdf"
      ),

      width = 10,
      height = 8

    )

    print(

      DimPlot(

        malignant_obj,

        reduction = "umap",

        group.by = "cell_type",

        raster = TRUE

      ) +

        ggtitle(
          "Cell-Type Composition of Malignant Candidates"
        )

    )

    dev.off()


    # --------------------------------------------------------
    # Malignant cells by sample
    # --------------------------------------------------------

    pdf(

      file.path(
        figure_dir,
        "20_malignant_cells_UMAP_sample.pdf"
      ),

      width = 10,
      height = 8

    )

    print(

      DimPlot(

        malignant_obj,

        reduction = "umap",

        group.by = "sample_id",

        raster = TRUE

      ) +

        ggtitle(
          "Malignant Candidates by Sample"
        )

    )

    dev.off()

  }

}


# ============================================================
# PART C
# MALIGNANT VS DIPLOID EPITHELIAL
# ============================================================


# ------------------------------------------------------------
# 12. Identify diploid epithelial reference
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("CREATING DIPLOID EPITHELIAL REFERENCE\n")
cat("====================================================\n")


reference_cells <- colnames(obj)[

  obj$epithelial_evidence &

  obj$copykat_prediction ==
    "diploid"

]


cat(

  "Diploid epithelial reference cells:",

  length(reference_cells),

  "\n"

)


if (length(reference_cells) == 0) {

  stop(
    "No diploid epithelial reference cells found."
  )

}


# ------------------------------------------------------------
# 13. Create DE object
# ------------------------------------------------------------

de_cells <- unique(

  c(
    malignant_cells,
    reference_cells
  )

)


de_obj <- subset(

  obj,

  cells = de_cells

)


cat("\n")
cat("DE object:\n")

cat(
  "Genes:",
  nrow(de_obj),
  "\n"
)

cat(
  "Cells:",
  ncol(de_obj),
  "\n"
)


# ------------------------------------------------------------
# 14. Define DE groups
# ------------------------------------------------------------

de_obj$malignant_comparison <-
  "Diploid_Epithelial"


de_obj$malignant_comparison[
  de_obj$malignant_candidate
] <-
  "CNV_Malignant_Candidate"


cat("\n")
cat("DE groups:\n")

print(

  table(
    de_obj$malignant_comparison
  )

)


# ------------------------------------------------------------
# 15. Set RNA assay
# ------------------------------------------------------------

DefaultAssay(
  de_obj
) <- "RNA"


# ------------------------------------------------------------
# 16. Check layers before joining
# ------------------------------------------------------------

cat("\n")
cat("RNA layers BEFORE JoinLayers:\n")

print(

  Layers(
    de_obj[["RNA"]]
  )

)


# ------------------------------------------------------------
# 17. Join Seurat v5 RNA layers
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("JOINING RNA LAYERS\n")
cat("====================================================\n")


de_obj <- JoinLayers(

  de_obj,

  assay = "RNA"

)


cat("\n")
cat("RNA layers AFTER JoinLayers:\n")

print(

  Layers(
    de_obj[["RNA"]]
  )

)


# ------------------------------------------------------------
# 18. Validate joined data layer
# ------------------------------------------------------------

de_layers <- Layers(
  de_obj[["RNA"]]
)


if (!"data" %in% de_layers) {

  stop(

    "\nJoined RNA data layer was not found.\n",

    "Available layers:\n",

    paste(
      de_layers,
      collapse = ", "
    ),

    "\n"

  )

}


# ------------------------------------------------------------
# 19. Run malignant vs diploid DE
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("MALIGNANT VS DIPLOID DIFFERENTIAL EXPRESSION\n")
cat("====================================================\n")


cat(
  "Comparison:\n"
)

cat(
  "CNV-supported malignant candidates",
  "vs",
  "diploid epithelial cells\n"
)


cat(

  "Malignant cells:",

  sum(
    de_obj$malignant_comparison ==
      "CNV_Malignant_Candidate"
  ),

  "\n"

)


cat(

  "Diploid epithelial cells:",

  sum(
    de_obj$malignant_comparison ==
      "Diploid_Epithelial"
  ),

  "\n"

)


cat(
  "Running Wilcoxon differential expression...\n"
)


de_results <- FindMarkers(

  object = de_obj,

  ident.1 =
    "CNV_Malignant_Candidate",

  ident.2 =
    "Diploid_Epithelial",

  group.by =
    "malignant_comparison",

  assay =
    "RNA",

  slot =
    "data",

  test.use =
    "wilcox",

  logfc.threshold =
    0.25,

  min.pct =
    0.10,

  only.pos =
    FALSE,

  verbose =
    TRUE

)


cat("\n")
cat(
  "Differential expression completed.\n"
)

cat(
  "Genes tested:",
  nrow(de_results),
  "\n"
)


# ------------------------------------------------------------
# 20. Standardize DE column names
# ------------------------------------------------------------

cat("\n")
cat("DE result columns:\n")

print(
  colnames(de_results)
)


if (
  "avg_log2FC" %in%
  colnames(de_results)
) {

  # Seurat already provides avg_log2FC

} else if (
  "avg_logFC" %in%
  colnames(de_results)
) {

  de_results$avg_log2FC <-
    de_results$avg_logFC

} else {

  stop(
    "\nCould not identify log fold-change column.\n"
  )

}


# ------------------------------------------------------------
# 21. Ensure percentage columns exist
# ------------------------------------------------------------

if (
  !"pct.1" %in%
  colnames(de_results)
) {

  de_results$pct.1 <-
    NA_real_

}


if (
  !"pct.2" %in%
  colnames(de_results)
) {

  de_results$pct.2 <-
    NA_real_

}


if (
  !"p_val_adj" %in%
  colnames(de_results)
) {

  stop(
    "p_val_adj column missing from DE results."
  )

}


# ------------------------------------------------------------
# 22. Add gene names
# ------------------------------------------------------------

de_results$gene <-
  rownames(
    de_results
  )


# ------------------------------------------------------------
# 23. Sort DE results
# ------------------------------------------------------------

de_results <-

  de_results %>%

  arrange(

    p_val_adj,

    desc(
      abs(
        avg_log2FC
      )
    )

  )


# ------------------------------------------------------------
# 24. Save complete DE results
# ------------------------------------------------------------

write.csv(

  de_results,

  file.path(
    result_dir,
    "20_malignant_vs_diploid_epithelial_DEG.csv"
  ),

  row.names = FALSE

)


# ------------------------------------------------------------
# 25. Significant malignant-associated genes
# ------------------------------------------------------------

significant_de <-

  de_results %>%

  filter(
    p_val_adj < 0.05
  ) %>%

  arrange(
    p_val_adj
  )


write.csv(

  significant_de,

  file.path(
    result_dir,
    "20_significant_malignant_associated_genes.csv"
  ),

  row.names = FALSE

)


cat("\n")

cat(
  "Significant genes:",
  nrow(significant_de),
  "\n"
)


# ------------------------------------------------------------
# 26. Top malignant-associated genes
# ------------------------------------------------------------

top_up <-

  de_results %>%

  filter(

    avg_log2FC > 0,

    p_val_adj < 0.05

  ) %>%

  arrange(

    p_val_adj

  ) %>%

  head(30)


cat("\n")
cat("Top malignant-associated genes:\n")

print(

  top_up[

    ,

    c(
      "gene",
      "avg_log2FC",
      "pct.1",
      "pct.2",
      "p_val_adj"
    )

  ]

)


write.csv(

  top_up,

  file.path(
    result_dir,
    "20_top_malignant_associated_genes.csv"
  ),

  row.names = FALSE

)


# ------------------------------------------------------------
# 27. Top downregulated genes
# ------------------------------------------------------------

top_down <-

  de_results %>%

  filter(

    avg_log2FC < 0,

    p_val_adj < 0.05

  ) %>%

  arrange(

    p_val_adj

  ) %>%

  head(30)


write.csv(

  top_down,

  file.path(
    result_dir,
    "20_top_downregulated_genes.csv"
  ),

  row.names = FALSE

)


# ============================================================
# PART D
# VOLCANO PLOT
# ============================================================


# ------------------------------------------------------------
# 28. Volcano data
# ------------------------------------------------------------

volcano_data <-

  de_results %>%

  mutate(

    significance = case_when(

      p_val_adj < 0.05 &
        avg_log2FC > 0.25
      ~ "Upregulated",

      p_val_adj < 0.05 &
        avg_log2FC < -0.25
      ~ "Downregulated",

      TRUE
      ~ "Not_significant"

    ),

    neg_log10_padj =

      -log10(

        pmax(
          p_val_adj,
          1e-300
        )

      )

  )


# ------------------------------------------------------------
# 29. Volcano plot
# ------------------------------------------------------------

pdf(

  file.path(
    figure_dir,
    "20_malignant_vs_diploid_volcano.pdf"
  ),

  width = 10,

  height = 8

)


print(

  ggplot(

    volcano_data,

    aes(

      x = avg_log2FC,

      y = neg_log10_padj

    )

  ) +

    geom_point(

      alpha = 0.5,

      size = 1

    ) +

    geom_vline(

      xintercept =
        c(
          -0.25,
          0.25
        ),

      linetype =
        "dashed"

    ) +

    geom_hline(

      yintercept =
        -log10(0.05),

      linetype =
        "dashed"

    ) +

    labs(

      title =
        "Malignant Candidates vs Diploid Epithelial Cells",

      x =
        "Average log2 fold change",

      y =
        "-log10 adjusted P value"

    ) +

    theme_classic()

)


dev.off()


# ============================================================
# PART E
# TOP MALIGNANT GENE DOTPLOT
# ============================================================


# ------------------------------------------------------------
# 30. Select top genes
# ------------------------------------------------------------

top_genes <-
  unique(
    top_up$gene
  )


top_genes <-

  top_genes[
    top_genes %in%
      rownames(
        de_obj
      )
  ]


top_genes <-
  head(
    top_genes,
    20
  )


# ------------------------------------------------------------
# 31. DotPlot
# ------------------------------------------------------------

if (
  length(top_genes) > 0
) {

  pdf(

    file.path(
      figure_dir,
      "20_top_malignant_genes_DotPlot.pdf"
    ),

    width = 12,

    height = 8

  )


  print(

    DotPlot(

      de_obj,

      features =
        top_genes,

      group.by =
        "malignant_comparison"

    ) +

      RotatedAxis() +

      ggtitle(
        "Top Malignant-Associated Genes"
      )

  )


  dev.off()

}


# ============================================================
# PART F
# LUMINAL VS BASAL-LIKE MALIGNANT CELLS
# ============================================================


# ------------------------------------------------------------
# 32. Identify malignant subtypes
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("LUMINAL VS BASAL-LIKE MALIGNANT CELLS\n")
cat("====================================================\n")


malignant_subtypes <- c(

  "Luminal_Epithelial",

  "Basal_Like_Epithelial"

)


subtype_cells <- colnames(

  malignant_obj

)[

  malignant_obj$cell_type %in%
    malignant_subtypes

]


cat(

  "Eligible malignant subtype cells:",

  length(
    subtype_cells
  ),

  "\n"

)


if (
  length(subtype_cells) > 20
) {


  # ----------------------------------------------------------
  # 33. Create subtype object
  # ----------------------------------------------------------

  subtype_obj <- subset(

    malignant_obj,

    cells =
      subtype_cells

  )


  subtype_obj$malignant_subtype <-

    as.character(

      subtype_obj$cell_type

    )


  cat("\n")
  cat("Subtype counts:\n")

  print(

    table(
      subtype_obj$malignant_subtype
    )

  )


  # ----------------------------------------------------------
  # 34. Set RNA assay
  # ----------------------------------------------------------

  DefaultAssay(
    subtype_obj
  ) <- "RNA"


  # ----------------------------------------------------------
  # 35. Check subtype layers
  # ----------------------------------------------------------

  cat("\n")
  cat("Subtype RNA layers BEFORE JoinLayers:\n")

  print(

    Layers(
      subtype_obj[["RNA"]]
    )

  )


  # ----------------------------------------------------------
  # 36. Join subtype RNA layers
  # ----------------------------------------------------------

  subtype_obj <-

    JoinLayers(

      subtype_obj,

      assay =
        "RNA"

    )


  cat("\n")
  cat("Subtype RNA layers AFTER JoinLayers:\n")

  print(

    Layers(
      subtype_obj[["RNA"]]
    )

  )


  subtype_layers <-

    Layers(
      subtype_obj[["RNA"]]
    )


  if (
    !"data" %in%
    subtype_layers
  ) {

    stop(
      "Joined subtype RNA data layer unavailable."
    )

  }


  # ----------------------------------------------------------
  # 37. Run subtype DE
  # ----------------------------------------------------------

  cat("\n")
  cat(
    "Running luminal vs basal-like DE...\n"
  )


  subtype_de <-

    FindMarkers(

      object =
        subtype_obj,

      ident.1 =
        "Luminal_Epithelial",

      ident.2 =
        "Basal_Like_Epithelial",

      group.by =
        "malignant_subtype",

      assay =
        "RNA",

      slot =
        "data",

      test.use =
        "wilcox",

      logfc.threshold =
        0.25,

      min.pct =
        0.10,

      only.pos =
        FALSE,

      verbose =
        TRUE

    )


  cat("\n")
  cat(
    "Luminal vs basal-like DE completed.\n"
  )

  cat(
    "Genes tested:",
    nrow(subtype_de),
    "\n"
  )


  # ----------------------------------------------------------
  # 38. Standardize subtype logFC
  # ----------------------------------------------------------

  cat("\n")
  cat("Subtype DE columns:\n")

  print(
    colnames(
      subtype_de
    )
  )


  if (
    "avg_log2FC" %in%
    colnames(subtype_de)
  ) {

    # Already available

  } else if (
    "avg_logFC" %in%
    colnames(subtype_de)
  ) {

    subtype_de$avg_log2FC <-
      subtype_de$avg_logFC

  } else {

    stop(
      "Could not identify subtype log fold-change column."
    )

  }


  # ----------------------------------------------------------
  # 39. Add subtype gene names
  # ----------------------------------------------------------

  subtype_de$gene <-
    rownames(
      subtype_de
    )


  # ----------------------------------------------------------
  # 40. Ensure subtype columns
  # ----------------------------------------------------------

  if (
    !"pct.1" %in%
    colnames(subtype_de)
  ) {

    subtype_de$pct.1 <-
      NA_real_

  }


  if (
    !"pct.2" %in%
    colnames(subtype_de)
  ) {

    subtype_de$pct.2 <-
      NA_real_

  }


  if (
    !"p_val_adj" %in%
    colnames(subtype_de)
  ) {

    stop(
      "Subtype DE p_val_adj column missing."
    )

  }


  # ----------------------------------------------------------
  # 41. Sort subtype results
  # ----------------------------------------------------------

  subtype_de <-

    subtype_de %>%

    arrange(

      p_val_adj,

      desc(
        abs(
          avg_log2FC
        )
      )

    )


  # ----------------------------------------------------------
  # 42. Save complete subtype DE
  # ----------------------------------------------------------

  write.csv(

    subtype_de,

    file.path(
      result_dir,
      "20_luminal_vs_basal_malignant_DEG.csv"
    ),

    row.names = FALSE

  )


  # ----------------------------------------------------------
  # 43. Significant subtype genes
  # ----------------------------------------------------------

  subtype_significant <-

    subtype_de %>%

    filter(
      p_val_adj < 0.05
    )


  write.csv(

    subtype_significant,

    file.path(
      result_dir,
      "20_luminal_vs_basal_significant_genes.csv"
    ),

    row.names = FALSE

  )


  cat("\n")

  cat(

    "Significant luminal-vs-basal genes:",

    nrow(
      subtype_significant
    ),

    "\n"

  )


  # ----------------------------------------------------------
  # 44. Top luminal-associated genes
  # ----------------------------------------------------------

  top_luminal <-

    subtype_de %>%

    filter(

      avg_log2FC > 0,

      p_val_adj < 0.05

    ) %>%

    arrange(
      p_val_adj
    ) %>%

    head(30)


  write.csv(

    top_luminal,

    file.path(
      result_dir,
      "20_top_luminal_malignant_genes.csv"
    ),

    row.names = FALSE

  )


  # ----------------------------------------------------------
  # 45. Top basal-associated genes
  # ----------------------------------------------------------

  top_basal <-

    subtype_de %>%

    filter(

      avg_log2FC < 0,

      p_val_adj < 0.05

    ) %>%

    arrange(
      p_val_adj
    ) %>%

    head(30)


  write.csv(

    top_basal,

    file.path(
      result_dir,
      "20_top_basal_malignant_genes.csv"
    ),

    row.names = FALSE

  )


  # ----------------------------------------------------------
  # 46. Print top subtype genes
  # ----------------------------------------------------------

  cat("\n")
  cat("Top luminal-associated genes:\n")

  print(

    top_luminal[

      ,

      c(
        "gene",
        "avg_log2FC",
        "pct.1",
        "pct.2",
        "p_val_adj"
      )

    ]

  )


  cat("\n")
  cat("Top basal-associated genes:\n")

  print(

    top_basal[

      ,

      c(
        "gene",
        "avg_log2FC",
        "pct.1",
        "pct.2",
        "p_val_adj"
      )

    ]

  )


} else {

  cat(
    "\nInsufficient malignant subtype cells ",
    "for luminal-vs-basal comparison.\n"
  )

}


# ============================================================
# PART G
# CANONICAL BREAST-CANCER MARKERS
# ============================================================


# ------------------------------------------------------------
# 47. Marker list
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("CANONICAL BREAST-CANCER MARKER ANALYSIS\n")
cat("====================================================\n")


candidate_markers <- c(

  "EPCAM",

  "KRT8",

  "KRT18",

  "KRT19",

  "KRT5",

  "KRT14",

  "KRT17",

  "MKI67",

  "TOP2A",

  "MMP9",

  "MMP11",

  "CEACAM6",

  "MUC1",

  "ERBB2",

  "ESR1",

  "PGR"

)


available_markers <-

  candidate_markers[

    candidate_markers %in%

      rownames(
        malignant_obj
      )

  ]


cat(

  "Available canonical markers:",

  length(
    available_markers
  ),

  "\n"

)


print(
  available_markers
)


# ------------------------------------------------------------
# 48. Canonical marker DotPlot
# ------------------------------------------------------------

if (
  length(available_markers) > 0
) {

  pdf(

    file.path(
      figure_dir,
      "20_malignant_canonical_marker_DotPlot.pdf"
    ),

    width = 12,

    height = 8

  )


  print(

    DotPlot(

      malignant_obj,

      features =
        available_markers,

      group.by =
        "cell_type"

    ) +

      RotatedAxis() +

      ggtitle(
        "Canonical Marker Expression in Malignant Candidates"
      )

  )


  dev.off()

}


# ------------------------------------------------------------
# 49. Marker expression summary
# ------------------------------------------------------------

marker_summary <- data.frame()


for (
  gene in available_markers
) {

  values <-

    FetchData(

      malignant_obj,

      vars =
        gene

    )[, 1]


  marker_summary <-

    rbind(

      marker_summary,

      data.frame(

        gene =
          gene,

        mean_expression =

          mean(
            values,
            na.rm = TRUE
          ),

        expressing_fraction =

          mean(
            values > 0,
            na.rm = TRUE
          )

      )

    )

}


write.csv(

  marker_summary,

  file.path(
    result_dir,
    "20_malignant_canonical_marker_summary.csv"
  ),

  row.names = FALSE

)


# ============================================================
# PART H
# SAVE OBJECT AND REPORT
# ============================================================


# ------------------------------------------------------------
# 50. Save malignant characterization object
# ------------------------------------------------------------

characterization_file <-

  file.path(

    result_dir,

    "GSE228499_malignant_characterization.rds"

  )


saveRDS(

  malignant_obj,

  characterization_file

)


# ------------------------------------------------------------
# 51. Final report
# ------------------------------------------------------------

report_file <-

  file.path(

    result_dir,

    "20_malignant_cell_characterization_summary.txt"

  )


report <- c(

  "MALIGNANT CELL CHARACTERIZATION SUMMARY",

  "=========================================",

  "",

  paste(
    "Total malignant candidates:",
    ncol(malignant_obj)
  ),

  "",

  "Cell-type composition:",

  paste(

    capture.output(

      print(
        celltype_summary
      )

    ),

    collapse = "\n"

  ),

  "",

  "Sample composition:",

  paste(

    capture.output(

      print(
        sample_summary
      )

    ),

    collapse = "\n"

  ),

  "",

  paste(

    "Genes tested in malignant vs diploid comparison:",

    nrow(de_results)

  ),

  "",

  paste(

    "Significant malignant-associated genes:",

    nrow(significant_de)

  ),

  "",

  "Primary comparison:",

  "CNV-supported malignant epithelial cells",

  "versus",

  "CopyKAT-diploid epithelial cells",

  "",

  "Secondary comparison:",

  "Luminal malignant epithelial cells",

  "versus",

  "Basal-like malignant epithelial cells",

  "",

  "Interpretation:",

  "Differentially expressed genes represent",

  "genes associated with the analyzed",

  "malignant populations.",

  "",

  "These results require biological",

  "validation and are not clinical",

  "diagnostic classifications."

)


writeLines(

  report,

  report_file

)


# ============================================================
# FINAL SUMMARY
# ============================================================

cat("\n")
cat("====================================================\n")
cat("STEP 20 COMPLETED\n")
cat("====================================================\n")


cat(

  "Malignant candidates:",

  ncol(
    malignant_obj
  ),

  "\n"

)


cat(

  "Genes tested:",

  nrow(
    de_results
  ),

  "\n"

)


cat(

  "Significant malignant-associated genes:",

  nrow(
    significant_de
  ),

  "\n"

)


cat("\n")
cat("Results directory:\n")

cat(
  result_dir,
  "\n"
)


cat("\n")
cat("Figures directory:\n")

cat(
  figure_dir,
  "\n"
)


cat("\n")
cat("====================================================\n")
cat("READY FOR BIOMARKER PRIORITIZATION\n")
cat("====================================================\n")