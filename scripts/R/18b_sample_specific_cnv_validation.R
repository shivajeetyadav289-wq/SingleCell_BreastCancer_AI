# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 18b_sample_specific_cnv_validation.R
#
# Purpose:
# Compare CopyKAT CNV architecture separately across samples
# with substantial aneuploid populations.
#
# Focus:
#   BC11
#   BC12
#   BC17
#
# IMPORTANT:
# - Does NOT rerun CopyKAT
# - Does NOT change CopyKAT predictions
# - Does NOT create final malignant labels
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})


# ------------------------------------------------------------
# 2. Paths
# ------------------------------------------------------------

project_dir <- path.expand(
  "~/Projects/SingleCell_BreastCancer_AI"
)

copykat_dir <- file.path(
  project_dir,
  "results",
  "cnv",
  "copykat"
)

integration_dir <- file.path(
  project_dir,
  "results",
  "cnv",
  "integration"
)

figure_dir <- file.path(
  project_dir,
  "figures",
  "cnv",
  "sample_specific"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 3. Input files
# ------------------------------------------------------------

cnv_file <- file.path(
  copykat_dir,
  "GSE228499_CopyKAT_CNAmat.rds"
)

integration_file <- file.path(
  integration_dir,
  "GSE228499_CNV_marker_integrated_results.csv"
)


# ------------------------------------------------------------
# 4. Header
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("SAMPLE-SPECIFIC CNV VALIDATION\n")
cat("====================================================\n")

cat(
  "Focus samples: BC11, BC12, BC17\n"
)

cat(
  "Using existing CopyKAT CNAmat\n"
)

cat(
  "CopyKAT will NOT be rerun.\n"
)


# ------------------------------------------------------------
# 5. Check files
# ------------------------------------------------------------

if (!file.exists(cnv_file)) {
  stop(
    "CopyKAT CNAmat file not found:\n",
    cnv_file
  )
}

if (!file.exists(integration_file)) {
  stop(
    "Integration file not found:\n",
    integration_file
  )
}


# ------------------------------------------------------------
# 6. Load data
# ------------------------------------------------------------

cat("\n")
cat("Loading CNV matrix...\n")

cnv <- readRDS(
  cnv_file
)

integration <- read.csv(
  integration_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "CNV dimensions:",
  nrow(cnv),
  "regions x",
  ncol(cnv) - 3,
  "cells\n"
)

cat(
  "Integration rows:",
  nrow(integration),
  "\n"
)


# ------------------------------------------------------------
# 7. Extract CNV matrix
# ------------------------------------------------------------

cnv_values <- as.matrix(
  cnv[
    ,
    4:ncol(cnv)
  ]
)

storage.mode(
  cnv_values
) <- "numeric"


# ------------------------------------------------------------
# 8. Chromosome information
# ------------------------------------------------------------

chrom <- toupper(
  trimws(
    as.character(
      cnv$chrom
    )
  )
)

chrom <- gsub(
  "^CHR",
  "",
  chrom
)

chrompos <- as.numeric(
  cnv$chrompos
)

canonical_chr <- c(
  as.character(1:22),
  "X",
  "Y"
)

keep_chr <- chrom %in%
  canonical_chr

cnv_values <- cnv_values[
  keep_chr,
  ,
  drop = FALSE
]

chrom <- chrom[
  keep_chr
]

chrompos <- chrompos[
  keep_chr
]


# ------------------------------------------------------------
# 9. Genomic ordering
# ------------------------------------------------------------

chr_order <- c(
  as.character(1:22),
  "X",
  "Y"
)

chr_factor <- factor(
  chrom,
  levels = chr_order
)

idx <- order(
  chr_factor,
  chrompos
)

cnv_values <- cnv_values[
  idx,
  ,
  drop = FALSE
]

chrom <- chrom[
  idx
]

chrompos <- chrompos[
  idx
]


# ------------------------------------------------------------
# 10. Cell ID normalization
# ------------------------------------------------------------

normalize_cell_id <- function(x) {

  x <- as.character(x)

  x <- gsub(
    "\\.([0-9]+)$",
    "-\\1",
    x
  )

  x

}


cnv_ids <- colnames(
  cnv_values
)

cnv_ids_normalized <- normalize_cell_id(
  cnv_ids
)

integration_ids <- normalize_cell_id(
  integration$cell_id
)


# ------------------------------------------------------------
# 11. Match metadata
# ------------------------------------------------------------

match_idx <- match(
  cnv_ids_normalized,
  integration_ids
)

matched <- !is.na(
  match_idx
)

cat(
  "Matched CNV cells:",
  sum(matched),
  "\n"
)

if (
  sum(matched) < 100
) {

  stop(
    "Insufficient metadata matching."
  )

}


# ------------------------------------------------------------
# 12. Metadata for CNV cells
# ------------------------------------------------------------

meta <- integration[
  match_idx,
  ,
  drop = FALSE
]

rownames(meta) <- cnv_ids


# ------------------------------------------------------------
# 13. Confirm CopyKAT predictions
# ------------------------------------------------------------

cat("\n")
cat("CopyKAT prediction summary:\n")

print(
  table(
    meta$copykat.pred,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# 14. Focus samples
# ------------------------------------------------------------

focus_samples <- c(
  "BC11",
  "BC12",
  "BC17"
)

focus_meta <- meta %>%

  filter(
    sample_id %in%
      focus_samples
  )


cat("\n")
cat("Focus sample cell counts:\n")

print(
  table(
    focus_meta$sample_id,
    focus_meta$copykat.pred
  )
)


# ------------------------------------------------------------
# 15. Calculate sample-specific profiles
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("CALCULATING SAMPLE-SPECIFIC CNV PROFILES\n")
cat("====================================================\n")


sample_profiles <- list()

sample_cell_counts <- data.frame()


for (
  s in focus_samples
) {

  ids <- rownames(
    meta
  )[
    meta$sample_id ==
      s &
    meta$copykat.pred ==
      "aneuploid"
  ]

  ids <- intersect(
    ids,
    colnames(cnv_values)
  )

  cat(
    s,
    ":",
    length(ids),
    "aneuploid cells\n"
  )

  if (
    length(ids) == 0
  ) {

    next

  }

  profile <- rowMeans(
    cnv_values[
      ,
      ids,
      drop = FALSE
    ],
    na.rm = TRUE
  )

  sample_profiles[[s]] <-
    profile

  sample_cell_counts <-
    rbind(
      sample_cell_counts,
      data.frame(
        sample_id = s,
        aneuploid_cells =
          length(ids)
      )
    )

}


# ------------------------------------------------------------
# 16. Combine profiles
# ------------------------------------------------------------

profile_matrix <- do.call(
  cbind,
  sample_profiles
)

colnames(
  profile_matrix
) <- names(
  sample_profiles
)


cat(
  "\nProfile matrix:",
  nrow(profile_matrix),
  "regions x",
  ncol(profile_matrix),
  "samples\n"
)


# ------------------------------------------------------------
# 17. Save profile matrix
# ------------------------------------------------------------

saveRDS(
  profile_matrix,
  file.path(
    integration_dir,
    "18b_sample_specific_CNV_profiles.rds"
  )
)


# ------------------------------------------------------------
# 18. Calculate chromosome-level sample profiles
# ------------------------------------------------------------

chromosome_profile <- data.frame(
  chromosome = chrom
)

for (
  s in colnames(
    profile_matrix
  )
) {

  chromosome_profile[[s]] <-
    profile_matrix[
      ,
      s
    ]

}


chromosome_summary <- chromosome_profile %>%

  group_by(
    chromosome
  ) %>%

  summarise(

    BC11 =
      mean(
        BC11,
        na.rm = TRUE
      ),

    BC12 =
      mean(
        BC12,
        na.rm = TRUE
      ),

    BC17 =
      mean(
        BC17,
        na.rm = TRUE
      ),

    .groups = "drop"

  )


# ------------------------------------------------------------
# 19. Save chromosome summary
# ------------------------------------------------------------

write.csv(
  chromosome_summary,
  file.path(
    integration_dir,
    "18b_sample_chromosome_CNV_summary.csv"
  ),
  row.names = FALSE
)


cat("\n")
cat("Chromosome-level CNV summary:\n")

print(
  chromosome_summary
)


# ------------------------------------------------------------
# 20. Calculate CNV burden
# ------------------------------------------------------------

burden <- data.frame()

for (
  s in colnames(
    profile_matrix
  )
) {

  values <- profile_matrix[
    ,
    s
  ]

  burden <- rbind(

    burden,

    data.frame(

      sample_id = s,

      mean_abs_CNV =
        mean(
          abs(values),
          na.rm = TRUE
        ),

      median_abs_CNV =
        median(
          abs(values),
          na.rm = TRUE
        ),

      gain_fraction =
        mean(
          values > 0.1,
          na.rm = TRUE
        ),

      loss_fraction =
        mean(
          values < -0.1,
          na.rm = TRUE
        ),

      net_CNV =
        mean(
          values,
          na.rm = TRUE
        )

    )

  )

}


cat("\n")
cat("Sample-specific CNV burden:\n")

print(
  burden
)


write.csv(
  burden,
  file.path(
    integration_dir,
    "18b_sample_specific_CNV_burden.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 21. Sample profile plot
# ------------------------------------------------------------

plot_data <- data.frame(
  chromosome = chrom
)

for (
  s in colnames(
    profile_matrix
  )
) {

  plot_data[[s]] <-
    profile_matrix[
      ,
      s
    ]

}


plot_long <- rbind(

  data.frame(
    chromosome =
      plot_data$chromosome,
    sample_id =
      "BC11",
    CNV =
      plot_data$BC11
  ),

  data.frame(
    chromosome =
      plot_data$chromosome,
    sample_id =
      "BC12",
    CNV =
      plot_data$BC12
  ),

  data.frame(
    chromosome =
      plot_data$chromosome,
    sample_id =
      "BC17",
    CNV =
      plot_data$BC17
  )

)


# Use genomic-region index to preserve chromosome order

plot_long$region_index <-
  rep(
    seq_len(
      nrow(plot_data)
    ),
    3
  )


# ------------------------------------------------------------
# 22. Plot sample CNV profiles
# ------------------------------------------------------------

pdf(
  file.path(
    figure_dir,
    "18b_sample_specific_CNV_profiles.pdf"
  ),
  width = 14,
  height = 8
)

print(

  ggplot(
    plot_long,
    aes(
      x = region_index,
      y = CNV
    )
  ) +

    geom_line(
      linewidth = 0.25
    ) +

    geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +

    facet_wrap(
      ~ sample_id,
      ncol = 1
    ) +

    labs(
      title =
        "Sample-Specific CNV Profiles in Aneuploid Cells",
      x =
        "Genomic position",
      y =
        "Mean CNV"
    ) +

    theme_classic()

)

dev.off()


# ------------------------------------------------------------
# 23. Heatmap of sample profiles
# ------------------------------------------------------------

cat("\n")
cat("Creating sample-specific CNV heatmap...\n")

# profile_matrix:
#   rows    = genomic regions
#   columns = samples
#
# image() requires:
#   nrow(z) == length(x)
#   ncol(z) == length(y)
#
# Therefore we keep the matrix in genomic-region x sample
# orientation and use x = genomic regions, y = samples.

heatmap_matrix <- profile_matrix

# Limit extreme values for visualization only.

heatmap_matrix[
  heatmap_matrix > 2
] <- 2

heatmap_matrix[
  heatmap_matrix < -2
] <- -2

# Remove rows containing non-finite values.

valid_heatmap_rows <- apply(
  heatmap_matrix,
  1,
  function(x) {
    all(is.finite(x))
  }
)

heatmap_matrix <- heatmap_matrix[
  valid_heatmap_rows,
  ,
  drop = FALSE
]

cat(
  "Heatmap matrix:",
  nrow(heatmap_matrix),
  "genomic regions x",
  ncol(heatmap_matrix),
  "samples\n"
)

# ------------------------------------------------------------
# Create heatmap
# ------------------------------------------------------------

pdf(
  file.path(
    figure_dir,
    "18b_sample_CNV_heatmap.pdf"
  ),
  width = 14,
  height = 5
)

par(
  mar = c(
    8,
    5,
    4,
    2
  )
)

image(
  x = seq_len(
    nrow(heatmap_matrix)
  ),

  y = seq_len(
    ncol(heatmap_matrix)
  ),

  z = heatmap_matrix,

  axes = FALSE,

  xlab =
    "Genomic regions",

  ylab =
    "Sample",

  main =
    "Mean CNV Profile: BC11 vs BC12 vs BC17"
)

# ------------------------------------------------------------
# Sample axis
# ------------------------------------------------------------

axis(
  2,

  at = seq_len(
    ncol(heatmap_matrix)
  ),

  labels =
    colnames(
      heatmap_matrix
    ),

  las = 2
)

# ------------------------------------------------------------
# Genomic-region axis
# ------------------------------------------------------------

axis(
  1,

  at = seq(
    1,
    nrow(heatmap_matrix),
    length.out = 10
  ),

  labels = FALSE
)

box()

dev.off()

cat(
  "Saved:",
  file.path(
    figure_dir,
    "18b_sample_CNV_heatmap.pdf"
  ),
  "\n"
)


# ------------------------------------------------------------
# 24. Compare samples
# ------------------------------------------------------------

sample_cor <- cor(
  profile_matrix,
  use = "pairwise.complete.obs",
  method = "pearson"
)


cat("\n")
cat("Correlation of sample CNV profiles:\n")

print(
  round(
    sample_cor,
    3
  )
)


write.csv(
  sample_cor,
  file.path(
    integration_dir,
    "18b_sample_CNV_profile_correlation.csv"
  )
)


# ------------------------------------------------------------
# 25. Save interpretation table
# ------------------------------------------------------------

interpretation <- data.frame(

  sample_id =
    burden$sample_id,

  aneuploid_cells =
    sample_cell_counts$aneuploid_cells[
      match(
        burden$sample_id,
        sample_cell_counts$sample_id
      )
    ],

  mean_abs_CNV =
    burden$mean_abs_CNV,

  median_abs_CNV =
    burden$median_abs_CNV,

  gain_fraction =
    burden$gain_fraction,

  loss_fraction =
    burden$loss_fraction,

  net_CNV =
    burden$net_CNV

)


write.csv(
  interpretation,
  file.path(
    integration_dir,
    "18b_sample_specific_CNV_interpretation.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 26. Final report
# ------------------------------------------------------------

report_file <- file.path(
  integration_dir,
  "18b_sample_specific_CNV_validation_summary.txt"
)


report <- c(

  "SAMPLE-SPECIFIC CNV VALIDATION",

  "==============================",

  "",

  "Samples evaluated:",

  "BC11, BC12, BC17",

  "",

  paste(
    "BC11 aneuploid cells:",
    sample_cell_counts$aneuploid_cells[
      sample_cell_counts$sample_id == "BC11"
    ]
  ),

  paste(
    "BC12 aneuploid cells:",
    sample_cell_counts$aneuploid_cells[
      sample_cell_counts$sample_id == "BC12"
    ]
  ),

  paste(
    "BC17 aneuploid cells:",
    sample_cell_counts$aneuploid_cells[
      sample_cell_counts$sample_id == "BC17"
    ]
  ),

  "",

  "CNV profile correlation:",

  paste(
    capture.output(
      print(
        round(
          sample_cor,
          3
        )
      )
    ),
    collapse = "\n"
  ),

  "",

  "Interpretation:",

  "Sample-specific CNV profiles were calculated",

  "from CopyKAT-predicted aneuploid cells.",

  "Similarity between samples should be evaluated",

  "using the chromosome-level profiles and heatmap.",

  "Consistent CNV architecture strengthens",

  "CNV-supported malignant-cell interpretation.",

  "Different architectures may indicate",

  "distinct malignant subpopulations.",

  "These are research annotations, not clinical calls."

)

writeLines(
  report,
  report_file
)


# ------------------------------------------------------------
# 27. Final summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("18B SAMPLE-SPECIFIC CNV VALIDATION COMPLETED\n")
cat("====================================================\n")

cat("\nResults:\n")

cat(
  file.path(
    integration_dir,
    "18b_sample_specific_CNV_burden.csv"
  ),
  "\n"
)

cat(
  file.path(
    integration_dir,
    "18b_sample_chromosome_CNV_summary.csv"
  ),
  "\n"
)

cat(
  file.path(
    integration_dir,
    "18b_sample_CNV_profile_correlation.csv"
  ),
  "\n"
)

cat(
  file.path(
    integration_dir,
    "18b_sample_specific_CNV_interpretation.csv"
  ),
  "\n"
)

cat(
  report_file,
  "\n"
)

cat("\nFigures:\n")

cat(
  file.path(
    figure_dir,
    "18b_sample_specific_CNV_profiles.pdf"
  ),
  "\n"
)

cat(
  file.path(
    figure_dir,
    "18b_sample_CNV_heatmap.pdf"
  ),
  "\n"
)

cat("\n")
cat("====================================================\n")
cat("READY FOR FINAL MALIGNANT-CELL INTERPRETATION\n")
cat("====================================================\n")