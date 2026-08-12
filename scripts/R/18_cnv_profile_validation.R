# ============================================================
# Project: AI-Assisted Discovery of Breast Cancer Biomarkers
# Dataset: GSE228499
# Script: 18_cnv_profile_validation.R
#
# Purpose:
# Validate CopyKAT CNV profiles for the epithelial
# tumor-candidate population.
#
# Focus:
#   - Aneuploid vs diploid CNV profiles
#   - BC11, BC12 and BC17
#   - Chromosome-level CNV patterns
#   - Sample-specific CNV signal
#
# IMPORTANT:
# This script does NOT rerun CopyKAT.
# This script does NOT change malignant classifications.
# ============================================================


# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})


# ------------------------------------------------------------
# 2. Project paths
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
  "profile_validation"
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

prediction_file <- file.path(
  copykat_dir,
  "GSE228499_CopyKAT_prediction.rds"
)


# ------------------------------------------------------------
# 4. Header
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("CNV PROFILE VALIDATION\n")
cat("====================================================\n")

cat(
  "Dataset: GSE228499\n"
)

cat(
  "Using existing CopyKAT CNAmat\n"
)

cat(
  "No CopyKAT rerun will be performed.\n"
)


# ------------------------------------------------------------
# 5. Check input files
# ------------------------------------------------------------

required_files <- c(
  cnv_file,
  integration_file,
  prediction_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (
  length(missing_files) > 0
) {

  stop(
    "Missing required files:\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )

}


# ------------------------------------------------------------
# 6. Load CNV matrix
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("LOADING COPYKAT CNV MATRIX\n")
cat("====================================================\n")

cnv <- readRDS(
  cnv_file
)

cat(
  "CNV dimensions:",
  nrow(cnv),
  "rows x",
  ncol(cnv),
  "columns\n"
)

cat(
  "CNV class:",
  class(cnv)[1],
  "\n"
)


# ------------------------------------------------------------
# 7. Validate CNV matrix
# ------------------------------------------------------------

if (
  !all(
    c(
      "chrom",
      "chrompos",
      "abspos"
    ) %in%
      colnames(cnv)[1:3]
  )
) {

  stop(
    "Expected CopyKAT chromosome annotation columns were not found."
  )

}


# ------------------------------------------------------------
# 8. Load integration table
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("LOADING INTEGRATED CELL INFORMATION\n")
cat("====================================================\n")

integration <- read.csv(
  integration_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat(
  "Integration rows:",
  nrow(integration),
  "\n"
)


# ------------------------------------------------------------
# 9. Load CopyKAT prediction
# ------------------------------------------------------------

prediction <- readRDS(
  prediction_file
)

cat(
  "CopyKAT prediction rows:",
  nrow(prediction),
  "\n"
)


# ------------------------------------------------------------
# 10. Extract CNV cell names
# ------------------------------------------------------------

cnv_cell_names <- colnames(
  cnv
)[
  4:ncol(cnv)
]

cat(
  "CNV cell profiles:",
  length(cnv_cell_names),
  "\n"
)


# ------------------------------------------------------------
# 11. Normalize cell IDs
# ------------------------------------------------------------

normalize_cell_id <- function(x) {

  x <- as.character(x)

  # R can convert "-" to "." when check.names = TRUE.
  # Restore the original Seurat-style suffix.

  x <- gsub(
    "\\.([0-9]+)$",
    "-\\1",
    x
  )

  x

}


cnv_ids_normalized <- normalize_cell_id(
  cnv_cell_names
)

integration_ids <- normalize_cell_id(
  integration$cell_id
)


# ------------------------------------------------------------
# 12. Match CNV cells to integration metadata
# ------------------------------------------------------------

match_index <- match(
  cnv_ids_normalized,
  integration_ids
)

matched <- !is.na(
  match_index
)

cat(
  "CNV cells matched to integration:",
  sum(matched),
  "\n"
)

cat(
  "CNV cells unmatched:",
  sum(!matched),
  "\n"
)


if (
  sum(matched) < 100
) {

  stop(
    "Too few CNV cells matched to integration metadata."
  )

}


# ------------------------------------------------------------
# 13. Build CNV metadata
# ------------------------------------------------------------

cnv_metadata <- integration[
  match_index,
]

rownames(
  cnv_metadata
) <- cnv_cell_names


# ------------------------------------------------------------
# 14. Check CopyKAT categories
# ------------------------------------------------------------

cat("\n")
cat("CopyKAT categories among CNV profiles:\n")

print(
  table(
    cnv_metadata$copykat.pred,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# 15. Extract expression/CNV matrix
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


cat(
  "Numeric CNV matrix:",
  nrow(cnv_values),
  "regions x",
  ncol(cnv_values),
  "cells\n"
)


# ------------------------------------------------------------
# 16. Remove problematic regions
# ------------------------------------------------------------

valid_rows <- is.finite(
  rowSums(
    cnv_values,
    na.rm = TRUE
  )
)

cnv_values <- cnv_values[
  valid_rows,
  ,
  drop = FALSE
]

chrom <- as.character(
  cnv$chrom[
    valid_rows
  ]
)

chrompos <- as.numeric(
  cnv$chrompos[
    valid_rows
  ]
)

abspos <- as.numeric(
  cnv$abspos[
    valid_rows
  ]
)


# ------------------------------------------------------------
# 17. Chromosome normalization
# ------------------------------------------------------------

chrom <- toupper(
  trimws(
    chrom
  )
)

chrom <- gsub(
  "^CHR",
  "",
  chrom
)


# Keep canonical chromosomes.

canonical_chromosomes <- c(
  as.character(1:22),
  "X",
  "Y"
)

keep_chr <- chrom %in%
  canonical_chromosomes

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

abspos <- abspos[
  keep_chr
]


# ------------------------------------------------------------
# 18. Order genomic regions
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

order_index <- order(
  chr_factor,
  chrompos,
  na.last = TRUE
)

cnv_values <- cnv_values[
  order_index,
  ,
  drop = FALSE
]

chrom <- chrom[
  order_index
]

chrompos <- chrompos[
  order_index
]


# ------------------------------------------------------------
# 19. Remove extreme CNV values for visualization
# ------------------------------------------------------------

# This is ONLY for plotting.
# Original CNV values remain untouched.

plot_values <- cnv_values

plot_values[
  plot_values > 2
] <- 2

plot_values[
  plot_values < -2
] <- -2


# ------------------------------------------------------------
# 20. Select aneuploid and diploid cells
# ------------------------------------------------------------

aneuploid_ids <- rownames(
  cnv_metadata
)[
  cnv_metadata$copykat.pred ==
    "aneuploid"
]

diploid_ids <- rownames(
  cnv_metadata
)[
  cnv_metadata$copykat.pred ==
    "diploid"
]

aneuploid_ids <- intersect(
  aneuploid_ids,
  colnames(plot_values)
)

diploid_ids <- intersect(
  diploid_ids,
  colnames(plot_values)
)

cat("\n")
cat(
  "Aneuploid CNV profiles:",
  length(aneuploid_ids),
  "\n"
)

cat(
  "Diploid CNV profiles:",
  length(diploid_ids),
  "\n"
)


# ------------------------------------------------------------
# 21. Randomly sample cells for heatmap
# ------------------------------------------------------------

set.seed(
  228499
)

max_cells_per_group <- 150

aneuploid_plot_ids <- sample(
  aneuploid_ids,
  size = min(
    max_cells_per_group,
    length(aneuploid_ids)
  )
)

diploid_plot_ids <- sample(
  diploid_ids,
  size = min(
    max_cells_per_group,
    length(diploid_ids)
  )
)

heatmap_ids <- c(
  aneuploid_plot_ids,
  diploid_plot_ids
)


# ------------------------------------------------------------
# 22. Reorder heatmap cells
# ------------------------------------------------------------

heatmap_metadata <- cnv_metadata[
  heatmap_ids,
  ,
  drop = FALSE
]

heatmap_metadata$group <- ifelse(
  heatmap_metadata$copykat.pred ==
    "aneuploid",
  "Aneuploid",
  "Diploid"
)

cell_order <- order(
  heatmap_metadata$group
)

heatmap_ids <- rownames(
  heatmap_metadata
)[
  cell_order
]


# ------------------------------------------------------------
# 23. Chromosome boundaries
# ------------------------------------------------------------

chrom_lengths <- table(
  chrom[
    rownames(
      plot_values
    ) %in%
      rownames(
        plot_values
      )
  ]
)

chrom_positions <- rle(
  chrom
)

chrom_ends <- cumsum(
  chrom_positions$lengths
)

chrom_starts <- c(
  1,
  head(
    chrom_ends,
    -1
  ) + 1
)

chrom_centers <- (
  chrom_starts +
  chrom_ends
) / 2

# ------------------------------------------------------------
# 24. CNV heatmap function
# ------------------------------------------------------------

make_cnv_heatmap <- function(
  matrix_data,
  cell_ids,
  metadata,
  output_file,
  plot_title
) {

  cat(
    "\nCreating heatmap:",
    plot_title,
    "\n"
  )

  mat <- matrix_data[
    ,
    cell_ids,
    drop = FALSE
  ]

  cat(
    "Original heatmap matrix:",
    nrow(mat),
    "regions x",
    ncol(mat),
    "cells\n"
  )

  # ----------------------------------------------------------
  # Remove rows with no variation
  # ----------------------------------------------------------

  row_sd <- apply(
    mat,
    1,
    sd,
    na.rm = TRUE
  )

  keep <- is.finite(
    row_sd
  ) &
    row_sd > 0

  mat <- mat[
    keep,
    ,
    drop = FALSE
  ]

  if (
    nrow(mat) < 10
  ) {

    cat(
      "Too few variable CNV regions.\n"
    )

    return(
      invisible(NULL)
    )

  }

  # ----------------------------------------------------------
  # Scale each genomic region
  # ----------------------------------------------------------

  mat_scaled <- t(
    scale(
      t(mat)
    )
  )

  # Replace possible NA/Inf values

  mat_scaled[
    !is.finite(mat_scaled)
  ] <- 0

  # Limit extreme values for visualization only

  mat_scaled[
    mat_scaled > 3
  ] <- 3

  mat_scaled[
    mat_scaled < -3
  ] <- -3


  # ----------------------------------------------------------
  # IMPORTANT:
  # image() expects:
  #
  #   length(x) = nrow(z)
  #   length(y) = ncol(z)
  #
  # Therefore transpose the matrix.
  # ----------------------------------------------------------

  plot_matrix <- t(
    mat_scaled
  )

  cat(
    "Plot matrix:",
    nrow(plot_matrix),
    "cells x",
    ncol(plot_matrix),
    "genomic regions\n"
  )

  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------

  pdf(
    output_file,
    width = 14,
    height = 8
  )

  par(
    mar = c(
      4,
      6,
      4,
      2
    )
  )

  image(
    x = seq_len(
      nrow(plot_matrix)
    ),

    y = seq_len(
      ncol(plot_matrix)
    ),

    z = plot_matrix,

    axes = FALSE,

    xlab = "Cells",

    ylab = "Genomic regions",

    main = plot_title
  )


  # ----------------------------------------------------------
  # Add cell axis
  # ----------------------------------------------------------

  axis(
    1,
    at = seq_len(
      nrow(plot_matrix)
    ),
    labels = FALSE
  )


  # ----------------------------------------------------------
  # Add genomic-region axis
  # ----------------------------------------------------------

  axis(
    2,
    at = seq(
      1,
      ncol(plot_matrix),
      length.out = 10
    ),
    labels = FALSE
  )


  # ----------------------------------------------------------
  # Add chromosome boundaries
  #
  # Since genomic regions are now the X dimension of the
  # original matrix but Y dimension after transpose.
  # ----------------------------------------------------------

  chromosome_positions <- rle(
    chrom
  )

  chromosome_ends <- cumsum(
    chromosome_positions$lengths
  )

  chromosome_starts <- c(
    1,
    head(
      chromosome_ends,
      -1
    ) + 1
  )

  chromosome_centers <- (
    chromosome_starts +
    chromosome_ends
  ) / 2


  abline(
    h = chromosome_ends[
      chromosome_ends <
        ncol(plot_matrix)
    ],
    lty = 3
  )


  # ----------------------------------------------------------
  # Add chromosome labels
  # ----------------------------------------------------------

  axis(
    2,

    at = chromosome_centers[
      chromosome_centers <=
        ncol(plot_matrix)
    ],

    labels = chromosome_positions$values[
      chromosome_centers <=
        ncol(plot_matrix)
    ],

    las = 2,

    cex.axis = 0.6
  )


  box()

  dev.off()

  cat(
    "Saved:",
    output_file,
    "\n"
  )

  invisible(
    mat_scaled
  )
}

# ------------------------------------------------------------
# 25. Create aneuploid vs diploid heatmap
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("CREATING ANEUPLOID VS DIPLOID CNV HEATMAP\n")
cat("====================================================\n")

make_cnv_heatmap(
  matrix_data = plot_values,
  cell_ids = heatmap_ids,
  metadata = heatmap_metadata,
  output_file = file.path(
    figure_dir,
    "CNV_heatmap_aneuploid_vs_diploid.pdf"
  ),
  plot_title =
    "CopyKAT CNV Profiles: Aneuploid vs Diploid"
)


# ------------------------------------------------------------
# 26. Select BC11, BC12 and BC17
# ------------------------------------------------------------

high_cnv_samples <- c(
  "BC11",
  "BC12",
  "BC17"
)

high_cnv_ids <- rownames(
  cnv_metadata
)[
  cnv_metadata$sample_id %in%
    high_cnv_samples &
  cnv_metadata$copykat.pred ==
    "aneuploid"
]

high_cnv_ids <- intersect(
  high_cnv_ids,
  colnames(plot_values)
)


# Sample equal numbers if necessary.

set.seed(
  228499
)

high_cnv_ids <- unlist(
  lapply(
    high_cnv_samples,
    function(s) {

      ids <- rownames(
        cnv_metadata
      )[
        cnv_metadata$sample_id ==
          s &
        cnv_metadata$copykat.pred ==
          "aneuploid"
      ]

      ids <- intersect(
        ids,
        colnames(plot_values)
      )

      if (
        length(ids) == 0
      ) {

        return(
          character(0)
        )

      }

      sample(
        ids,
        size = min(
          100,
          length(ids)
        )
      )

    }
  )
)


cat(
  "BC11/BC12/BC17 CNV profiles:",
  length(high_cnv_ids),
  "\n"
)


# ------------------------------------------------------------
# 27. High-aneuploid sample heatmap
# ------------------------------------------------------------

if (
  length(high_cnv_ids) > 0
) {

  make_cnv_heatmap(
    matrix_data = plot_values,
    cell_ids = high_cnv_ids,
    metadata = cnv_metadata[
      high_cnv_ids,
      ,
      drop = FALSE
    ],
    output_file = file.path(
      figure_dir,
      "CNV_heatmap_BC11_BC12_BC17.pdf"
    ),
    plot_title =
      "CopyKAT CNV Profiles: BC11, BC12 and BC17"
  )

}


# ------------------------------------------------------------
# 28. Calculate mean CNV by sample
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("CALCULATING SAMPLE-LEVEL CNV PROFILES\n")
cat("====================================================\n")


sample_profiles <- list()

samples <- sort(
  unique(
    cnv_metadata$sample_id
  )
)

for (
  s in samples
) {

  ids <- rownames(
    cnv_metadata
  )[
    cnv_metadata$sample_id ==
      s &
    cnv_metadata$copykat.pred ==
      "aneuploid"
  ]

  ids <- intersect(
    ids,
    colnames(cnv_values)
  )

  if (
    length(ids) == 0
  ) {

    next

  }

  sample_profiles[[s]] <-
    rowMeans(
      cnv_values[
        ,
        ids,
        drop = FALSE
      ],
      na.rm = TRUE
    )

}


# ------------------------------------------------------------
# 29. Sample-level profile table
# ------------------------------------------------------------

sample_profile_matrix <- do.call(
  cbind,
  sample_profiles
)

cat(
  "Sample profile dimensions:",
  nrow(sample_profile_matrix),
  "regions x",
  ncol(sample_profile_matrix),
  "samples\n"
)


# ------------------------------------------------------------
# 30. Save sample profile
# ------------------------------------------------------------

saveRDS(
  sample_profile_matrix,
  file.path(
    integration_dir,
    "18_sample_mean_CNV_profiles.rds"
  )
)


# ------------------------------------------------------------
# 31. Calculate mean absolute CNV by sample
# ------------------------------------------------------------

sample_cnv_summary <- data.frame()

for (
  s in colnames(
    sample_profile_matrix
  )
) {

  values <- sample_profile_matrix[
    ,
    s
  ]

  sample_cnv_summary <- rbind(

    sample_cnv_summary,

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

      positive_CNV_fraction =
        mean(
          values > 0.1,
          na.rm = TRUE
        ),

      negative_CNV_fraction =
        mean(
          values < -0.1,
          na.rm = TRUE
        )

    )

  )

}


cat("\nSample-level CNV summary:\n")

print(
  sample_cnv_summary
)


write.csv(
  sample_cnv_summary,
  file.path(
    integration_dir,
    "18_sample_CNV_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 32. Plot sample-level CNV burden
# ------------------------------------------------------------

pdf(
  file.path(
    figure_dir,
    "CNV_profile_by_sample.pdf"
  ),
  width = 10,
  height = 6
)

print(

  ggplot(
    sample_cnv_summary,
    aes(
      x = sample_id,
      y = mean_abs_CNV
    )
  ) +

    geom_col() +

    geom_text(
      aes(
        label =
          round(
            mean_abs_CNV,
            3
          )
      ),
      vjust = -0.4
    ) +

    labs(
      title =
        "Mean Absolute CNV Signal in Aneuploid Cells",
      x = "Sample",
      y = "Mean absolute CNV"
    ) +

    theme_classic()

)

dev.off()



# ------------------------------------------------------------
# 33. Chromosome-level summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("CHROMOSOME-LEVEL CNV SUMMARY\n")
cat("====================================================\n")

# Use aneuploid cells only

all_aneuploid_ids <- intersect(
  aneuploid_ids,
  colnames(cnv_values)
)

if (length(all_aneuploid_ids) == 0) {

  stop(
    "No aneuploid cells available for chromosome-level analysis."
  )

}

aneuploid_mean <- rowMeans(
  cnv_values[
    ,
    all_aneuploid_ids,
    drop = FALSE
  ],
  na.rm = TRUE
)

chromosome_summary <- data.frame(
  chromosome = chrom,
  mean_CNV = aneuploid_mean
) %>%

  group_by(
    chromosome
  ) %>%

  summarise(

    mean_CNV =
      mean(
        mean_CNV,
        na.rm = TRUE
      ),

    mean_abs_CNV =
      mean(
        abs(mean_CNV),
        na.rm = TRUE
      ),

    .groups = "drop"

  )

print(
  chromosome_summary
)

write.csv(
  chromosome_summary,
  file.path(
    integration_dir,
    "18_chromosome_CNV_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 34. Chromosome-level plot
# ------------------------------------------------------------

pdf(
  file.path(
    figure_dir,
    "CNV_chromosome_level_profile.pdf"
  ),
  width = 12,
  height = 6
)

print(

  ggplot(
    chromosome_summary,
    aes(
      x = factor(
        chromosome,
        levels = chr_order
      ),
      y = mean_CNV
    )
  ) +

    geom_col() +

    geom_hline(
      yintercept = 0,
      linetype = "dashed"
    ) +

    labs(
      title =
        "Mean CNV Signal Across Chromosomes",
      x = "Chromosome",
      y = "Mean CNV"
    ) +

    theme_classic()

)

dev.off()


# ------------------------------------------------------------
# 35. Sample CNV interpretation table
# ------------------------------------------------------------

high_cnv_samples <- c(
  "BC11",
  "BC12",
  "BC17"
)

sample_cnv_summary$high_CNV_signal <-
  sample_cnv_summary$sample_id %in%
  high_cnv_samples


write.csv(
  sample_cnv_summary,
  file.path(
    integration_dir,
    "18_CNV_profile_validation_table.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 36. Validation report
# ------------------------------------------------------------

report_file <- file.path(
  integration_dir,
  "18_cnv_profile_validation_summary.txt"
)

report_lines <- c(

  "CNV PROFILE VALIDATION SUMMARY",

  "==============================",

  "",

  paste(
    "Total CopyKAT CNV profiles:",
    ncol(cnv_values)
  ),

  paste(
    "Genomic regions analyzed:",
    nrow(cnv_values)
  ),

  paste(
    "Aneuploid CNV profiles:",
    length(aneuploid_ids)
  ),

  paste(
    "Diploid CNV profiles:",
    length(diploid_ids)
  ),

  "",

  "Samples with substantial CopyKAT aneuploid populations:",

  paste(
    high_cnv_samples,
    collapse = ", "
  ),

  "",

  "Sample-level CNV summary:",

  paste(
    capture.output(
      print(
        sample_cnv_summary
      )
    ),
    collapse = "\n"
  ),

  "",

  "Interpretation:",

  "CopyKAT aneuploidy provides CNV-based evidence.",

  "Coherent chromosome-scale CNV patterns strengthen",

  "the interpretation of aneuploid epithelial cells",

  "as malignant-cell candidates.",

  "CNV evidence is integrated with epithelial and",

  "tumor-candidate identity.",

  "These classifications are research annotations",

  "and are not clinical diagnostic calls."

)

writeLines(
  report_lines,
  report_file
)


# ------------------------------------------------------------
# 37. Final summary
# ------------------------------------------------------------

cat("\n")
cat("====================================================\n")
cat("CNV PROFILE VALIDATION COMPLETED\n")
cat("====================================================\n")

cat(
  "CNV profiles:",
  ncol(cnv_values),
  "\n"
)

cat(
  "Genomic regions:",
  nrow(cnv_values),
  "\n"
)

cat(
  "Aneuploid profiles:",
  length(aneuploid_ids),
  "\n"
)

cat(
  "Diploid profiles:",
  length(diploid_ids),
  "\n"
)

cat("\nResults:\n")

cat(
  file.path(
    integration_dir,
    "18_sample_CNV_summary.csv"
  ),
  "\n"
)

cat(
  file.path(
    integration_dir,
    "18_chromosome_CNV_summary.csv"
  ),
  "\n"
)

cat(
  file.path(
    integration_dir,
    "18_CNV_profile_validation_table.csv"
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
    "CNV_heatmap_aneuploid_vs_diploid.pdf"
  ),
  "\n"
)

cat(
  file.path(
    figure_dir,
    "CNV_heatmap_BC11_BC12_BC17.pdf"
  ),
  "\n"
)

cat(
  file.path(
    figure_dir,
    "CNV_profile_by_sample.pdf"
  ),
  "\n"
)

cat(
  file.path(
    figure_dir,
    "CNV_chromosome_level_profile.pdf"
  ),
  "\n"
)

cat("\n")
cat("====================================================\n")
cat("READY FOR CNV-BASED MALIGNANT CELL INTERPRETATION\n")
cat("====================================================\n")

