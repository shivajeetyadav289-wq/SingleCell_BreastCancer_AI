library(Seurat)

project_dir <- "~/Projects/SingleCell_BreastCancer_AI"

input_file <- file.path(
  project_dir,
  "data",
  "processed",
  "GSE228499_singlets_clean.rds"
)

output_file <- file.path(
  project_dir,
  "data",
  "processed",
  "GSE228499_integration_prepared.rds"
)

cat("\nLoading clean singlet object...\n")

old_obj <- readRDS(input_file)

cat("Dimensions:\n")
print(dim(old_obj))

cat("\nCells per sample:\n")
print(table(old_obj$sample_id))


# ------------------------------------------------------------
# Join existing Seurat layers
# ------------------------------------------------------------

cat("\nJoining existing layers...\n")

old_obj <- JoinLayers(
  old_obj,
  assay = "RNA"
)

cat("\nLayers after JoinLayers:\n")
print(Layers(old_obj[["RNA"]]))


# ------------------------------------------------------------
# Extract ONLY counts
# ------------------------------------------------------------

cat("\nExtracting counts...\n")

counts <- LayerData(
  old_obj,
  assay = "RNA",
  layer = "counts"
)

cat("Counts dimensions:\n")
print(dim(counts))


# ------------------------------------------------------------
# Extract metadata
# ------------------------------------------------------------

metadata <- old_obj@meta.data

metadata <- metadata[
  colnames(counts),
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# Create NEW Seurat object
# ------------------------------------------------------------

cat("\nCreating fresh Seurat object...\n")

obj <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata,
  project = "GSE228499"
)

cat("\nFresh RNA layers:\n")
print(Layers(obj[["RNA"]]))

cat("\nFresh dimensions:\n")
print(dim(obj))


# ------------------------------------------------------------
# Split by sample
# ------------------------------------------------------------

cat("\nSplitting by sample_id...\n")

obj[["RNA"]] <- split(
  obj[["RNA"]],
  f = obj$sample_id
)

cat("\nLayers after splitting:\n")
print(Layers(obj[["RNA"]]))


# ------------------------------------------------------------
# Normalize
# ------------------------------------------------------------

cat("\nNormalizing...\n")

obj <- NormalizeData(
  obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)


# ------------------------------------------------------------
# Variable features
# ------------------------------------------------------------

cat("\nFinding variable features...\n")

obj <- FindVariableFeatures(
  obj,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = TRUE
)


# ------------------------------------------------------------
# Scale
# ------------------------------------------------------------

cat("\nScaling...\n")

obj <- ScaleData(
  obj,
  features = VariableFeatures(obj),
  verbose = TRUE
)


# ------------------------------------------------------------
# PCA
# ------------------------------------------------------------

cat("\nRunning PCA...\n")

obj <- RunPCA(
  obj,
  features = VariableFeatures(obj),
  npcs = 30,
  verbose = TRUE
)


# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

saveRDS(
  obj,
  output_file
)


cat("\n")
cat("====================================================\n")
cat("INTEGRATION PREPARATION COMPLETED\n")
cat("====================================================\n")

cat("Genes:", nrow(obj), "\n")
cat("Cells:", ncol(obj), "\n")
cat("Samples:", length(unique(obj$sample_id)), "\n")
cat("HVGs:", length(VariableFeatures(obj)), "\n")

cat("\nCells per sample:\n")
print(table(obj$sample_id))

cat("\nFinal RNA layers:\n")
print(Layers(obj[["RNA"]]))

cat("\nSaved:\n")
cat(output_file, "\n")

cat("\nREADY FOR RPCA INTEGRATION\n")