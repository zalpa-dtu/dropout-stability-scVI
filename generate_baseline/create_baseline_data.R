# Create a real single-cell reference dataset in .h5ad format.

# The script loads a Seurat object from an RDS file, filters out rare cell groups,
# selects highly variable genes, extracts the corresponding raw count matrix,
# and saves the processed data together with cell annotations as an .h5ad file.

# Read the settings from the command line
args <- commandArgs(trailingOnly = TRUE)
datatype    = args[1] # reference dataset type
dataname    = args[2] # dataset name used in the output filename
datapath    = args[3] # path to the input Seurat .rds object
ngenes      = as.numeric(as.character(args[4])) # number of highly variable genes to use for simulation
seed       = as.numeric(as.character(args[5])) # random seed for reproducibility
outputpath  = args[6] # directory where the .h5ad output file should be written

library(Matrix)
library(anndata)
library(Seurat)

set.seed(seed)

data_seurat_raw <- readRDS(datapath)

# Choose the annotation field used as the grouping label and the assay based on the dataset type
if (datatype == "Other") {
  label_name <- "cell_line"
  assay = "originalexp"
} else if (datatype %in% c("TS", "TM")) {
  label_name <- "cell_ontology_class"
  assay = "RNA"
} else {
  stop("Unsupported datatype. Please use 'Other', 'TS', or 'TM'.")
}

if (datatype %in% c("Other", "TM")) {
  layer = "data"
  layer2 = "counts"
} else if (datatype == "TS") {
  layer = "log_normalized"
  layer2 = "decontXcounts"
}

# Keep only cell groups that make up at least 2% of the dataset
props <- prop.table(table(data_seurat_raw[[label_name]][, 1]))
keep_groups <- names(props[props >= 0.02])

# Filter out rare groups before feature selection
data_seurat <- subset(
  data_seurat_raw,
  subset = !!as.name(label_name) %in% keep_groups
)

# Select the top highly variable genes for downstream experiments
data_seurat <- FindVariableFeatures(
  data_seurat,
  nfeatures = ngenes,
  selection.method = "vst",
  assay = assay,
  layer = layer
)

# Extract and inspect group labels
labels <- as.character(data_seurat@meta.data[[label_name]])
cat("Number of groups:", length(unique(labels)), "\n")
print(table(labels))

hvg <- VariableFeatures(data_seurat)

# Extract the raw count matrix for the selected genes (cells x genes)
X <- t(as.matrix(GetAssayData(data_seurat, assay = assay, layer = layer2)[hvg, ]))

# Store cell metadata and gene metadata
obs_df <- data_seurat@meta.data
obs_df$group <- labels
var_df <- data.frame(row.names = hvg)

# Compute the fraction of zero entries for the output filename
true_zero_fraction <- sum(X == 0) / length(X)
true_zero_pct <- round(100 * true_zero_fraction, 2)
true_zeros <- gsub("\\.", "p", sprintf("%.2f", true_zero_pct))

# The AnnData object
ann_true <- AnnData(
  X = Matrix(X, sparse = TRUE),
  obs = obs_df,
  var = var_df
)

# Save the processed real expression dataset as .h5ad
write_h5ad(
  ann_true,
  paste0(
    outputpath,
    "real_true_gene_expr_", datatype, "_", dataname, "_",
    ngenes, "genes_",
    true_zeros, "pct0_", seed, "seed.h5ad"
  )
)
