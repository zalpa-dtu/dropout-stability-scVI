# Generate a Splatter-based true expression dataset and save it as .h5ad.

# The script simulates grouped single-cell counts, extracts the true count matrix and group labels, 
# computes the dropout, and saves the result in .h5ad format.

# Read simulation settings from the command line
args <- commandArgs(trailingOnly = TRUE)
ncells      = as.numeric(as.character(args[1])) # number of cells to simulate
ngenes      = as.numeric(as.character(args[2])) # number of highly variable genes to use for simulation
seed       = as.numeric(as.character(args[3])) # random seed for reproducibility
outputpath  = args[4] # directory where the .h5ad output file should be written

library(splatter)
library(SingleCellExperiment)
library(anndata)
library(Matrix)

# Fix the random seed for reproducibility
set.seed(seed)

# Simulate grouped single-cell expression data with splatter
sim <- splatSimulate(
  nGenes = ngenes,
  batchCells = ncells,
  group.prob = rep(1 / 8, 8),
  method = "groups",
  dropout.type = "experiment",
  sparsify = FALSE,
  verbose = FALSE
)

assayNames(sim)

# In SingleCellExperiment assays are stored as genes x cells
true_counts <- assay(sim, "TrueCounts")
labels <- as.character(colData(sim)$Group)

# Compute the fraction of zero entries for the output filename
zero_fraction <- sum(true_counts == 0) / length(true_counts)
zero_pct <- round(100 * zero_fraction, 2)
zeros <- gsub("\\.", "p", sprintf("%.2f", zero_pct))

# The AnnData object
ann_true <- AnnData(
  X = Matrix(t(true_counts), sparse = TRUE),
  obs = data.frame(group = labels,
                   row.names = colnames(true_counts)),
  var = data.frame(gene = rownames(true_counts),
                   row.names = rownames(true_counts))
)

# Save the Splat true expression dataset as .h5ad
write_h5ad(
  ann_true,
  paste0(outputpath, "splat_true_counts_", ngenes, "genes_", ncells, "cells_", zeros, "pct0_", seed, "seed.h5ad")
)