# Our own SPARSim-like scRNA-seq data simulation

args <- commandArgs(trailingOnly = TRUE)
true_expr <- args[1]
outputpath <- args[2]
convert_to_counts <- as.logical(args[3])
max_lib_size <- as.numeric(args[4])
dataname <- args[5]
datatype <- args[6]
library_size <- as.numeric(unlist(strsplit(args[7], ",")))

library(Matrix)
library(anndata)
library(extraDistr) # for rmvhyper

set.seed(110)

true_expr <- read_h5ad(true_expr)
X <- true_expr$X

# Ensure cells x genes
n_obs <- nrow(as.data.frame(true_expr$obs))
n_var <- nrow(as.data.frame(true_expr$var))
if (nrow(X) == n_var && ncol(X) == n_obs) X <- t(X)

# Convert SPARSim intensities to integer pseudo-counts
if (convert_to_counts) {
  input_fragment_lib_size <- max_lib_size * 1e2
  digits <- floor(log10(input_fragment_lib_size)) + 1
  new_libsize <- 10^digits
  cell_sums <- Matrix::rowSums(X)
  cell_sums[cell_sums == 0] <- 1
  X_int <- round(X * (new_libsize / cell_sums))
  X_int[X_int < 0] <- 0
  X <- X_int
}

# Number of cells
n_cells <- nrow(X)

# For each library size
for (l in library_size) {
  # Create a list to store sampled columns for each cell
  sampled_cols <- vector("list", n_cells)

  # For each cell
  for (cell_idx in seq_len(n_cells)) {
    # Get the counts for the current cell
    # counts_vec <- as.numeric(X_int[cell_idx, ])
    counts_vec <- as.numeric(X[cell_idx, ])
    # Sample counts using multivariate hypergeometric distribution
    total <- sum(counts_vec)
    k <- min(l, total) # number of draws cannot exceed total counts
    sampled_cols[[cell_idx]] <- as.integer(rmvhyper(1, counts_vec, k))
  }
  # Combine sampled columns into a matrix
  sampled_mat <- do.call(rbind, sampled_cols)

  # Calculate sparsity
  zero_fraction <- sum(sampled_mat == 0) / length(sampled_mat)
  zero_pct <- round(100 * zero_fraction, 2)
  zeros <- gsub("\\.", "p", sprintf("%.2f", zero_pct))

  # Create AnnData object for observed counts
  ann_observed <- AnnData(
    X = Matrix(sampled_mat, sparse = TRUE),
    obs = as.data.frame(true_expr$obs),
    var = as.data.frame(true_expr$var)
  )

  write_h5ad(ann_observed, paste0(outputpath, "sparsim_observed_counts_", datatype, "_", dataname, "_2000genes_", l, "library_", zeros, "pct0_110seed.h5ad"))
}