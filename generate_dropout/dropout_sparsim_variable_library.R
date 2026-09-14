# Our own SPARSim-like scRNA-seq data simulation with variable library sizes

args <- commandArgs(trailingOnly = TRUE)
true_expr <- args[1]
outputpath <- args[2]
convert_to_counts <- as.logical(args[3])
max_lib_size <- as.numeric(args[4])
dataname <- args[5]
datatype <- args[6]
library_mean <- as.numeric(unlist(strsplit(args[7], ",")))

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

library_variance <- (0.2 * library_mean)^2
min_library_size <- 1

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
for (i in seq_along(library_mean)) {
  mean_lib <- library_mean[i]
  var_lib <- library_variance[i]

  # Create a list to store sampled columns for each cell
  sampled_cols <- vector("list", n_cells)

  cell_library_sizes <- integer(n_cells)

  # For each cell
  for (cell_idx in seq_len(n_cells)) {
    # Get the counts for the current cell
    counts_vec <- as.numeric(X[cell_idx, ])

    # Sample counts using multivariate hypergeometric distribution
    total <- sum(counts_vec)

    # Add some variability to the library size for each cell
    cell_library_size <- round(rnorm(
      n = 1,
      mean = mean_lib,
      sd = sqrt(var_lib)
    ))

    # Ensure library size is at least the minimum to avoid issues with rmvhyper
    cell_library_size <- max(cell_library_size, min_library_size)
    k <- as.integer(min(cell_library_size, total)) # number of draws cannot exceed total counts
    cell_library_sizes[cell_idx] <- k
    sampled_cols[[cell_idx]] <- as.integer(rmvhyper(1, counts_vec, k))
  }
  # Combine sampled columns into a matrix
  sampled_mat <- do.call(rbind, sampled_cols)

  # Calculate sparsity
  zero_fraction <- sum(sampled_mat == 0) / length(sampled_mat)
  zero_pct <- round(100 * zero_fraction, 2)
  zeros <- gsub("\\.", "p", sprintf("%.2f", zero_pct))

  # Make sure the is no "." in the variable library name for file naming
  var_lib_label <- gsub("\\.", "p", as.character(var_lib))

  obs_df <- as.data.frame(true_expr$obs)
  obs_df$sampled_library_size <- cell_library_sizes

  # Create AnnData object for observed counts
  ann_observed <- AnnData(
    X = Matrix(sampled_mat, sparse = TRUE),
    obs = obs_df,
    var = as.data.frame(true_expr$var)
  )

  write_h5ad(ann_observed, paste0(outputpath, "sparsim_observed_counts_", datatype, "_", dataname, "_2000genes_", mean_lib, "meanlibrary_", var_lib_label, "varlibrary_", zeros, "pct0_110seed.h5ad"))
}