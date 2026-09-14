# Our own Splat-inspired scRNA-seq data simulation

args <- commandArgs(trailingOnly = TRUE)
true_expr <- args[1]
outputpath <- args[2]
convert_to_counts <- as.logical(args[3])
max_lib_size <- as.numeric(args[4])
dataname <- args[5]
datatype <- args[6]
dropout_midpoints <- as.numeric(unlist(strsplit(args[7], ",")))

library(Matrix)
library(anndata)

set.seed(110)

# Load true gene expression data
true_expr <- read_h5ad(true_expr)
# Extract the expression matrix
X <- true_expr$X

# Ensure cells x genes
obs_df <- as.data.frame(true_expr$obs)
var_df <- as.data.frame(true_expr$var)
n_obs <- nrow(obs_df)
n_var <- nrow(var_df)
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

add_splat_dropout <- function(true_counts, dropout_mid, dropout_shape) {
  lambda <- true_counts
  lambda[lambda < 0] <- 0

  dropout_prob <- 1 / (1 + exp(-dropout_shape * (log(lambda) - dropout_mid)))
  dropout_prob[lambda == 0] <- 1

  keep_indicator <- matrix(
    rbinom(length(dropout_prob), size = 1, prob = 1 - as.numeric(dropout_prob)),
    nrow = nrow(true_counts),
    ncol = ncol(true_counts)
  )
  
  true_counts * keep_indicator
}

dropout_shape <- -1 # this is a default value

for (dropout_mid in dropout_midpoints) {
  observed_mat <- add_splat_dropout(
    true_counts = X,
    dropout_mid = dropout_mid,
    dropout_shape = dropout_shape
  )

  zero_fraction <- sum(observed_mat == 0) / length(observed_mat)
  zero_pct <- round(100 * zero_fraction, 2)
  zeros <- gsub("\\.", "p", sprintf("%.2f", zero_pct))
  mid_str <- gsub("\\.", "p", sprintf("%.2f", dropout_mid))
  shape_str <- gsub("-", "m", gsub("\\.", "p", sprintf("%.2f", dropout_shape)))

  ann_observed <- AnnData(
    X = Matrix(observed_mat, sparse = TRUE),
    obs = obs_df,
    var = var_df
  )

  write_h5ad(
    ann_observed,
    paste0(
      outputpath,
      "splat_observed_counts_", datatype, "_", dataname, "_2000genes_",
      mid_str, "mid_",
      shape_str, "shape_",
      zeros, "pct0_110seed.h5ad"
    )
  )
}