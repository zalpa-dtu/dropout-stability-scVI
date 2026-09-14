# Our own ESCO-inspired scRNA-seq technical noise simulation (similar to Splat)

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

add_esco_dropout <- function(true_counts, dropout_mid, dropout_shape) {
  lambda <- true_counts
  lambda[lambda < 0] <- 0

  # First scale every cell to the median library size
  # This corresponds to cell.normmeans in the ESCO code in escoSimZeroInflate:
  # median(colSums(cell.means)) * t(t(cell.means) / colSums(cell.means))
  # Here the matrix is cells x genes, so we use rowSums instead of colSums
  # https://github.com/JINJINT/ESCO/blob/master/R/esco-simulate.R#L1165
  cell_sums <- Matrix::rowSums(lambda)
  median_cell_sum <- median(cell_sums)

  if (median_cell_sum == 0 || is.na(median_cell_sum)) {
    stop("Median cell sum is zero; cannot normalize ESCO dropout means.")
  }

  size_factors <- cell_sums / median_cell_sum
  size_factors[size_factors == 0] <- 1
  norm_means <- Matrix::Diagonal(x = 1 / size_factors) %*% lambda

  # Then apply the ESCO dropout formula to get the dropout probabilities (same as Splat)
  dropout_prob <- 1 / (1 + exp(-dropout_shape * (log(norm_means) - dropout_mid)))
  dropout_prob[norm_means == 0] <- 1

  # ESCO correction for the Splat bias:

  # This is the term Pr{Y_hat_gc = 0} in equation (13) in https://pmc.ncbi.nlm.nih.gov/articles/instance/8388018/bin/btab116_supplementary_data.pdf
  # Note dpois(0, lambda) = exp(-lambda)
  poisson_zero_prob <- exp(-norm_means)
  # Probability that the simulated true count is nonzero:
  # The denominator in equation (13) in https://pmc.ncbi.nlm.nih.gov/articles/instance/8388018/bin/btab116_supplementary_data.pdf
  poisson_nonzero_prob <- 1 - poisson_zero_prob

  keep_prob <- (1 - dropout_prob) / poisson_nonzero_prob
  keep_prob[keep_prob > 1] <- 1
  keep_prob[keep_prob < 0] <- 0
  keep_prob[is.na(keep_prob)] <- 1

  keep_indicator <- matrix(
    rbinom(length(keep_prob), size = 1, prob = as.numeric(keep_prob)),
    nrow = nrow(true_counts),
    ncol = ncol(true_counts)
  )
  
  # return the observed counts after dropout
  true_counts * keep_indicator
}

# ESCO dropout parameters
dropout_shape <- -1 # this is a default value

for (dropout_mid in dropout_midpoints) {
  observed_mat <- add_esco_dropout(
    true_counts = X,
    dropout_mid = dropout_mid,
    dropout_shape = dropout_shape
  )

  # Calculate the percentage of zeros in the observed matrix
  zero_fraction <- sum(observed_mat == 0) / length(observed_mat)
  zero_pct <- round(100 * zero_fraction, 2)
  zeros <- gsub("\\.", "p", sprintf("%.2f", zero_pct))
  mid_str <- gsub("\\.", "p", sprintf("%.2f", dropout_mid))
  shape_str <- gsub("-", "m", gsub("\\.", "p", sprintf("%.2f", dropout_shape)))

  # Create an AnnData object for the observed data
  ann_observed <- AnnData(
    X = Matrix(observed_mat, sparse = TRUE),
    obs = obs_df,
    var = var_df
  )

  # Save the observed data to an .h5ad file
  write_h5ad(
    ann_observed,
    paste0(
      outputpath,
      "esco_observed_counts_", datatype, "_", dataname, "_2000genes_",
      mid_str, "mid_",
      shape_str, "shape_",
      zeros, "pct0_110seed.h5ad"
    )
  )
}