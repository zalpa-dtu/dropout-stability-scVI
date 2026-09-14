# Random scRNA-seq technical noise simulation

args <- commandArgs(trailingOnly = TRUE)
true_expr <- args[1]
outputpath <- args[2]
convert_to_counts <- as.logical(args[3])
max_lib_size <- as.numeric(args[4])
dataname <- args[5]
datatype <- args[6]
propzeros <- as.numeric(unlist(strsplit(args[7], ",")))

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

add_random_dropout <- function(true_counts, propZeroes) {
    zeroId <- matrix(1, nrow = nrow(true_counts), ncol = ncol(true_counts))

    samp <- sample(
        1:length(zeroId),
        floor(length(zeroId) * propZeroes)
    )

    zeroId[samp] <- 0

    # Do not label already-zero entries as newly introduced dropout
    tc_mat <- as.matrix(true_counts)
    zeroId[tc_mat == 0] <- 1
    samp <- samp[!samp %in% which(tc_mat == 0)]

    counts_dropout <- tc_mat * zeroId

    return(list(
        counts = counts_dropout,
        dropout_mask = zeroId,
        added_dropout_positions = samp
    ))
}

for (prop in propzeros) {
    observed_res <- add_random_dropout(
        true_counts = X,
        propZeroes = prop
    )

    # Calculate the percentage of zeros in the observed matrix
    zero_fraction <- sum(as.matrix(observed_res$counts) == 0) / length(observed_res$counts)
    zero_pct <- round(100 * zero_fraction, 2)
    zeros <- gsub("\\.", "p", sprintf("%.2f", zero_pct))
    prop_str <- gsub("\\.", "p", sprintf("%.2f", prop))

    # Create an AnnData object for the observed data
    ann_observed <- AnnData(
        X = Matrix(observed_res$counts, sparse = TRUE),
        obs = obs_df,
        var = var_df
    )

    # Save the observed data to an .h5ad file
    write_h5ad(
        ann_observed,
        paste0(
        outputpath,
        "random_observed_counts_", datatype, "_", dataname, "_2000genes_",
        prop_str, "prop_",
        zeros, "pct0_110seed.h5ad"
        )
    )
}