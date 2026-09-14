# Our own ZINB-WaVE-inspired scRNA-seq data simulation

args <- commandArgs(trailingOnly = TRUE)
true_expr <- args[1]
outputpath <- args[2]
convert_to_counts <- as.logical(args[3])
max_lib_size <- as.numeric(args[4])
dataname <- args[5]
datatype <- args[6]
dropout_alpha0 <- as.numeric(unlist(strsplit(args[7], ",")))

library(Matrix)
library(anndata)
library(zinbwave)
library(matrixStats)
library(magrittr)
library(ggplot2)
library(sparseMatrixStats)

# Register BiocParallel Serial Execution
BiocParallel::register(BiocParallel::SerialParam())

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

inv_logit <- function(x) 1 / (1 + exp(-x))

n <- nrow(X)
J <- ncol(X)
M <- 1
L <- 1
K <- 2 # We chose this K based on this: https://bioconductor.org/packages/release/bioc/vignettes/zinbwave/inst/doc/intro.html

scale_or_zero <- function(x) {
  x <- as.numeric(x)
  if (sd(x) == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

# Cell-level covariate:
X_cell <- matrix(scale_or_zero(log1p(Matrix::rowSums(X))), ncol = M) # n x M

# Gene-level covariate:
V_gene <- matrix(scale_or_zero(log1p(Matrix::colMeans(X))), ncol = L) # J x L

# Offset
O_pi <- matrix(0, nrow = n, ncol = J)

X_gc_counts <- t(X) # genes x cells

se <- SummarizedExperiment(
  assays = list(counts = X_gc_counts)
)

zinb <- zinbFit(se, K=2,  X = X_cell, V = V_gene, epsilon=nrow(se))

beta_pi  <- zinbwave::getBeta_pi(zinb)
gamma_pi <- zinbwave::getGamma_pi(zinb)
# Unknown cell-level covariates
alpha_pi <- zinbwave::getAlpha_pi(zinb)
W        <- zinbwave::getW(zinb)

# Checks
cat("X:", dim(X), "\n")
cat("X_gc_counts:", dim(X_gc_counts), "\n")
cat("X_cell:", dim(X_cell), "\n")
cat("beta_pi:", dim(beta_pi), "\n")
cat("V_gene:", dim(V_gene), "\n")
cat("gamma_pi:", dim(gamma_pi), "\n")
cat("W:", dim(W), "\n")
cat("alpha_pi:", dim(alpha_pi), "\n")
cat("O_pi:", dim(O_pi), "\n")

# Base structured part of logit(pi_ij)
eta_base <- X_cell %*% beta_pi + t(V_gene %*% gamma_pi) + W %*% alpha_pi + O_pi

for (a0 in dropout_alpha0) {
  eta <- a0 + eta_base
  pi_mat <- inv_logit(eta)

  dropout_mask <- matrix(
    rbinom(length(pi_mat), 1, c(pi_mat)),
    nrow = n,
    ncol = J
  )

  observed_mat <- X * (1 - dropout_mask)

  zero_fraction <- sum(observed_mat == 0) / length(observed_mat)
  zero_pct <- round(100 * zero_fraction, 2)
  zeros <- gsub("\\.", "p", sprintf("%.2f", zero_pct))
  a0_str <- gsub("\\.", "p", as.character(a0))

  ann_observed <- AnnData(
    X = Matrix(observed_mat, sparse = TRUE),
    obs = obs_df,
    var = var_df
  )

  write_h5ad(
    ann_observed,
    paste0(
      outputpath,
      "zinbwave_observed_counts_", datatype, "_", dataname, "_2000genes_",
      a0_str, "alpha0_",
      zeros, "pct0_110seed.h5ad"
    )
  )
}
