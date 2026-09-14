# Pavel et al. scRNA-seq technical noise simulation

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
true_expr <- args[1]
outputpath <- args[2]
convert_to_counts <- as.logical(args[3])
max_lib_size <- as.numeric(args[4])
dataname <- args[5]
datatype <- args[6]
dropout_rates <- as.numeric(unlist(strsplit(args[7], ",")))

library(Matrix)
library(anndata)

set.seed(110)

true_expr <- read_h5ad(true_expr)
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

# Labels
labels <- as.character(obs_df$group)

# Alisa's function to add dropouts
add_dropouts_geneexpression <- function(exp_mat, labels, rate = 0.1) {
    exp_out <- matrix(0, nrow = nrow(exp_mat), ncol = ncol(exp_mat))
    rownames(exp_out) <- rownames(exp_mat)
    colnames(exp_out) <- colnames(exp_mat)

    for (cl in unique(labels)) {
        cat("class", cl, "\n")

        idx <- which(labels == cl)
        exp2 <- exp_mat[idx, , drop = FALSE]

        # normalize each gene within the current class across cells
        gene_min <- apply(exp2, 2, min)
        gene_max <- apply(exp2, 2, max)
        gene_range <- gene_max - gene_min
        gene_range[gene_range == 0] <- 1

        normalized_df <- sweep(exp2, 2, gene_min, "-")
        normalized_df <- sweep(normalized_df, 2, gene_range, "/")
        normalized_df <- (1 - normalized_df) * rate

        cnt <- 0

        for (c in seq_len(nrow(exp2))) {
            t <- normalized_df[c, ]
            tv <- exp2[c, ]

            tgenes <- numeric(length(tv))

            for (g in seq_along(t)) {
                f <- runif(1)
                valt <- t[g]
                val <- tv[g]

                if (f <= valt) {
                    val <- 0
                }

                tgenes[g] <- abs(val)
            }

            exp_out[idx[c], ] <- tgenes
            cnt <- cnt + 1

            if (cnt %% 100 == 0) {
                cat(cnt, nrow(exp2), "\n")
            }
        }
    }

    exp_out
}

for (rate in dropout_rates) {
  observed_mat <- add_dropouts_geneexpression(X, labels, rate = rate)

  zero_fraction <- sum(observed_mat == 0) / length(observed_mat)
  zero_pct <- round(100 * zero_fraction, 2)
  zeros <- gsub("\\.", "p", sprintf("%.2f", zero_pct))

  rate2 <- gsub("\\.", "p", sprintf("%.2f", rate))

  ann_observed <- AnnData(
    X = Matrix(observed_mat, sparse = TRUE),
    obs = obs_df,
    var = var_df
  )

  write_h5ad(
    ann_observed,
    paste0(
      outputpath,
      "dropout_observed_counts_", datatype, "_", dataname, "_2000genes_",
      rate2, "rate_",
      zeros, "pct0_110seed.h5ad"
    )
  )
}
