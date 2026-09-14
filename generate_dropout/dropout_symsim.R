# Our own SymSim-like scRNA-seq data simulation

args <- commandArgs(trailingOnly = TRUE)
true_expr <- args[1]
outputpath <- args[2]
convert_to_counts <- as.logical(args[3])
max_lib_size <- as.numeric(args[4])
dataname <- args[5]
datatype <- args[6]
depth_mean <- as.numeric(unlist(strsplit(args[7], ",")))

library(Matrix)
library(anndata)
library(SymSim)

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

# SymSim expects genes x cells
true_counts <- t(X)

# Cell metadata for True2ObservedCounts
labels <- obs_df$group
# Create a data frame with the population labels
cell_meta <- data.frame(pop = labels)
# Set row names of cell_meta to match the cell columns of true_counts
rownames(cell_meta) <- colnames(true_counts)

# Estimate gene length from dataset
data(gene_len_pool)
gene_len = sample(gene_len_pool, nrow(true_counts), replace = FALSE)

print ("Simulating observed counts with SymSim noise...")
# Add SymSim noise
for(dm in depth_mean){

    observed_counts_res = True2ObservedCounts(true_counts,
                                            cell_meta,
                                            protocol="UMI",
                                            gene_len=gene_len,
                                            depth_mean = dm,
                                            depth_sd = 1e3)

    print (paste0("Done simulating observed counts with depth mean ", dm))
    # Calculate sparsity
    observed_mat <- t(observed_counts_res$counts)
    zero_fraction <- sum(observed_mat == 0) / length(observed_mat)
    zero_pct <- round(100 * zero_fraction, 2)
    zeros <- gsub("\\.", "p", sprintf("%.2f", zero_pct))

    # Create AnnData object for observed counts
    ann_observed <- AnnData(
        X = Matrix(observed_mat, sparse = TRUE),
        obs = obs_df,
        var = var_df
    )
    write_h5ad(ann_observed, paste0(outputpath, "symsim_observed_counts_", datatype, "_", dataname, "_2000genes_", dm,"depth_", zeros, "pct0_110seed.h5ad"))
}