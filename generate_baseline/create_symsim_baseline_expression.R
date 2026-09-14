# Generate a SymSim-based true expression dataset and save it as .h5ad.

# The script simulates raw counts for a specified number of cells and genes,
# extracts the simulated cell group labels, computes the dropout, and 
# saves the result in .h5ad format.

# Read simulation settings from the command line
args <- commandArgs(trailingOnly = TRUE)
ncells     = as.numeric(as.character(args[1]))
ngenes     = as.numeric(as.character(args[2]))
seed       = as.numeric(as.character(args[3]))
outputpath = args[4]

library(anndata)
library(Matrix)
library(SymSim)
# Define the phylogenetic tree used by SymSim
phyla <- Phyla5()

symsim_simulatedata = function(ncells=1000,ngenes=2000,seed=123,outputpath=outputpath){
  
  ngroups = length(phyla$tip.label)
  
  # Generate true counts
  true_counts_res = SimulateTrueCounts(ncells_total=ncells, 
                                       min_popsize=floor(ncells/ngroups), 
                                       i_minpop=2, 
                                       ngenes=ngenes, 
                                       nevf=40,
                                       evf_type="discrete", 
                                       n_de_evf=20, 
                                       vary="s", 
                                       Sigma=0.5, 
                                       phyla=phyla, 
                                       randseed=seed, 
                                       gene_effects_sd=1, 
                                       gene_effect_prob = 0.2)
  
  # Convert the simulated counts to a data frame for zero-fraction calculation
  true_count_data = data.frame(t(true_counts_res$counts))
  labels = true_counts_res$cell_meta$pop

  # Compute the fraction of zero entries
  zero_fraction <- sum(true_count_data == 0) / (nrow(true_count_data) * ncol(true_count_data))
  zero_pct <- round(100 * zero_fraction, 2)
  zeros <- gsub("\\.", "p", sprintf("%.2f", zero_pct))

  # The AnnData object
  ann_true <- AnnData(
    X = Matrix(t(true_counts_res$counts),sparse=T),
    obs = data.frame(group = labels,
                     row.names = rownames(true_count_data)),
    var = data.frame(col.names = colnames(true_count_data))
  )
  # Save the SymSim true expression dataset as .h5ad
  write_h5ad(ann_true, paste0(outputpath,"symsim_true_counts_",ngenes,"genes_",ncells,"cells_",zeros,"pct0_",seed,"seed.h5ad"))
}


symsim_simulatedata(ncells=ncells,ngenes=ngenes,seed=seed,outputpath=outputpath)