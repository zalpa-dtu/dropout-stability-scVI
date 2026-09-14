# Generate a SPARSim-based true expression dataset and save it as .h5ad.

# This script uses a real single-cell reference dataset to estimate SPARSim
# simulation parameters, then generates a synthetic true gene expression matrix. 

# Read simulation settings from the command line
args <- commandArgs(trailingOnly = TRUE)
datatype    = args[1] # reference dataset type, currently expected to be "TS"
dataname    = args[2] # dataset name used in the output filename
datapath    = args[3] # path to the input Seurat .rds object
ngenes      = as.numeric(as.character(args[4])) # number of highly variable genes to use for simulation
seed       = as.numeric(as.character(args[5])) # random seed for reproducibility
outputpath  = args[6] # directory where the .h5ad output file should be written

library(SPARSim)
library(Seurat)
library(Matrix)
library(anndata)

# TS - Tabula Sapiens
if(datatype == 'TS'){
  # Read rds data in R
  data_seurat_raw <- readRDS(datapath)
  label_name = 'cell_ontology_class'
  
  # Create proportions table
  props <- prop.table(table(data_seurat_raw[[label_name]][,1]))
  # Keep groups with at least 2% representation
  keep_groups <- names(props[props >= 0.02])
  
  # Subset data
  data_seurat <- subset(data_seurat_raw, subset = !!as.name(label_name) %in% keep_groups)
  
  # Find variable features - top ngenes that are highly variable
  data_seurat_ <- FindVariableFeatures(
    data_seurat,
    nfeatures = ngenes,
    selection.method = 'vst',
    assay = 'RNA',
    layer = 'log_normalized' 
  )
  
  # Get labels and experimental conditions *(cell type)*
  labels <- droplevels(data_seurat@meta.data$cell_ontology_class)
  labels <- as.character(labels)
  experimental_condition = split(seq_along(labels), labels) # this creates a list of cell indices for each cell type

  # Get HVG and expression data
  hvg <- VariableFeatures(data_seurat_) # this return a character vector of gene names
  expr_hvg_raw <- as.matrix(GetAssayData(data_seurat_, layer = "decontXcounts")[hvg, ])
  expr_hvg_norm <- as.matrix(exp(GetAssayData(data_seurat_, layer = "log_normalized")[hvg, ])-1) # SPARSim assumes non-logaritmized data (usually scran normalization which is non-log https://www.sc-best-practices.org/preprocessing_visualization/normalization.html)
}

# Retrieve old cell names
old_cellnames = Cells(data_seurat)

# Modify colnames to include group information
colnames(expr_hvg_raw) = paste0(labels,'/',old_cellnames)
colnames(expr_hvg_norm) = paste0(labels,'/',old_cellnames)

# Estimate parameters for simulation
sparsim_param = SPARSim_estimate_parameter_from_data(raw_data = expr_hvg_raw,
                                                     norm_data = expr_hvg_norm,
                                                     conditions = experimental_condition)

simulate_data_sparsim = function(params,
                                 outputpath,
                                 ngenes,
                                 seed,
                                 datatype,
                                 dataname){
  
  params_local <- params
  
    
  # Simulate data
  sim_result = SPARSim_simulation(dataset_parameter = params_local,
                                  gene_expr_simulation_seed = seed,
                                  count_data_simulation_seed = seed + 1)
  
  # New labels derived from cell names
  # Take part before '/' as group label
  labels_new = sub("/.*", "", colnames(sim_result$gene_matrix))

  # Compute the fraction of zero entries for the output filename
  gene_zero_fraction <- sum(sim_result$gene_matrix == 0) / length(sim_result$gene_matrix)
  gene_zero_pct <- round(100 * gene_zero_fraction, 2)
  gene_zeros <- gsub("\\.", "p", sprintf("%.2f", gene_zero_pct))
  
  # The AnnData object
  ann_true <- AnnData(
    X = Matrix(t(sim_result$gene_matrix), sparse = TRUE),
    obs = data.frame(group = labels_new,
                      row.names = rownames(t(sim_result$gene_matrix))),
    var = data.frame(row.names = colnames(t(sim_result$gene_matrix)))
  )
  # Save the SPARSim true expression dataset as .h5ad
  write_h5ad(ann_true, paste0(outputpath, "sparsim_true_gene_expr_", datatype, '_', dataname, '_', ngenes, "genes_", gene_zeros, "pct0_", seed, 'seed.h5ad'))
}

simulate_data_sparsim(sparsim_param,
                    outputpath,
                    ngenes,
                    seed,
                    datatype,
                    dataname)