# Cluster Accuracy and Stability Under Dropouts through Variational Inference in scRNA-seq Analysis
We want to compare embedding approaches for single-cell gene expression data using real and synthetic true expression sources with controlled noise generation. Current literature mostly uses PCA (maybe for good reason), we want to try if we get a better embedding using deep learning approaches like VAEs. Most popular (and basic) is the [scVI model](https://docs.scvi-tools.org/en/stable/api/reference/scvi.model.SCVI.html#scvi.model.SCVI).

## Environment setup

This project uses one Conda environment for Python and `renv` for R packages.

### Python

Create the Python environment from `environment.yml`:

```bash
conda env create --prefix ./.conda-env -f environment.yml
```

Activate it from the repository root:

```bash
conda activate "$PWD/.conda-env"
```

Register it as a Jupyter kernel:

```bash
python -m ipykernel install --user --name dropout-project --display-name "Python (dropout-project)"
```

### R

This project uses `renv` to record the R package environment. The lockfile was created with R 4.5.2.

Install R 4.5.2 or a close compatible version, then restore the R packages from the lockfile.

Start R from the repository root:

```bash
R
```

Then inside R, run:

```r
renv::restore()
q()
```

## Repository Overview

### Baseline expression data
The scripts used to generate baseline expression datasets can be found in the `generate_baseline` directory. All baseline expression datasets are saved in AnnData format (`.h5ad`) and are used as the input reference for generating dropout versions in later steps of the workflow.
- `create_real_data.R` - create a real single-cell reference dataset (TS Ear, TS Skin, and Lung).
- `create_symsim_true_expression.R` - generate a SymSim-based baseline expression dataset.
- `create_splat_true_expression.R` - generate a Splatter-based baseline expression dataset.
- `create_sparsim_true_expression.R` - generate a SPARSim-based baseline expression dataset.

### Noise generation
The scripts used to add technical dropout can be found in the `generate_dropout` directory. All datasets are saved in AnnData format (`.h5ad`). 
- `dropout_alisa.R`
- `dropout_esco.R`
- `dropout_random.R`
- `dropout_sparsim_variable_library.R`
- `dropout_sparsim.R`
- `dropout_splat.R`
- `dropout_symsim.R`
- `dropout_zinbwave.R`

### Pipelines
Both single-cell analysis pipelines are implemented in `pipelines.py`. Additionally, the file contains the parent class of both pipelines "EmbeddingPipeline", a dataclass "PipelineResults" to store the experimental results and a wrapper function "timeit_store" to track the computing times for each step of the pipeline.

### Helper Files
1. `data_handling.py` - functions for loading datasets and finding input files.
2. `analysis.py` - utilities for parsing result filenames and collecting output files.
3. `plotting.py` - helper functions for formatting and generating summary plots.

### Job Scripts
- `scripts/submit_generate_baseline.sh`: script for running baseline files 
- `scripts/submit_generate_zeros.sh`: script for running dropout generation files 
- `scripts/submit_pipelines_gpu.sh`: runs experiment for chosen baseline expression and dropout generation method

### Notebooks
- `notebooks/AnalysisDropout.ipynb`: loads results, computes metrics, and generates plots.
- `notebooks/AnalysisRuntime.ipynb`: loads results from all the experiments, and plot runtime information

## Workflow
1. Simulate baseline expression data (`generate_baseline`).
2. Simulate observed count matrices using different noise models and sparsity levels (`scripts/submit_generate_zeros.sh`).
3. Copy the corresponding baseline expression file into each dropout-generated dataset folder.
4. Run embedding experiments using PCA and scVI (`scripts/submit_pipelines_gpu.sh`).
5. Analyze results in `notebooks/AnalysisDropout.ipynb` and `notebooks/AnalysisRuntime.ipynb`.

## Outputs
- .h5ad format (count matrices)
- .pkl format (pipeline results)
- plots, like UMAP visualizations
