#!/bin/bash

# Usage:
#   bash scripts/submit_generate_baseline.sh SCRIPT OUTPUTPATH [SCRIPT_ARGUMENTS...]

SCRIPT="$1"
OUTPUTPATH="$2"
shift 2

JOB_DIR="${OUTPUTPATH%/}/jobs"
DATA_DIR="${OUTPUTPATH%/}/"

mkdir -p "${JOB_DIR}"
mkdir -p "${DATA_DIR}"

# Override the default project directory with:
# PROJECT_DIR="/path/to/project" bash scripts/submit_generate_zeros.sh ...
PROJECT_DIR="${PROJECT_DIR:-/dtu-compute/digitstem/zuzanna/zuzanna_dropout_project}"
PYTHON_ENV="${PROJECT_DIR}/.conda-env"

# R module used on the DTU cluster. Override with:
# R_MODULE="R/4.5.2" bash scripts/submit_generate_zeros.sh ...
R_MODULE="${R_MODULE:-R/4.5.2-mkl2025update2}"

JOB_NAME="GenerateBaseline_${SCRIPT}"

bsub -q hpc \
     -o ${JOB_DIR}/Output_%J.out \
     -e ${JOB_DIR}/Error_%J.err \
     -J "${JOB_NAME}" \
     -n 1 \
     -R "span[hosts=1]" \
     -R "rusage[mem=16GB]" \
     -W 24:00 \
     "module load ${R_MODULE} && \
      export RETICULATE_PYTHON=${PYTHON_ENV}/bin/python && \
      Rscript generate_baseline/${SCRIPT}.R $* ${DATA_DIR}"