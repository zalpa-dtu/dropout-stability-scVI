#!/bin/bash

# Usage:
# bash scripts/submit_generate_zeros.sh \
#   SCRIPT DATAPATH OUTPUTPATH CONVERT_TO_COUNTS \
#   MAX_LIB_SIZE DATANAME DATATYPE PARAMETERS

SCRIPT="$1"
DATAPATH="$2"
OUTPUTPATH="$3"
CONVERT_TO_COUNTS="$4"
MAX_LIB_SIZE="$5"
DATANAME="$6"
DATATYPE="$7"
PARAMETERS="$8"

JOB_DIR="${OUTPUTPATH%/}/jobs"
DATA_DIR="${OUTPUTPATH%/}/Data"

mkdir -p "${JOB_DIR}"
mkdir -p "${DATA_DIR}"

# Override the default project directory with:
# PROJECT_DIR="/path/to/project" bash scripts/submit_generate_zeros.sh ...
PROJECT_DIR="${PROJECT_DIR:-/dtu-compute/digitstem/zuzanna/zuzanna_dropout_project}"
PYTHON_ENV="${PROJECT_DIR}/.conda-env"

# R module used on the DTU cluster. Override with:
# R_MODULE="R/4.5.2" bash scripts/submit_generate_zeros.sh ...
R_MODULE="${R_MODULE:-R/4.5.2-mkl2025update2}"

bsub -q hpc \
     -o ${JOB_DIR}/Output_${DATANAME}_%J.out \
     -e ${JOB_DIR}/Error_${DATANAME}_%J.err \
     -J "GenerateZeros_${DATANAME}" \
     -n 1 \
     -R "span[hosts=1]" \
     -R "rusage[mem=16GB]" \
     -W 24:00 \
     "module load ${R_MODULE} && \
      export RETICULATE_PYTHON=${PYTHON_ENV}/bin/python && \
      Rscript generate_dropout/${SCRIPT}.R ${DATAPATH} ${DATA_DIR}/ ${CONVERT_TO_COUNTS} ${MAX_LIB_SIZE} ${DATANAME} ${DATATYPE} ${PARAMETERS}"