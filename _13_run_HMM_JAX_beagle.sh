#!/bin/bash -l
#PBS -q short -l nodes=1:ppn=1,walltime=2:00:00
module load R/3.3.2
cd /hpcdata/dgatti/Dog/Scripts/CanFamHmm
R --no-save --args ${INPUT_DIR} ${MODEL_DIR} ${OUTPUT_DIR} ${WINDOW} ${DOG} < _13_run_HMM_JAX.R > _13_run_HMM_JAX_${DOG}.Rout

