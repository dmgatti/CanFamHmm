#!/bin/bash -l
#PBS -q short -l nodes=1:ppn=1,walltime=2:00:00
module load R/3.3.2
cd /hpcdata/dgatti/Dog/Scripts/CanFamHmm
R --no-save --args ${INPUT_DIR} ${OUTPUT_DIR} ${WINDOW} ${DOG} < _8_run_HMM.R > _8_run_HMM_${DOG}.Rout

