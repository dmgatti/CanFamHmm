#!/bin/bash -l
#PBS -q short -l nodes=1:ppn=1,mem=16gb,walltime=3:00:00
module load R/3.4.1
cd /projects/dgatti/Dog/Scripts/CanFamHmm
R --no-save --args ${INPUT_DIR} ${MODEL_DIR} ${OUTPUT_DIR} ${WINDOW} ${TPF} ${DOG} < HMM_ref_samples.R > HMM_ref_samples_${DOG}.Rout

