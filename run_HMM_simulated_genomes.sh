#!/bin/bash -l
#PBS -q short -l nodes=1:ppn=1,mem=8gb,walltime=3:59:00
module load R/3.4.1
cd /projects/dgatti/Dog/Scripts/CanFamHmm

R --no-save --args ${INPUT_DIR} ${MODEL_DIR} ${OUTPUT_DIR} ${WINDOW} ${TPF} ${DOG} < run_HMM_simulated_genomes.R > run_HMM_simulated_genomes_${DOG}.Rout

