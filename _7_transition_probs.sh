#!/bin/bash -l
#PBS -l nodes=1:ppn=20,walltime=12:00:00
module load R/3.3.2
cd /hpcdata/dgatti/Dog/Scripts/CanFamHmm
R CMD BATCH --no-save _7_transition_probs.R
