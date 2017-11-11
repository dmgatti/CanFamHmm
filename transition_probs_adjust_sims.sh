#!/bin/bash -l
#PBS -q batch -l nodes=1:ppn=1,walltime=96:00:00
module load R/3.4.1
cd /projects/dgatti/Dog/Scripts/CanFamHmm
WIN=1
R --no-save --args ${WIN} < transition_probs_adjust_sims.R > transition_probs_adjust_sims_${WIN}.Rout

