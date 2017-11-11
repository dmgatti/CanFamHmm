#!/bin/bash -l
#PBS -l nodes=1:ppn=1,mem=128gb,walltime=24:00:00
module load R/3.4.1
cd /projects/dgatti/Dog/Scripts/CanFamHmm
WIN=10
R --no-save --args ${WIN} < transition_probs_adjust.R > transition_probs_adjust_${WIN}.Rout

