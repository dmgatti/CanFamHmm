#!/bin/bash -l
#PBS -l nodes=1:ppn=20,walltime=12:00:00
cd /hpcdata/dgatti/Dog/Scripts/CanFamHmm
module load R/3.4.1

R CMD BATCH --no-save _3_merge_reference.R

