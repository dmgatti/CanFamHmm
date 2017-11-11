#!/bin/bash -l
#PBS -l nodes=1:ppn=1,mem=8GB,walltime=1:00:00
module load R/3.3.2
cd /projects/dgatti/Dog/Scripts/CanFamHmm
R --no-save --args ${NUMSAMPLES} ${CHUNK} ${PHASED} ${REFDIR} ${MAPDIR} ${OUTDIR} < simulate_genomes_engine.R \
   > simulate_genomes_engine_${CHUNK}.Rout

