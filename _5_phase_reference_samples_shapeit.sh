#!/bin/bash -l
#PBS -l nodes=1:ppn=20,walltime=8:00:00

# Run SHAPEIT on one chromosome.
module load shapeit/2.r837

#########################################
# Arguments to pass in
# CHR: the chromosome to phase.
# INPUT: the full path and filename prefix to the unphased BED/BIM/FAM
#               files. We will append "chr" and the chromosome number to this.
# OUTPUT: the full path and filename prefix 
#########################################

# Number of threads.
NUM_BURN=7
NUM_PRUNE=8
NUM_MAIN=20
STATES=100
NE=10000
WINDOW=5

NTHREADS=10

INPUT_PREFIX=${INPUT}chr${CHR}
OUTPUT_PREFIX=${OUTPUT}chr${CHR}

shapeit --input-bed ${INPUT_PREFIX}.bed ${INPUT_PREFIX}.bim ${INPUT_PREFIX}.fam \
        --output-max ${OUTPUT_PREFIX}_phased.haps.gz ${OUTPUT_PREFIX}_phased.sample \
        --output-graph ${OUTPUT_PREFIX}.graph \
        --thread ${NTHREADS}

