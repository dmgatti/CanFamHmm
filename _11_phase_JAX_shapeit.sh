#!/bin/bash -l
#PBS -l nodes=1:ppn=1,walltime=1:00:00

# Run SHAPEIT on one chromosome.
module load R/3.3.2
module load shapeit/2.r837

#########################################
# Variables to be set by the calling script.
# CHR=1
# REFINDIR: Path to the input directory for reference files.
# REFOUTDIR: Path to the output directory for reference files.
# MAPDIR: Path to the map files.
# DATADIR: Path to the JAX data files.
# OUTDIR: Output path to the JAX phased data files.
#########################################

SRC_DIR=/projects/dgatti/Dog/Scripts/CanFamHmm/

# Convert the SHAPEIT refrence files from HAP/SAMPLE format to reference format.
R --no-save --args ${CHR} ${REFINDIR} ${REFOUTDIR} < ${SRC_DIR}convert_shapeit2refhaps.R

# Write out a list of SNPs to exclude with 100% no-calls in the JAX data.
# SHAPEIT throws an error in this case.
R --no-save --args ${CHR} ${DATADIR} < ${SRC_DIR}find_missing_snps.R


# Prefix for phased reference files in HAP/LEGEND/SAMPLE format, created by 
# convert_shapeit2refhaps.R
REF_PREFIX=${REFOUTDIR}ref_chr${CHR}

DATA_PREFIX=${DATADIR}jax_chr${CHR}

MAP_FILE=${MAPDIR}jax_chr${CHR}.map

OUT_PREFIX=${OUTDIR}jax_chr${CHR}_phased

# Number of burn in iterations (--burn).
NBURN=7
# Number of pruning in iterations (--prune).
NPRUNE=8
# Number of main in iterations (--main).
NMAIN=20
# Number of allowed states at each marker (--states).
NSTATES=100
# Window size in Mb (--window).
NSTATES=100
# Effective population size (--effective-size).
NE=15000
# Nubmer of parallel threads.
NTHREADS=1


# Run SHAPEIT
shapeit --input-bed ${DATA_PREFIX}.bed ${DATA_PREFIX}.bim ${DATA_PREFIX}.fam \
        --input-ref ${REF_PREFIX}.haps.gz ${REF_PREFIX}.legend ${REF_PREFIX}.sample \
        --output-max ${OUT_PREFIX} \
        --exclude-snp ${DATA_PREFIX}.exclude \
        --thread ${NTHREADS}

# Convert the phased JAX haplotype files to *.rds files for input into the HMM.
R --no-save --args ${CHR} ${OUTDIR} ${MAPDIR} < ${SRC_DIR}convert_shapeit2rds.R

#shapeit -check \
#        -B ${DATA_PREFIX} \
#        --input-ref ${REF_PREFIX}.haps.gz ${REF_PREFIX}.legend ${REF_PREFIX}.sample \
#        --exclude-snp ${DATA_PREFIX}.exclude \
#        --output-log gwas.alignments


