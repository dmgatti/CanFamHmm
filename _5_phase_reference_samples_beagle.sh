#!/bin/bash -l
#PBS -l nodes=1:ppn=20,walltime=8:00:00

# Run Beagle on one chromosome.
module load java/1.8.0_73
module load samtools/1.3.1

#########################################
# NOTE: CHR is set by the calling script.
#CHR=1
#########################################

BEAGLE_JAR=/home/dgatti/software/beagle/beagle.21Jan17.6cc.jar
WD=/hpcdata/dgatti/Dog/ReferenceData/PureBreedsHQ/
VCF_FILE=PureBreedsHQ_chr${CHR}.vcf.gz
MAP_FILE=PureBreedsHQ_chr${CHR}.map
OUT_FILE=PureBreedsHQ_chr${CHR}.phased
# Beagle manual recommends ~ 5cM.
WINDOW=500
# Beagle manual recommends ~ 0.5cM.
OVERLAP=50
# Number of phasing iterations.
NITER=10
# Number of threads.
NTHREADS=20
# Effective population size.
NE=10000

cd ${WD}

# Run Beagle
java -Xss5m -Xmx200g -jar ${BEAGLE_JAR} gt=$VCF_FILE out=$OUT_FILE \
  map=$MAP_FILE nthreads=$NTHREADS ne=$NE window=$WINDOW overlap=$OVERLAP \
  niterations=$NITER ibd=true

# Index the VCF file.
tabix -f -p vcf ${OUT_FILE}.vcf.gz

