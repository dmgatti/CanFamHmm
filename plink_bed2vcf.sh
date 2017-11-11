#!/bin/bash -l
#PBS -l nodes=1:ppn=1,mem=16gb,walltime=8:00:00
# Convert the JAX samples from PLINK BED to VCF files.
module load plink/1.9

REFDIR=/hpcdata/dgatti/Dog/ReferenceData/PureBreedsHQ/
WORKDIR=/hpcdata/dgatti/Dog/GeneseekData/Cleaned/

cd ${WORKDIR}

for i in {1..39}
do

  BEDFILE=`ls *.bed | grep chr${i}.bed`
  PREFIX=${BEDFILE/.bed/}

  plink --bfile ${PREFIX} --dog --recode vcf --out ${PREFIX}

  gzip --force ${PREFIX}.vcf

done
