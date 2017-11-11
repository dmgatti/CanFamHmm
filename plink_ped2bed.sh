#!/bin/bash -l
# Convert the JAX samples from PLINK PED to BED files and zip up the PED files.
# Use alleles from the reference data.
module load plink/1.9

REFDIR=/hpcdata/dgatti/Dog/ReferenceData/PureBreedsHQ/
WORKDIR=/hpcdata/dgatti/Dog/GeneseekData/Cleaned/

cd ${WORKDIR}

for i in {1..39}
do

  PEDFILE=`ls *.ped | grep chr${i}.ped`
  PREFIX=${PEDFILE/.ped/}
  BIMFILE=`ls ${REFDIR}*.bim | grep chr${i}.bim`

  plink --file ${PREFIX} --dog --make-bed --a2-allele ${BIMFILE} 6 2 --out ${PREFIX}
  gzip --force ${PREFIX}.ped

done
