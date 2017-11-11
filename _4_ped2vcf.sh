#!#/bin/bash
# Convert PLINK PED input files into VCF format and gzip.
# Also write out the PLINK files in binary format and gzip the PED files.
module load plink/1.9

INPUT_DIR=/hpcdata/dgatti/Dog/ReferenceData/PureBreedsHQ
OUTPUT_DIR=/hpcdata/dgatti/Dog/ReferenceData/PureBreedsHQ
FILE_PREFIX=PureBreedsHQ_chr

for CHR in `seq 1 39`
do

  plink --file ${INPUT_DIR}/${FILE_PREFIX}${CHR} --dog --recode vcf --out ${OUTPUT_DIR}/${FILE_PREFIX}${CHR}
  gzip --force ${OUTPUT_DIR}/${FILE_PREFIX}${CHR}.vcf

  plink --file ${INPUT_DIR}/${FILE_PREFIX}${CHR} --dog --make-bed --out ${OUTPUT_DIR}/${FILE_PREFIX}${CHR}
  gzip --force ${OUTPUT_DIR}/${FILE_PREFIX}${CHR}.ped

done


