#!/bin/bash

module load plink/1.9

cd /projects/dgatti/Dog/ReferenceData/PureBreedsHQ

for f in `ls *.bed`
do
  prefix=${f/.bed/}
  plink --bfile ${prefix} --dog --distance square ibs --out ${prefix}
  plink --bfile ${prefix} --dog --make-rel square --out ${prefix}
  plink --bfile ${prefix} --dog --homozyg group --out ${prefix}
done



