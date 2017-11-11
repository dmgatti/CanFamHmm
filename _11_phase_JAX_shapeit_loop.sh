#!/bin/bash -l
cd /projects/dgatti/Dog/Scripts/CanFamHmm

REFIN=/projects/dgatti/Dog/ReferenceData/PureBreedsHQ/shapeit/
REFOUT=/projects/dgatti/Dog/ReferenceData/PureBreedsHQ/shapeit/ref/
MAP=/projects/dgatti/Dog/ReferenceData/PureBreedsHQ/
DATA=/hpcdata/dgatti/Dog/GeneseekData/Cleaned/
OUT=/projects/dgatti/Dog/GeneseekData/Cleaned/shapeit/

for i in `seq 1 39`
do
  qsub -v CHR=$i,REFINDIR=${REFIN},REFOUTDIR=${REFOUT},MAPDIR=${MAP},DATADIR=${DATA},OUTDIR=${OUT} _11_phase_JAX_shapeit.sh
  sleep 1s
done
