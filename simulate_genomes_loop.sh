#!/bin/bash -l
cd /projects/dgatti/Dog/Scripts/CanFamHmm

N=1000

# Reference sample data.
RD=/projects/dgatti/Dog/ReferenceData/PureBreedsHQ/shapeit/
# Map file directory.
MD=/projects/dgatti/Dog/ReferenceData/PureBreedsHQ/
# Output file directory.
OD_BASE=/projects/dgatti/Dog/SimulatedData/PureBreedHQ/phased/shapeit/

for CHUNK in 1 2 4 8 16
do
  OD=/projects/dgatti/Dog/SimulatedData/PureBreedHQ/phased/shapeit/Mb_${CHUNK}/
  qsub -v NUMSAMPLES=${N},CHUNK=${CHUNK},PHASED=TRUE,REFDIR=${RD},MAPDIR=${MD},OUTDIR=${OD} simulate_genomes_engine.sh
  sleep 1s
done
