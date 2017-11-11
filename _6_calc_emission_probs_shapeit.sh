#!/bin/bash -l
#PBS -q short -l nodes=1:ppn=1,walltime=2:00:00
module load R/3.4.1

cd /projects/dgatti/Dog/Scripts/CanFamHmm

WIN=5
INDIR=/projects/dgatti/Dog/ReferenceData/PureBreedsHQ/shapeit/
OUTDIR=/projects/dgatti/Dog/InputData/PureBreedsAllHQ/shapeit/win${WIN}/
MAPDIR=/projects/dgatti/Dog/ReferenceData/PureBreedsHQ/
OBSDIR=/projects/dgatti/Dog/GeneseekData/Cleaned/shapeit/

R --no-save --args ${INDIR} ${OUTDIR} ${MAPDIR} ${OBSDIR} ${WIN} < _6_calc_emission_probs_shapeit.R > _6_calc_emission_probs_shapeit.Rout

# Copy the map files and rename them.
cd ${MAPDIR}
for i in $(ls *.map)
do
  cp $i ${OUTDIR}${i/PureBreedsHQ_/}
done
