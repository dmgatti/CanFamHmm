#!/bin/bash -l
cd /projects/dgatti/Dog/Scripts/CanFamHmm
for i in `seq 1 39`
do
  qsub -v CHR=$i,INPUT=/projects/dgatti/Dog/ReferenceData/PureBreedsHQ/PureBreedsHQ_,OUTPUT=/projects/dgatti/Dog/ReferenceData/PureBreedsHQ/shapeit/PureBreedsHQ_ \
       _5_phase_reference_samples_shapeit.sh
done
