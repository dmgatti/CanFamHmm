#!/bin/bash -l
cd /hpcdata/dgatti/Dog/Scripts/CanFamHmm
for i in `seq 1 39`
do
  qsub -v CHR=$i _5_phase_reference_samples_beagle.sh
done
