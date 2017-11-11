#!/bin/bash -l
cd /hpcdata/dgatti/Dog/Scripts/CanFamHmm
for i in `seq 1 39`
do
  qsub -v CHR=$i _11_phase_JAX_beagle.sh
done
