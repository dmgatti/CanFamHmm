#!/bin/bash -l
#PBS -q short -l nodes=1:ppn=1,mem=24gb,walltime=3:59:00
module load R/3.4.1

cd /projects/dgatti/Dog/Scripts/CanFamHmm

WIN=1
SIM_DIR=/hpcdata/dgatti/Dog/SimulatedData/PureBreedHQ/phased/shapeit/multi/
HMM_DIR=/hpcdata/dgatti/Dog/SimulatedData/PureBreedHQ/HMM/shapeit/win${WIN}/
MAP_DIR=/projects/dgatti/Dog/InputData/PureBreedsAllHQ/shapeit/win${WIN}/
OUT_DIR=/hpcdata/dgatti/Dog/SimulatedData/PureBreedHQ/HMM/shapeit/
TPF=-36

R --no-save --args ${SIM_DIR} ${HMM_DIR} ${MAP_DIR} ${OUT_DIR} ${WIN} ${TPF} < evaluate_simulations.R > evaluate_simulations_${WIN}.Rout

