#!/bin/bash -l
# Simulate genomes on one chromosome, run the HMM and evaluate the quality of 
# the reconstructions.

# Load in required software.
module load R/3.3.2
module load beagle-lib/2.1.2

# Input variables.
# Observation directory. Simulated observations will go here.
OBS_DIR=/projects/dgatti/Dog/SimulatedData/PureBreedHQ/phased/win1/Mb_1/
# Model directory. Contains emission probs, transition probs and map.
MODEL_DIR=/data/dgatti/Dog/InputData/PureBreedsAllHQ/
# Phased reference directory. Contains phased VCF reference files.
REF_DIR=/data/dgatti/Dog/ReferenceData/PureBreedsHQ/
# Output directory.
OUT_DIR=/data/dgatti/Dog/HMM/sim_phased_window1/
# Script directory.
SCRIPT_DIR=/data/dgatti/Dog/Scripts/CanFamHmm/
# Number of genomes to simulate.
N=100
# Chromosome.
CHR=1
# Number of SNPs in window.
WINDOW=1
# TRUE if the allele calls are swapped to disrupt the phasing, FALSE otherwise.
SCRAMBLED=FALSE

# Simulate genomes.
R --no-save --args ${CHR} ${N} ${SCRAMBLED} ${REF_DIR} ${OBS_DIR} < \
  ${SCRIPT_DIR}simulate_genomes_engine.R > \
  ${SCRIPT_DIR}simulate_genomes_engine_${CHR}.Rout

# Phase the simulated genomes.
# TBD use Beagle with phased reference data.

# Reconsruct genomes using HMM.
for DOG in `seq 1 ${N}`
do
  R --no-save --args ${OBS_DIR} ${MODEL_DIR} ${OUT_DIR} ${WINDOW} ${DOG} < \
    ${SCRIPT_DIR}_8_run_HMM.R > \
    ${SCRIPT_DIR}_8_run_HMM_${DOG}.Rout
done


# Evaluate the proportion of each genome that was reconstructed correctly.
R --no-save --args ${OBS_DIR} ${OUT_DIR} ${WINDOW} < evaluate_simulations.R \
  > evaluate_simulations.Rout

