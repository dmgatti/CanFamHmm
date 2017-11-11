# Simulated dogs HMM.
cd /projects/dgatti/Dog/Scripts/CanFamHmm

rm *.Rout
rm *.sh.e*
rm *.sh.o*

# SNP window size.
WIN=5
# Input directory.
IN_DIR=/hpcdata/dgatti/Dog/SimulatedData/PureBreedHQ/phased/shapeit/multi/win${WIN}/
# Model directory.
MOD_DIR=/projects/dgatti/Dog/InputData/PureBreedsAllHQ/shapeit/win${WIN}/
# Output directory.
OUT_DIR=/hpcdata/dgatti/Dog/SimulatedData/PureBreedHQ/HMM/shapeit/win${WIN}/
# Transition probability factor.
TRANSFACT=-91


for I in `seq 1 2000`
do
  qsub  -v INPUT_DIR=${IN_DIR},MODEL_DIR=${MOD_DIR},OUTPUT_DIR=${OUT_DIR},WINDOW=${WIN},TPF=${TRANSFACT},DOG=${I} run_HMM_simulated_genomes.sh
  sleep 1s
done
