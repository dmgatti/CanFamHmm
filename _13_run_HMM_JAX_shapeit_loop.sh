# Parameters:
# win = 1, trans.factor = -36
# win = 3, trans.factor = -24
# win = 5, trans.factor = -91
# win = 7, trans.factor = -17
# win = 10, trans.factor = -200

cd /projects/dgatti/Dog/Scripts/CanFamHmm

# SNP window size.
WIN=10
# Input directory.
IN_DIR=/projects/dgatti/Dog/GeneseekData/Cleaned/shapeit/win${WIN}/
# Model directory.
MOD_DIR=/projects/dgatti/Dog/InputData/PureBreedsAllHQ/shapeit/win${WIN}/
# Output directory.
OUT_DIR=/projects/dgatti/Dog/HMM/HMM_v7_SHAPEIT/win${WIN}/
# Transition probability factor.
TRANSFACT=-200

for I in `seq 1 214`
do
  qsub -v INPUT_DIR=${IN_DIR},MODEL_DIR=${MOD_DIR},OUTPUT_DIR=${OUT_DIR},WINDOW=${WIN},TPF=${TRANSFACT},DOG=${I} _13_run_HMM_JAX_shapeit.sh
  sleep 1
done
