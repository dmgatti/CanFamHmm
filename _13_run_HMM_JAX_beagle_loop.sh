cd /data/dgatti/Dog/Scripts/CanFamHmm
# Input directory.
IN_DIR=/projects/dgatti/Dog/GeneseekData/Cleaned/
# Model directory.
MOD_DIR=/hpcdata/dgatti/Dog/InputData/PureBreedsAllHQ/
# Output directory.
OUT_DIR=/projects/dgatti/Dog/HMM/HMM_v6/win1/
# Set SNP window size.
WIN=1

for I in `seq 1 190`
do
  qsub -v INPUT_DIR=${IN_DIR},MODEL_DIR=${MOD_DIR},OUTPUT_DIR=${OUT_DIR},WINDOW=${WIN},DOG=${I} _13_run_HMM_JAX.sh
done
