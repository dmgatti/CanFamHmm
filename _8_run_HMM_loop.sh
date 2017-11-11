cd /data/dgatti/Dog/Scripts/CanFamHmm
# Input directory.
IN_DIR=/hpcdata/dgatti/Dog/InputData/PureBreedsAllHQ/
# Model directory.
MOD_DIR=/hpcdata/dgatti/Dog/InputData/PureBreedsAllHQ/
# Output directory.
OUT_DIR=/hpcdata/dgatti/Dog/HMM/HMM_PureBreedHQ_win10/
# Set SNP window size.
WIN=10

for I in `seq 1 190`
do
  qsub -v INPUT_DIR=${IN_DIR},MODEL_DIR=${MOD_DIR},OUTPUT_DIR=${OUT_DIR},WINDOW=${WIN},DOG=${I} _8_run_HMM.sh
done
