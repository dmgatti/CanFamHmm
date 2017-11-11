cd /projects/dgatti/Dog/Scripts/CanFamHmm

rm *.Rout
rm *.sh.e*
rm *.sh.o*

# SNP window size.
WIN=10
# Input directory.
IN_DIR=/projects/dgatti/Dog/ReferenceData/PureBreedsHQ/shapeit/win${WIN}/
# Model directory.
MOD_DIR=/projects/dgatti/Dog/InputData/PureBreedsAllHQ/shapeit/win${WIN}/
# Output directory.
OUT_DIR=/projects/dgatti/Dog/HMM/HMM_PureBreedHQ/win${WIN}/
# Transition probability factor.
TRANSFACT=-91

for I in `seq 1 10 10000`
do
  qsub -v INPUT_DIR=${IN_DIR},MODEL_DIR=${MOD_DIR},OUTPUT_DIR=${OUT_DIR},WINDOW=${WIN},TPF=${TRANSFACT},DOG=${I} HMM_ref_samples.sh
  sleep 1s
done
