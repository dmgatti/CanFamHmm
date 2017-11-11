# Check that the emission probs and obs with the given window 
# are correct.
options(stringsAsFactors = F)

window = 10

ref.dir = "/projects/dgatti/Dog/ReferenceData/PureBreedsHQ/shapeit/"
hmm.dir = "/projects/dgatti/Dog/InputData/PureBreedsAllHQ/shapeit/win10/"
obs.dir = "/projects/dgatti/Dog/GeneseekData/Cleaned/shapeit/"
obs.win.dir = "/projects/dgatti/Dog/GeneseekData/Cleaned/shapeit/win10/"

chr = 1

# Read in the ref data.
ref  = readRDS(paste0(ref.dir, "PureBreedsHQ_chr1_phased.rds"))
emis = readRDS(paste0(hmm.dir, "emission_probs_chr1.rds"))

# Read in the obs data.
obs = readRDS(paste0(obs.dir, "jax_chr1_phased.rds"))
obs.win = readRDS(paste0(obs.win.dir, "jax_chr1_phased.rds"))

ref = ref[,colnames(obs)]

# Calculate the first window values.
vec = 2^(0:(window-1))

ref.hap = matrix(0, nrow(ref), ncol(ref) - window)
for(i in 1:ncol(ref.hap)) {
  ref.hap[,i] = ref[,i:(window + i - 1)] %*% vec
}

hap = matrix(0, nrow(obs), ncol(obs) - window)
for(i in 1:ncol(hap)) {
  hap[,i] = obs[,i:(window + i - 1)] %*% vec
}

un.ref.hap = apply(ref.hap, 2, unique)
un.hap = apply(hap, 2, unique)

m = mapply(function(x,y) { mean(x %in% y) }, x = un.hap, y = un.ref.hap)


