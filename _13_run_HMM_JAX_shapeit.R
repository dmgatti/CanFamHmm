########################################################################
# Haplotype window HMM.
# Run the HMM on one dog on one chromosome.
# Parameters:
# win = 1, trans.factor = -29
# win = 3, trans.factor = -24
# win = 5, trans.factor = -21
# win = 7, trans.factor = -17
# win = 10, trans.factor = -12
# Daniel Gatti
# dan.gatti@jax.org
# May, 12, 2017
########################################################################
options(stringsAsFactors = F)
library(Rcpp)

src.dir    = "/projects/dgatti/Dog/Scripts/CanFamHmm/"
sourceCpp(paste0(src.dir, "HMM.cpp"))

window = 10
obs.dir  = paste0("/projects/dgatti/Dog/GeneseekData/Cleaned/shapeit/win",
           window, "/")
emis.dir  = paste0("/projects/dgatti/Dog/InputData/PureBreedsAllHQ/shapeit/win",
            window, "/")
output.dir = paste0("/projects/dgatti/Dog/HMM/HMM_v7_SHAPEIT/win",
             window, "/")
dog = 1

# Get command line arguments.
# Arg 1 is input file directory.
# Arg 2 is emission file directory.
# Arg 3 is output file directory.
# Arg 4 is window size.
# Arg 5 is transition probability factor.
# Arg 6 is dog index.
args = commandArgs(trailingOnly = TRUE)

if(length(args) != 6) {
  stop(paste("6 arguments required. obs.dir, emission.dir, output.dir,",
      "window, trans prob factor & dog index."))
}

obs.dir  = args[1]
emis.dir   = args[2]
output.dir = args[3]
window = as.numeric(args[4])
trans.factor = as.numeric(args[5])
dog = as.numeric(args[6])

# Set the transtion probs multiplier. This is a negative number that will
# be subtracted from the cM values on the log scale.

probs = as.list(1:39)
names(probs) = 1:39

vit = as.list(1:39)
names(vit) = 1:39

for(chr in 1:39) {

  print(paste("CHR", chr))

  # Read in observations.
# NEED TO STANDARDIZE INPUT FILE NAMES.
  obs = readRDS(paste0(obs.dir, "jax_chr", chr, "_phased.rds"))

  # Read in emission probs.
# NEED TO STANDARDIZE EMISSION FILE NAMES.
  emis = readRDS(paste0(emis.dir, "emission_probs_chr", chr, ".rds"))
  emis = log(emis)
  breeds = unique(sub("_[0-9]+$", "", rownames(emis)))

  # Read in marker map.
# NEED TO STANDARDIZE MAP FILE NAMES.
  map = read.delim(paste0(obs.dir, "jax_chr", chr, ".map"),
        sep = "", header = F)

  # NOTE: for SHAPEIT, we have to subset the markers.
  mkr = intersect(colnames(obs), colnames(emis))
  obs  = obs[,mkr]
  emis = emis[,mkr]
  map  = map[map[,2] %in% mkr,]

  stopifnot(colnames(obs) == colnames(emis))

#  probs[[chr]] = forwardBackward(obs = as.integer(obs[dog,]), emis = emis, 
#                 map = map[,3], mapfactor = trans.factor, window = window)
#  probs[[chr]] = exp(probs[[chr]])
#  dimnames(probs[[chr]]) = list(breeds, colnames(emis))

  # Add one to the Viterbi result because C++ is 0-based and R is 1-based.
  vit[[chr]] = viterbi(obs = as.integer(obs[dog,]), emis = emis,
               map = map[,3], mapfactor = trans.factor, window = window)
  vit[[chr]][,1] = vit[[chr]][,1] + 1
  dimnames(vit[[chr]]) = list(map[,2], c("breed", "prob"))

} # for(chr)

saveRDS(probs, file = paste0(output.dir, rownames(obs)[dog], "_HMM_probs.rds"))
saveRDS(vit, file = paste0(output.dir, rownames(obs)[dog], "_HMM_viterbi.rds"))


#probs = forwardBackward(obs = as.integer(jax[,dog]), emis = emis, 
#                 map = map[,3] * 1e-6)
#probs = exp(probs)
#dimnames(probs) = list(breeds, colnames(emis))
#plot(apply(probs, 2, which.max))
#sort(rowMeans(probs))

#probs[[chr]] = probs[[chr]][nrow(probs[[chr]]):1,]

#par(plt = c(0.1, 0.9, 0.06, 0.95))
#image(1:ncol(probs[[chr]]), 1:nrow(probs[[chr]]), t(probs[[chr]]), 
#      yaxt = "n", ann = F)
#mtext(side = 2, line = 0.2, at = 1:nrow(probs[[chr]]), text = 
#     rownames(probs[[chr]]), las = 2, cex = 0.5)
#mtext(side = 4, line = 0.2, at = 1:nrow(probs[[chr]]), text = 
#     rownames(probs[[chr]]), las = 2, cex = 0.5)
#abline(h = 1:nrow(probs[[chr]]) + 0.5, col = "grey80")

