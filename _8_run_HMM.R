########################################################################
# Haplotype window HMM.
# Run the HMM on one dog on one chromosome.
# Daniel Gatti
# dan.gatti@jax.org
# May, 12, 2017
########################################################################
options(stringsAsFactors = F)
library(Rcpp)

src.dir    = "/data/dgatti/Dog/Scripts/CanFamHmm/"
sourceCpp(paste0(src.dir, "HMM.cpp"))

#obs.dir  = "/data/dgatti/Dog/SimulatedData/PureBreedsHQ/"
#emis.dir  = "/data/dgatti/Dog/InputData/PureBreedsAllHQ/"
#output.dir = "/data/dgatti/Dog/ReferenceData/PureBreedsHQ/"
#window = 1
#dog = 1

# Get command line arguments.
# Arg 1 is input file directory.
# Arg 2 is emission file directory.
# Arg 3 is output file directory.
# Arg 4 is window size.
# Arg 5 is dog index.
args = commandArgs(trailingOnly = TRUE)
obs.dir  = args[1]
emis.dir   = args[2]
output.dir = args[3]
window = as.numeric(args[4])
dog = as.numeric(args[5])

probs = as.list(1:39)
names(probs) = 1:39

vit = as.list(1:39)
names(vit) = 1:39

for(chr in 1:39) {

  print(paste("CHR", chr))

  # Read in observations.
# NEED TO STANDARDIZE INPUT FILE NAMES.
  obs = readRDS(paste0(obs.dir, "simulated_chr", chr, "_window",
                window, ".rds"))

  # Read in emission probs.
# NEED TO STANDARDIZE EMISSION FILE NAMES.
  emis = readRDS(paste0(emis.dir, "emission_probs_chr", chr, "_window", 
                 window, ".rds"))
  emis = log(emis)
  breeds = unique(sub("_[0-9]+$", "", rownames(emis)))

  # Read in marker map.
# NEED TO STANDARDIZE MAP FILE NAMES.
  map = read.delim(paste0(emis.dir, "chr", chr, ".map"),
        sep = "", header = F)

  probs[[chr]] = forwardBackward(obs = as.integer(obs[dog,]), emis = emis, 
                 map = map[,3] * 1e-8, window = window)
  probs[[chr]] = exp(probs[[chr]])
  dimnames(probs[[chr]]) = list(breeds, colnames(emis))

  # Add one to the Viterbi result because C++ is 0-based and R is 1-based.
  vit[[chr]] = viterbi(as.integer(obs[dog,]), emis = emis,
               map = map[,3] * 1e-8, window = window) + 1

} # for(chr)

saveRDS(probs, file = paste0(output.dir, colnames(obs)[dog], "_HMM_probs.rds"))
saveRDS(vit, file = paste0(output.dir, colnames(obs)[dog], "_HMM_viterbi.rds"))


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

