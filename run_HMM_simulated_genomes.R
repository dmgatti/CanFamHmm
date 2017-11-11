###############################################################################
# Haplotype window HMM for simulated samples.
# Run the HMM on one dog.
# Daniel Gatti
# dan.gatti@jax.org
# June 14, 2017
###############################################################################
options(stringsAsFactors = F)
library(Rcpp)

src.dir    = "/projects/dgatti/Dog/Scripts/CanFamHmm/"
sourceCpp(paste0(src.dir, "HMM.cpp"))

# Get command line arguments.
# Arg 1 is input file directory.
# Arg 2 is emission file directory.
# Arg 3 is output file directory.
# Arg 4 is window size.
# Arg 5 is transition probability factor.
# Arg 6 is dog index.
args = commandArgs(trailingOnly = TRUE)

print(paste("length(args) =", length(args)))

if(length(args) != 6) {
  stop(paste("6 arguments required. obs.dir, emission.dir, output.dir,",
      "window, trans prob factor & dog index."))
}

obs.dir    = args[1]
emis.dir   = args[2]
output.dir = args[3]
window = as.numeric(args[4])
trans.factor = as.numeric(args[5])
dog = as.numeric(args[6])

# window = 10
# obs.dir = paste0("/hpcdata/dgatti/Dog/SimulatedData/PureBreedHQ/phased/shapeit/multi/win", window ,"/")
# emis.dir = paste0("/projects/dgatti/Dog/InputData/PureBreedsAllHQ/shapeit/win", window ,"/")
# output.dir = paste0("/hpcdata/dgatti/Dog/SimulatedData/PureBreedHQ/HMM/shapeit/win", window, "/")
# trans.factor = -32
# dog = 1

print(paste("DOG:", dog))

#probs = as.list(1:39)
#names(probs) = 1:39

vit = as.list(1:39)
names(vit) = 1:39

for(chr in 1:39) {

  print(paste("CHR", chr))

  # Read in simulated observations.
  obs = readRDS(paste0(obs.dir, "simulated_chr", chr, "_phased.rds"))
  rownames(obs) = make.unique(rownames(obs))

  # Read in emission probs.
  emis = readRDS(paste0(emis.dir, "emission_probs_chr", chr, ".rds"))
  emis = log(emis)
  breeds = unique(sub("_[0-9]+$", "", rownames(emis)))

  # Read in marker map.
  map = read.delim(paste0(emis.dir, "chr", chr, ".map"),
                   sep = "", header = F)

  # NOTE: for SHAPEIT, we have to subset the markers.
  common.markers = intersect(colnames(obs), colnames(emis))
  obs  = obs[,common.markers]
  emis = emis[,common.markers]
  map  = map[map[,2] %in% common.markers,]

  stopifnot(colnames(obs) == colnames(emis))

#  probs[[chr]] = forwardBackward(obs = as.integer(obs[dog,]), emis = emis, 
#                 map = map[,3], mapfactor = trans.factor, window = window)
#  probs[[chr]] = exp(probs[[chr]])
#  dimnames(probs[[chr]]) = list(breeds, colnames(emis))

  # Add one to the Viterbi result because C++ is 0-based and R is 1-based.
  vit[[chr]] = viterbi(obs = as.integer(obs[dog,]), emis = emis,
               map = map[,3], mapfactor = trans.factor, window = window)
  vit[[chr]][,1] = vit[[chr]][,1] + 1
  dimnames(vit[[chr]]) = list(colnames(emis), c("breed", "prob"))

  # Reverse Viterbi.
#  tmp = viterbi(obs = rev(as.integer(obs[dog,])), emis = emis[,ncol(emis):1],
#                map = map[,3], mapfactor = trans.factor, window = window)
#  tmp[,1] = tmp[,1] + 1

} # for(chr)

#saveRDS(probs, file = paste0(output.dir, rownames(obs)[dog],
#        "_HMM_probs_win", window, ".rds"))
saveRDS(vit, file = paste0(output.dir, rownames(obs)[dog],
        "_viterbi_win", window, ".rds"))


#par(plt = c(0.15, 0.99, 0.03, 0.97))
#image(1:ncol(probs[[chr]]), 1:nrow(probs[[chr]]), t(probs[[chr]]), yaxt = "n", ann = F)
#mtext(side = 2, at = 1:nrow(probs[[chr]]), line = 0.2, las = 2, text = rownames(probs[[chr]]))

