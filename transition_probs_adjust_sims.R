# Compare different transition prob. factors and compare to the true simulations.
options(stringsAsFactors = F)
library(Rcpp)
library(BiocParallel)

src.dir    = "/projects/dgatti/Dog/Scripts/CanFamHmm/"
sourceCpp(paste0(src.dir, "HMM.cpp"))

#window = 1
args = commandArgs(trailingOnly = TRUE)
window = as.numeric(args[1])

obs.dir    = paste0("/projects/dgatti/Dog/SimulatedData/PureBreedHQ/phased/shapeit/multi/win", window, "/")
emis.dir   = paste0("/projects/dgatti/Dog/InputData/PureBreedsAllHQ/shapeit/win", window, "/")
map.dir    = paste0("/projects/dgatti/Dog/InputData/PureBreedsAllHQ/shapeit/win", window, "/")
truth.dir  = "/projects/dgatti/Dog/SimulatedData/PureBreedHQ/phased/shapeit/multi/"
output.dir = paste0("/projects/dgatti/Dog/HMM/HMM_v7_SHAPEIT/win", window, "/")

fig.dir = "/projects/dgatti/Dog/Figures/window_size_sim/"

# Read in observations.
obs.files = dir(path = obs.dir, pattern = "_phased.rds$", full.names = T)

##########
#obs.files = obs.files[grep("chr1_phased.rds$", obs.files)]

obs = lapply(obs.files, readRDS)
chr = as.numeric(gsub(paste0("^", obs.dir, "/simulated_chr|\\_phased.rds$"), "", obs.files))
names(obs) = chr
obs = obs[order(chr)]

# Read in emission probs.
emis.files = dir(path = emis.dir, pattern = "emission_probs_chr", full.names = T)

##########
#emis.files = emis.files[grep("chr1\\.rds$", emis.files)]

emis = lapply(emis.files, readRDS)
emis = lapply(emis, log)
breeds = unique(sub("_[0-9]+$", "", rownames(emis[[1]])))
chr = as.numeric(gsub(paste0("^", emis.dir, "/emission_probs_chr|\\.rds$"), "", emis.files))
names(emis) = chr
emis = emis[order(chr)]

# Read in true simulation files.
truth.files = dir(path = truth.dir, pattern = "^truth", full.names = T)

#########
#truth.files = truth.files[grep("chr1.rds$", truth.files)]

truth = lapply(truth.files, readRDS)
chr = as.numeric(gsub(paste0("^", truth.dir, "/truth_chr|\\.rds$"), "", truth.files))
names(truth) = chr
truth = truth[order(chr)]


# Read in marker map.
map.files = dir(path = map.dir, pattern = "\\.map$", full.names = T)

#########
#map.files = map.files[grep("chr1\\.map$", map.files)]

map = lapply(map.files, read.delim, sep = "", header = F)
chr = as.numeric(gsub(paste0("^", map.dir, "/chr|\\.map$"), "", map.files))
names(map) = chr
map = map[order(chr)]

stopifnot(names(obs) == names(emis))
stopifnot(names(obs) == names(truth))
stopifnot(names(obs) == names(map))

# NOTE: for SHAPEIT, we have to subset the markers.
for(i in 1:length(obs)) {
  mkr = intersect(colnames(obs[[i]]), colnames(emis[[i]]))
  obs[[i]]  = obs[[i]][,mkr]
  emis[[i]] = emis[[i]][,mkr]
  truth[[i]] = truth[[i]][,mkr]
  map[[i]]  = map[[i]][map[[i]][,2] %in% mkr,]
} # for(i)

stopifnot(sapply(obs, ncol) == sapply(emis, ncol))
stopifnot(sapply(obs, ncol) == sapply(truth, ncol))
stopifnot(sapply(obs, ncol) == sapply(map, nrow))


# Set the transtion probs multiplier. This is a negative number that will
# be subtracted from the cM values on the log scale.
trans.factor = seq(-1, -200, -10)
num.dogs = nrow(obs[[1]])

# Subset the dogs, if desired.
if(num.dogs < nrow(obs[[1]])) {
  obs = lapply(obs, function(z) { z[1:num.dogs,] })
} # if(num.dogs < nrow(obs[[1]]))

dog.haps = rownames(obs[[1]])

num.recomb = matrix(0, length(trans.factor), length(dog.haps), dimnames =
             list(trans.factor, dog.haps))
true.num.recomb = rep(0, length(dog.haps))
names(true.num.recomb) = dog.haps
prop.correct = matrix(0, length(trans.factor), length(dog.haps), dimnames =
               list(trans.factor, dog.haps))

for(dog in 1:nrow(obs[[1]])) {

  # True number of recombinations.
  tnr = sapply(truth, function(z) { sum(abs(diff(z[dog,])) > 0) })
  true.num.recomb[dog] = sum(tnr) 

  # Run models with different trans.factors.
  for(j in 1:length(trans.factor)) {

    print(paste(dog, j))

    t1 = proc.time()[3]

    vit = vector("list", 39)

    for(chr in 1:39) {

#    vit = vector("list", 1)

#    for(chr in 1) {

#  probs[[j]] = forwardBackward(obs = as.integer(obs[dog,]), emis = emis, 
#               map = map[,3], mapfactor = trans.factor[j], window = window)
#  probs[[j]] = exp(probs[[j]])
#  dimnames(probs[[j]]) = list(breeds, colnames(emis))

      # Add one to the Viterbi result because C++ is 0-based and R is 1-based.
      vit[[chr]] = viterbi(obs = as.integer(obs[[chr]][dog,]), emis = emis[[chr]],
                    map = map[[chr]][,3], mapfactor = trans.factor[j], window = window)
      vit[[chr]][,1] = vit[[chr]][,1] + 1

    } # for(chr)

    # Split up the viterbit results into breeds calls (column 1) and probs
    # (column 2).
    v1 = lapply(vit, function(z) { z[,1] })
    v2 = lapply(vit, function(z) { z[,2] })

    # Number of recombinations.
    nr = sapply(v1, function(z) { sum(abs(diff(z)) > 0) })
    num.recomb[j,dog] = sum(nr)

    # Comparison to truth.
    prop.correct[j,dog] = mean(unlist(mapply(function(x, y) { x == y[dog,] }, x = v1, y = truth)))

    print(proc.time()[3] - t1)

  } # for(j)

} # for(dog)

save(num.recomb, prop.correct, file = paste0(fig.dir, "trans_prob_sims_window", window, ".Rdata"))

png(paste0(fig.dir, "trans_prob_sims_window", window, ".png"), width = 2400, 
    height = 1200, res = 200)
layout(matrix(1:2, 1, 2))
par(las = 1)
boxplot(t(num.recomb) - true.num.recomb, main = "Num. Recomb. - True Num. Recomb.")
boxplot(t(prop.correct), main = "Proportion of Chr 1 Correct")
dev.off()


