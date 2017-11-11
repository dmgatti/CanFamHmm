# Adjust the transition probs and look at the viterbi and F-B output.
options(stringsAsFactors = F)
library(Rcpp)

src.dir    = "/projects/dgatti/Dog/Scripts/CanFamHmm/"
sourceCpp(paste0(src.dir, "HMM.cpp"))

#window = 1
args = commandArgs(trailingOnly = TRUE)
window = as.numeric(args[1])

obs.dir    = paste0("/projects/dgatti/Dog/SimulatedData/PureBreedHQ/phased/shapeit/multi/win", window, "/")
emis.dir   = paste0("/projects/dgatti/Dog/InputData/PureBreedsAllHQ/shapeit/win", window, "/")
map.dir    = paste0("/projects/dgatti/Dog/InputData/PureBreedsAllHQ/shapeit/win", window, "/")
output.dir = paste0("/projects/dgatti/Dog/HMM/HMM_v7_SHAPEIT/win", window, "/")

fig.dir = "/projects/dgatti/Dog/Figures/window_size_sim/"

# Read in observations.
obs.files = dir(path = obs.dir, pattern = "_phased.rds$", full.names = T)

##########
obs.files = obs.files[grep("chr1_phased.rds$", obs.files)]

obs = lapply(obs.files, readRDS)
chr = as.numeric(gsub(paste0("^", obs.dir, "/simulated_chr|\\_phased.rds$"), "", obs.files))
names(obs) = chr
obs = obs[order(chr)]

# Read in emission probs.
emis.files = dir(path = emis.dir, pattern = "emission_probs_chr", full.names = T)

##########
emis.files = emis.files[grep("chr1\\.rds$", emis.files)]

emis = lapply(emis.files, readRDS)
emis = lapply(emis, log)
breeds = unique(sub("_[0-9]+$", "", rownames(emis[[1]])))
chr = as.numeric(gsub(paste0("^", emis.dir, "/emission_probs_chr|\\.rds$"), "", emis.files))
names(emis) = chr
emis = emis[order(chr)]

# Read in marker map.
map.files = dir(path = map.dir, pattern = "\\.map$", full.names = T)

#########
map.files = map.files[grep("chr1\\.map$", map.files)]

map = lapply(map.files, read.delim, sep = "", header = F)
chr = as.numeric(gsub(paste0("^", map.dir, "/chr|\\.map$"), "", map.files))
names(map) = chr
map = map[order(chr)]

stopifnot(names(obs) == names(emis))
stopifnot(names(obs) == names(map))

# NOTE: for SHAPEIT, we have to subset the markers.
for(i in 1:length(obs)) {
  mkr = intersect(colnames(obs[[i]]), colnames(emis[[i]]))
  obs[[i]]  = obs[[i]][,mkr]
  emis[[i]] = emis[[i]][,mkr]
  map[[i]]  = map[[i]][map[[i]][,2] %in% mkr,]
} # for(i)

stopifnot(sapply(obs, ncol) == sapply(emis, ncol))
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
num.prop90 = matrix(0, length(trans.factor), length(dog.haps), dimnames =
             list(trans.factor, dog.haps))
model.prob = matrix(0, length(trans.factor), length(dog.haps), dimnames =
             list(trans.factor, dog.haps))

for(dog in 1:nrow(obs[[1]])) {

  # Run models with different trans.factors.
  for(j in 1:length(trans.factor)) {

    print(paste(dog, j))

    t1 = proc.time()[3]

#    vit = vector("list", 39)

#    for(chr in 1:39) {

    vit = vector("list", 1)

    for(chr in 1) {

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

    tbl = sort(table(unlist(v1)), decreasing = T)
    # Nubmer of breeds that make up 90% of the genome.
    num.prop90[j,dog] = min(which(cumsum(tbl) / sum(tbl) >= 0.9))
    # Proportion of genome with prob > 0.95.
    model.prob[j,dog] = mean(unlist(v2) > 0.95)

    print(proc.time()[3] - t1)

  } # for(j)

} # for(dog)

save(num.recomb, num.prop90, model.prob, file = paste0(fig.dir, "trans_prob_window", window, ".Rdata"))

png(paste0(fig.dir, "trans_prob_window", window, ".png"), width = 2000, 
    height = 800, res = 128)
layout(matrix(1:3, 1, 3))
par(las = 1)
boxplot(t(num.recomb), main = "Num. Recomb.")
boxplot(t(num.prop90), main = "Num. Breeds covering 90% of Chr 1")
boxplot(t(model.prob), main = "Prop. of Chr 1 w/ Prob. > 0.99.")
dev.off()

