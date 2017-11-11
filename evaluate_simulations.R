# Given the directory containing the simulations, evaluate the quality
# of the simulations by comparing to the truth.
options(stringsAsFactors = F)
library(tidyverse)

# Arg 1: simulated genome directory.
# Arg 2: HMM output directory.
# Arg 3: Map directory.
# Arg 4: Output directory.
# Arg 5: SNP window size.
# Arg 6: Trans. prob. factor.
args = commandArgs(trailingOnly = T)

window = 5
sim.dir = "/hpcdata/dgatti/Dog/SimulatedData/PureBreedHQ/phased/shapeit/multi/"
hmm.dir = paste0("/hpcdata/dgatti/Dog/SimulatedData/PureBreedHQ/HMM/shapeit/win", window, "/")
map.dir = paste0("/projects/dgatti/Dog/InputData/PureBreedsAllHQ/shapeit/win", window, "/")
out.dir = "/hpcdata/dgatti/Dog/SimulatedData/PureBreedHQ/HMM/shapeit/"
trans.fact = -91

sim.dir = args[1]
hmm.dir = args[2]
map.dir = args[3]
out.dir = args[4]
window  = as.numeric(args[5])
trans.fact = as.numeric(args[6])

# Get the truth files.
truth.files = dir(path = sim.dir, pattern = "^truth", full.names = T)

# Get the HMM viterbi files.
vit.files = dir(path = hmm.dir, pattern = "^sim", full.names = T)


# Read in the truth files.
truth = vector("list", length(truth.files))
for(i in 1:length(truth.files)) {

  truth[[i]] = readRDS(truth.files[i])

} # for(i)

# Sort truth files by chr.
gr = unlist(gregexpr("chr", truth.files))
truth.chr = substring(truth.files, gr, gr + 4)
truth.chr = as.numeric(sub("^chr|\\.$", "", truth.chr))
names(truth) = truth.chr
truth = truth[order(truth.chr)]


# Read in all of the simulated Viterbi files and put them in large matrices,
# one per chromosome.

vit = vector("list", length(truth.files))
names(vit) = names(truth.files)

prob = vector("list", length(truth.files))
names(prob) = names(truth.files)

# Get the sample IDs.
st = unlist(gregexpr("sim", vit.files))
en = unlist(gregexpr("viterbi", vit.files)) - 2
samples = substring(vit.files, st, en)

# Read in one sample to get the marker IDs.
tmp = readRDS(vit.files[i])

# Create the sample matrices for each chromosome.
for(i in 1:length(vit)) {

  vit[[i]]  = matrix(0, length(samples), ncol(truth[[i]]), dimnames = 
              list(samples, colnames(truth[[i]])))
  prob[[i]] = matrix(0, nrow(vit[[i]]), ncol(vit[[i]]), dimnames = 
              dimnames(vit[[i]]))

  vit[[i]]  = vit[[i]][,rownames(tmp[[i]])]
  prob[[i]] = prob[[i]][,rownames(tmp[[i]])]
  truth[[i]] = truth[[i]][,rownames(tmp[[i]])]

} # for(i)

# Verify that the nubmer of marker on each chromosome is the same 
# between Viterbi and truth.
stopifnot(sapply(vit, ncol) == sapply(truth, ncol))

# Read in the Viterbi files.
for(i in 1:length(vit.files)) {

  print(i)

  tmp = readRDS(vit.files[i])

  # Get the breed on each chr.
  vit.breed = sapply(tmp, function(z) { z[,1] })
  vit.prob  = sapply(tmp, function(z) { z[,2] })

  # Place breed call and probability in the list.
  for(j in 1:length(vit)) {

    vit[[j]][samples[i],]  = vit.breed[[j]]
    prob[[j]][samples[i],] = vit.prob[[j]]

  } # for(j)

} # for(i)

# Keep only the samples in Viterbi (to accomodate failed runs on the cluster).
m = match(rownames(vit[[1]]), rownames(truth[[1]]))
for(i in 1:length(truth)) {

  truth[[i]] = truth[[i]][m,]
  stopifnot(rownames(vit[[i]]) == rownames(truth[[i]]))

} # for(i)

stopifnot(sapply(truth,nrow) == sapply(vit,nrow))


# For each sample, calculate the proportion of correct calls.
total = sum(sapply(truth, ncol))
result = rep(0, length(samples))
names(result) = samples

for(i in 1:length(vit)) {

  result = result + rowSums(truth[[i]] == vit[[i]])

} # for(i)
result = result / total

saveRDS(result, file = paste0(out.dir, "breed_comp_win", window, ".rds"))

png(paste0(out.dir, "sim_prop_correct_win", window ,".png"), width = 1000, height = 1000, 
    res = 200)
ggplot(data.frame(Breed_Match = result), aes(Breed_Match)) + geom_density() +
   labs(xlab = "Proportion Correct", title = 
   paste0("Proportion of Simulated Genomes Correct: Window = ", window))
dev.off()


# Read in the maps.
map.files = dir(path = map.dir, pattern = "map$", full.names = T)
map = lapply(map.files, read.delim, sep = "", header = F)
names(map) = gsub(paste0(map.dir, "/chr|\\.map$"), "", map.files)
map = map[order(as.numeric(names(map)))]

for(i in 1:length(map)) {

  map[[i]] = map[[i]][match(colnames(truth[[i]]), map[[i]][,2]),]
  stopifnot(map[[i]][,2] == colnames(truth[[i]]))

} # for(i)


# Get the % correct by chunk size.
res = as.list(rep(0, 100))
names(res) = 1:100

for(i in 1:nrow(vit[[1]])) {

  print(i)

  vt = lapply(vit,   function(z) { z[i,] })
  tr = lapply(truth, function(z) { z[i,] })

  for(j in 1:length(tr)) {

    dz = diff(c(-1, tr[[j]]))
    dz = c(which(dz != 0), nrow(map[[j]]))
    blocks = cbind(map[[j]], tr[[j]], vt[[j]], tr[[j]] == vt[[j]])
    blocks[,1] = c(rep(1:(length(dz) - 1), diff(dz)), length(dz) - 1)
    blocks[,4] = blocks[,4] * 1e-6

    spl = split(blocks, blocks[,1])
    eq = sapply(spl, function(z) { 
           c(round(z[nrow(z),4] - z[1,4]), m = mean(z[,7]))
         })
    if(!is.matrix(eq)) {
      eq = matrix(eq, nrow = 2)
    }
    eq = eq[,eq[1,] > 0, drop = F]

    for(k in 1:ncol(eq)) {
      res[[eq[1,k]]] = c(res[[eq[1,k]]], eq[2,k])
    } # for(k)

  } # for(j)

} # for(i)

png(paste0(out.dir, "prop_correct_by_chunk_size_window", window,"_tf", abs(trans.fact) ,".png"),
    width = 2400, height = 1200, res = 200)
boxplot(res[1:30], main = paste("Window =", window, "TF =", trans.fact), col = "grey80")
dev.off()

