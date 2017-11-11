################################################################################################
# Gather the HMM results for each breed and look for mis-assignments that are in common across
# replicates within a breed.
# Daniel Gatti
# dan.gatti@jax.org
# June 22, 2017
################################################################################################
options(stringsAsFactors = FALSE)
setwd("/data/dgatti/Dog/")
library(tidyverse)

window = 1
hmm.dir = paste0("/projects/dgatti/Dog/HMM/HMM_PureBreedHQ/win", window, "/")
fam.dir = "/projects/dgatti/Dog/ReferenceData/PureBreedsHQ/"
out.dir = "/projects/dgatti/Dog/Analysis/PureBreedsHQ/"

# Read in the Chr 38 data to get the breed information.
fam = read.delim(paste0(fam.dir, "PureBreedsHQ_chr38.fam"), header = F, sep = " ")
rownames(fam) = paste(fam[,1], fam[,2], sep = "_")

# Get the unique breeds.
breeds = sort(unique(fam[,1]))

# For each breeds, read in the viterbi results.
setwd(hmm.dir)
files = dir(path = hmm.dir, pattern = "viterbi")

# Read in one file to get the dimensions.
f = readRDS(files[1])
len = length(unlist(f))

# Read in the markers.
map.files = dir(path = fam.dir, pattern = "map$", full.names = TRUE)
chrlen = rep(0, length(map.files))
markers = NULL
for(i in 1:length(map.files)) {

  map = read.delim(map.files[i], header = F)
  chrlen[i] = nrow(map)
  markers = rbind(markers, map)

} # for(i)
markers = markers[order(markers[,1]),]
names(chrlen) = gsub("^/projects/dgatti/Dog/ReferenceData/PureBreedsHQ//PureBreedsHQ_chr|\\.map$",
                     "", map.files)
chrlen = chrlen[order(as.numeric(names(chrlen)))]
chrmid = chrlen / 2 + cumsum(c(0, chrlen[-length(chrlen)]))

breaks = 0:100/100
col = rev(gray(breaks[-length(breaks)]))

# Read in the files and calculate the percentage correct.
vit.files = dir(pattern = "viterbi_win[0-9]+.rds$")
mean.correct = rep(0, length(vit.files))

for(i in 1:length(vit.files)) {

  print(i)

  vit = readRDS(vit.files[[i]])
  vit = unlist(sapply(vit, function(z) { z[,1] }))
  tbl = table(vit) / length(vit)
  mean.correct[i] = max(tbl)

} # for(i)

png(paste0("/projects/dgatti/Dog/Figures/ref_het/ref_correct_win", window, ".png"),
    width = 2400, height = 1200, res = 300)
ggplot(data.frame(Mean_Correct = mean.correct), aes(Mean_Correct)) +
    geom_density() +
    labs(title = paste0("Pure Breeds, Proportion Correct, Window = ", window))
dev.off()







# For each breed, read in all of the haplotypes and tabulate the proportion of the genome
# that haplotypes at the correct breed.
for(i in 1:length(breeds)) {

  print(breeds[i])

  breed.files = dir(pattern = paste0("^", breeds[i], ".+viterbi"))

  mat = matrix(0, nrow = length(breeds), ncol = len, dimnames = 
        list(sort(breeds), markers[,2]))

  for(j in 1:length(breed.files)) {

    f = readRDS(breed.files[j])
    f = unlist(f)
    pos.mat = cbind(f, 1:ncol(mat))
    mat[pos.mat] = mat[pos.mat] + 1

  } # for(j)

  mat = mat / length(breed.files)

  png(paste0(out.dir, breeds[i], "_HMM_mean_call.png"), height = 2400, width = 4000, res = 128)

  layout(matrix(1:2, 2, 1), heights = c(0.8, 0.2))
  par(plt = c(0.11, 0.89, 0.05, 0.92))
  image(1:ncol(mat), 1:nrow(mat), t(mat), axes = F, ann = F, main = breeds[i],
        breaks = breaks, col = col)
  mtext(side = 1, line = 1, at = chrmid, text = names(chrlen), font = 2, cex = 2)
  mtext(side = 3, line = 0.5, at = chrmid, text = names(chrlen), font = 2, cex = 2)
  mtext(side = 3, line = 2.5,  text = breeds[i], font = 2, cex = 2)
  mtext(side = 2, line = 0.25, at = 1:nrow(mat), text = rownames(mat), las = 2)
  mtext(side = 4, line = 0.25, at = 1:nrow(mat), text = rownames(mat), las = 2)
  abline(v = cumsum(chrlen), col = "grey80")
  abline(h = 1:length(breeds) + 0.5, col = "grey80")

  wh = which.max(rowSums(mat))
  par(plt = c(0.11, 0.89, 0.05, 0.99))
  plot(mat[wh,], type = "l", las = 2, xaxt = "n", xlab = "", ylab = "Proportion",
      xaxs = "i")
  abline(v = cumsum(chrlen), col = "grey80")

  dev.off()

  saveRDS(mat, file = paste0(out.dir, breeds[i], "_HMM_mean_call.rds"))

} # for(i)

