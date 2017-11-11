################################################################################
# Calculate an initial transition probability matrix using a kinship matrix.
# Daniel Gatti
# dan.gatti@jax.org
# June 19, 2017
################################################################################
options(stringsAsFactors = F)
library(BiocParallel)

input.dir  = "/hpcdata/dgatti/Dog/ReferenceData/PureBreedsHQ/"
output.dir = "/hpcdata/dgatti/Dog/InputData/PureBreedsAllHQ/"

# Function to read in the genotypes on one chromosome and calculate a
# kinship matrix.
calc_kinship = function(geno.file) {

  geno = readRDS(geno.file)
  geno = t(geno)
  K = matrix(0, ncol(geno), ncol(geno), dimnames = list(
      colnames(geno), colnames(geno)))
  for(i in 1:ncol(geno)) {

    if(i %% 100 == 0) print(paste(i, "of", ncol(geno)))
    K[i,] = colMeans(geno[,i] == geno)

  } # for(i)

  return(K)

} # calc_kinship()

geno.files = dir(path = input.dir, pattern = ".phased.rds$", 
             full.names = T)
chr = strsplit(geno.files, "chr")
chr = sapply(chr, "[", 2)
chr = strsplit(chr, "\\.")
chr = as.numeric(sapply(chr, "[", 1))
geno.files = geno.files[order(chr)]
chr = sort(chr)

setwd(output.dir)

param = MulticoreParam(workers = 5, progressbar = TRUE, log = TRUE)

K = bplapply(X = geno.files, FUN = calc_kinship, BPPARAM = param)
names(K) = chr

for(i in 1:length(K)) {

  saveRDS(K[[i]], file = paste0(output.dir, "kinship_chr", chr[i], ".rds"))

} # for(i)


# Make kinship plots.
for(i in 1:length(chr)) {

  print(i)

  K = readRDS(file = paste0(output.dir, "kinship_chr", chr[i], ".rds"))

  png(paste0(output.dir, "kinship_chr", chr[i], ".png"), width = 3000, height = 3000)
  breaks = 0:100/100
  col = heat.colors(length(breaks) - 1)
  image(1:nrow(K), 1:ncol(K), K, ann = F, axes = F)
  s = seq(1, nrow(K), 5)
  mtext(side = 1, line = 0.2, at = s, text = rownames(K)[s], las = 2)
  mtext(side = 2, line = 0.2, at = s, text = rownames(K)[s], las = 2)
  mtext(side = 3, line = 0.2, at = s, text = rownames(K)[s], las = 2)
  mtext(side = 4, line = 0.2, at = s, text = rownames(K)[s], las = 2)
  dev.off()

} # for(i)


# Condense sample kinship matrix down to a breed kinship matrix, one per chromosome.
fam = read.delim(dir(path = input.dir, pattern = "\\.fam$", full.names = T)[1],
      sep = "", header = F)
breeds = sort(unique(fam[,1]))

for(i in 1:length(chr)) {

  print(i)

  K = readRDS(file = paste0(output.dir, "kinship_chr", chr[i], ".rds"))
  breedK = matrix(0, length(breeds), length(breeds), dimnames = list(breeds, breeds))

  for(j in 1:length(breeds)) {

    b1 = grep(paste0("^", breeds[j]), rownames(K))
    for(k in 1:length(breeds)) {

      b2 = grep(paste0("^", breeds[k]), rownames(K))
      breedK[j,k] = mean(K[b1,b2][upper.tri(K[b1,b2])])

    } # for(k)
  } # for(j)

  saveRDS(breedK, paste0(output.dir, "kinship_breed_chr", chr[i], ".rds"))
  
} # for(i)



# Working with PLINK relationship matrix.
breedD = matrix(0, length(breeds), length(breeds), dimnames = list(breeds, breeds))

for(j in 1:length(breeds)) {

    b1 = grep(paste0("^", breeds[j]), rownames(d))
    for(k in 1:length(breeds)) {

      b2 = grep(paste0("^", breeds[k]), rownames(d))
      if(j == k) {
        breedD[j,k] = mean(d[b1,b2][upper.tri(d[b1,b2])])
      } else {
        breedD[j,k] = mean(d[b1,b2])
      }

    } # for(k)
  } # for(j)

