################################################################################
# Draw and save the karyograms for the JAX dogs.
# Daniel Gatti
# dan.gatti@jax.org
# July 31, 2017
################################################################################
options(stringsAsFactors = FALSE)

src.dir = "/projects/dgatti/Dog/Scripts/CanFamHmm/"
source(paste0(src.dir, "dog_karyogram.R"))

window = 10

hmm.dir = paste0("/projects/dgatti/Dog/HMM/HMM_v7_SHAPEIT/win", window, "/")
map.dir = paste0("/projects/dgatti/Dog/GeneseekData/Cleaned/shapeit/win", window,
          "/")
fam.dir = "/projects/dgatti/Dog/GeneseekData/Cleaned/"
ref.dir = "/projects/dgatti/Dog/ReferenceData/PureBreedsHQ/shapeit/"

setwd(hmm.dir)

# Get the sex information for each dog.
fam.files = dir(path = fam.dir, pattern = "fam$", full.names = T)
fam = read.delim(fam.files[1], header = F, sep = " ")
sex = rep("F", nrow(fam))
sex[fam[,5] == 1] = "M"
names(sex) = fam[,1]
sex = sex[order(names(sex))]

# Get the breed order from the reference.
sam = read.delim(paste0(ref.dir, "PureBreedsHQ_chr38_phased.sample"), sep = "")[-1,]
breeds = sort(unique(sam[,1]))

# Read in the map files.
map.files = dir(path = map.dir, pattern = "map$", full.names = T)
markers = lapply(map.files, read.delim, header = F, sep = "")
chr = strsplit(sub("\\.map$", "", map.files), "chr")
chr = as.numeric(sapply(chr, "[", 2))
names(markers) = chr
markers = markers[order(chr)]
map = lapply(markers, function(z) { z[,4] * 1e-6 })
for(i in 1:length(map)) {
  names(map[[i]]) = markers[[i]][,2]
} # for(i)

# Get the viterbi files.
viterbi.files = dir(pattern = "viterbi.rds$")

dog.names = unique(gsub("_(A|B)_HMM_viterbi.rds$", "", viterbi.files))

# Subset the markers to accomodate the sliding window.
hap1 = readRDS(viterbi.files[1])
stopifnot(names(hap1) == names(map))
for(i in 1:length(hap1)) {

  map[[i]] = map[[i]][rownames(hap1[[i]])]

} # for(i)
stopifnot(sapply(map, length) == sapply(hap1, nrow))

for(d in dog.names) {

  print(d)

  gr = grep(paste0("^", d, "_(A|B)"), viterbi.files)

  stopifnot(length(gr) == 2)

  # Read in the haplotype files.
  hap1 = readRDS(viterbi.files[gr[1]])
  hap2 = readRDS(viterbi.files[gr[2]])

  # For now, strip out the probs in column 2, we may use them later.
  probA = lapply(hap1, function(z) { z[,2] })
  probB = lapply(hap2, function(z) { z[,2] })
  hapA  = lapply(hap1, function(z) { z[,1] })
  hapB  = lapply(hap2, function(z) { z[,1] })

  # Replace low frequency breed haplotypes.
#  retval = simplify_haps(hapA = hapA, hapB = hapB)
#  hap1 = retval$hapA
#  hap2 = retval$hapB

  # Reassemble haps and probs.
#  hap1 = mapply(cbind, hap1, probA)
#  hap2 = mapply(cbind, hap2, probB)

  png(paste0(d, "_karyogram.png"), width = 2000, height = 2600, res = 200)
  dog_karyogram(hapA = hap1, hapB = hap2, map = map, breeds = breeds,
                title = d, sex = sex[d])
  dev.off()

  png(paste0(d, "_barplot.png"), width = 2000, height = 1200, res = 200)
  dog_breed_barplot(hapA = hap1, hapB = hap2, breeds = breeds, title = d)
  dev.off()

} # for(i)

