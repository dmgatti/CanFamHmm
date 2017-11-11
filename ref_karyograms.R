################################################################################
# Draw and save the karyograms for the reference dogs.
# Daniel Gatti
# dan.gatti@jax.org
# July 31, 2017
################################################################################
options(stringsAsFactors = FALSE)

src.dir = "/projects/dgatti/Dog/Scripts/CanFamHmm/"
source(paste0(src.dir, "dog_karyogram.R"))

window = 10

hmm.dir = paste0("/projects/dgatti/Dog/HMM/HMM_PureBreedHQ/win", window, "/")
map.dir = paste0("/projects/dgatti/Dog/InputData/PureBreedsAllHQ/shapeit/win", window, "/")
fam.dir = paste0("/projects/dgatti/Dog/ReferenceData/PureBreedsHQ/")

setwd(hmm.dir)

# Get the sex information for each dog.
fam.files = dir(path = fam.dir, pattern = "fam$", full.names = T)
fam = read.delim(fam.files[1], header = F, sep = " ")
sex = rep("F", nrow(fam))
sex[fam[,5] == 1] = "M"
names(sex) = fam[,1]
sex = sex[order(names(sex))]

# Get the breed order.
breeds = sort(unique(fam[,1]))

# Get the viterbi files.
viterbi.files = dir(pattern = paste0("viterbi_win", window, ".rds$"))

# Read in the map files.
map.files = dir(path = map.dir, pattern = "map$", full.names = T)
markers = lapply(map.files, read.delim, header = F, sep = "")
chr = strsplit(sub("\\.map$", "", map.files), "chr")
chr = as.numeric(sapply(chr, "[", 2))
names(markers) = chr
markers = markers[order(chr)]

# Subset the map to match the markers used by the model.
hap1 = readRDS(viterbi.files[1])
for(i in 1:length(markers)) {
  markers[[i]] = markers[[i]][markers[[i]][,2] %in% rownames(hap1[[i]]),]
} # for(i)
stopifnot(sapply(hap1, nrow) == sapply(markers, nrow))

map = lapply(markers, function(z) { z[,4] * 1e-6 })

dog.names = unique(gsub(paste0("_viterbi_win", window,".rds$"), "", viterbi.files))
dog.names = unique(sub("_(A|B)\\.", "_(A|B).", dog.names))
dog.names = unique(sub("_(A|B)\\.", "_(A|B).", dog.names))
dog.names = unique(sub("_(A|B)$", "_(A|B)", dog.names))


for(d in dog.names) {

  print(d)

  gr = grep(paste0("^", d, "_viterbi_win", window, ".rds$"), viterbi.files)

  stopifnot(length(gr) == 2)

  # Read in the haplotype files.
  hap1 = readRDS(viterbi.files[gr[1]])
  hap2 = readRDS(viterbi.files[gr[2]])

  png(paste0(sub("_\\(A\\|B\\)", "", d), "_karyogram.png"), width = 2000, height = 1000, res = 128)
  dog_karyogram(hapA = hap1, hapB = hap2, map = map, breeds = breeds,
                 title = d, sex = "F")
  dev.off()

} # for(i)

