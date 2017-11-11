################################################################################
# SHAPEIT can't handle SNPs with 100% missing data in the JAX data set.
# Write out a list of marker positions with SNPs that SHAPEIT should exclude
# from the phasing.
# Daniel Gatti
# dan.gatti@jax.org
# Aug 2, 2017
################################################################################
options(stringsAsFactors = F)

# Arguments:
# CHR: the current chromosome.
# DATA_DIR: path to the directory where the JAX files are.

args = commandArgs(trailingOnly = TRUE)

if(length(args) != 2) {
  stop(paste("2 arguments required: CHR and REF_DIR"))
}

chr = args[1]
data.dir = args[2]

#chr = 1
#data.dir = "/projects/dgatti/Dog/GeneseekData/Cleaned/"

# Get the input PED and MAP filenames.
ped.file = dir(path = data.dir, pattern = paste0("chr", chr, ".ped.gz$"),
           full.names = T)
if(length(ped.file) != 1) {
  stop(paste("PED file for CHR", chr, "not found in data.dir", data.dir))
}

map.file = dir(path = data.dir, pattern = paste0("chr", chr, ".map$"),
              full.names = T)
if(length(map.file) != 1) {
  stop(paste("MAP file for CHR", chr, "not found in data.dir", data.dir))
}

# Read in the data.
ped = read.delim(ped.file, header = F, sep = "")
map = read.delim(map.file, header = F, sep = "")

stopifnot(2 * nrow(map) == ncol(ped) - 6)

# Get the columns in PED with missing data and convert over to the 
# corresponding MAP file rows. The PED file has 6 header columns and 2
# columns for each marker.
wh = which(colMeans(ped == 0) > 0.9999)
wh = wh - 6 # 6 header columns.
wh = wh[wh > 0]
wh = matrix(wh, nrow = 2)
wh = wh[2,] / 2

write(map[wh,4], file = paste0(data.dir, "jax_chr", chr, ".exclude"), 
      ncolumns = 1, sep = "\n")


