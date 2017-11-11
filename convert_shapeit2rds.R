################################################################################
# Convert the phased JAX haplotype files to *.rds files for the HMM.
# Also write out new MAP files containing the reduced SNP set.
#
# Daniel Gatti
# dan.gatti@jax.org
# Aug 2, 2017
################################################################################
options(stringsAsFactors = F)

# Arguments:
# CHR: the current chromosome.
# DATA_DIR: path to the directory where the JAX phased hap files are.
# MAP_DIR: path to the directory where the MAP files are.
args = commandArgs(trailingOnly = TRUE)

if(length(args) != 3) {
  stop(paste("3 arguments required: CHR, DATA_DIR and MAP_DIR"))
}

chr = args[1]
data.dir = args[2]
map.dir  = args[3]

#chr = 1
#data.dir = "/projects/dgatti/Dog/GeneseekData/Cleaned/shapeit/"
#map.dir = "/projects/dgatti/Dog/GeneseekData/Cleaned/"

# Get the input HAP and SAMPLE filenames.
hap.file = dir(path = data.dir, pattern = paste0("chr", chr, "_phased.haps$"),
           full.names = T)
if(length(hap.file) != 1) {
  stop(paste("HAP file for CHR", chr, "not found in data.dir", data.dir))
}

map.file = dir(path = map.dir, pattern = paste0("chr", chr, ".map$"),
              full.names = T)
if(length(map.file) != 1) {
  stop(paste("MAP file for CHR", chr, "not found in map.dir", map.dir))
}

sample.file = dir(path = data.dir, pattern = paste0("chr", chr, "_phased.sample$"),
              full.names = T)
if(length(sample.file) != 1) {
  stop(paste("SAMPLE file for CHR", chr, "not found in data.dir", data.dir))
}

# Read in the data.
hap = read.delim(hap.file, header = F, sep = "")
map = read.delim(map.file, header = F, sep = "")
sample = read.delim(sample.file, header = T, sep = "")
sample = sample[-1,]

stopifnot(2 * nrow(sample) == ncol(hap) - 5)

print("Writing REFERENCE files...")

# HAP file.
hdr = hap[,1:5]
hap = hap[,-(1:5)]
dimnames(hap) = list(hdr[,2], paste(rep(sample[,1], each = 2), LETTERS[1:2],
                sep = "_"))
new.hap.filename = sub("haps$", "rds", hap.file)
saveRDS(t(hap), file = new.hap.filename)

# MAP file.
new.map.filename = paste0(data.dir, "jax_chr", chr, ".map")
map = map[map[,2] %in% rownames(hap),]
write.table(map, file = new.map.filename, sep = " ", row.names = F,
            col.names = F, quote = F)
