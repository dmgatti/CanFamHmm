################################################################################
# Intersect the JAX and reference markers and create PLINK files for each 
# chromosome.
# Daniel Gatti
# dan.gatti@jax.org
# July 16, 2017
################################################################################
options(stringsAsFactors = F)
setwd("/projects/dgatti/Dog/")
library(ggplot2)

ref.dir = "/projects/dgatti/Dog/ReferenceData/PureBreedsHQ/"
in.dir  = "/projects/dgatti/Dog/GeneseekData/Cleaned/"
out.dir = "/projects/dgatti/Dog/GeneseekData/Cleaned/"
source.dir = "/projects/dgatti/Dog/Scripts/CanFamHmm/"

source(paste0(source.dir, "write_plink.R"))

# Read in the JAX data.
jax = readRDS(paste0(in.dir, "jax_cleaned.rds"))
jax.map = readRDS(paste0(in.dir, "jax_markers.rds"))
jax.hdr = jax[,1:6]
jax = jax[,-(1:6)]

# Make sample IDs unique.
jax.hdr[,2] = make.names(jax.hdr[,2])
jax.hdr[,2] = make.unique(jax.hdr[,2])
jax.hdr[,2] = sub("\\.1$", "2", jax.hdr[,2])
jax.hdr[,1] = jax.hdr[,2]

# Read in the Reference pure breed HQ markers.
map.files = dir(path = ref.dir, pattern = "map$", full.names = T)
map = lapply(map.files, read.delim, header = F)
map.chr = gsub("^/projects/dgatti/Dog/ReferenceData/PureBreedsHQ//PureBreedsHQ_|\\.map$",
          "", map.files)
names(map) = map.chr
chr.num = as.numeric(sub("chr", "", map.chr))
map = map[order(chr.num)]

# Split up the JAX data and markers by chromosome and write them out in PLINK 
# format.
for(i in 1:length(map)) {

  print(i)

  jax.ss = jax[,map[[i]][,2]]
  stopifnot(ncol(jax.ss) == nrow(map[[i]]))

  write_plink(fam = jax.hdr, data = jax.ss, map = map[[i]], 
              file.prefix = paste0(out.dir, "jax_", names(map)[i]))

} # for(i)

stop()

########################################################################
# Now run plink_ped2bed.sh on these files and then run the code below.
########################################################################

# Synch up the reference and alternate alleles in the reference and JAX data.
ref.bim.files = dir(path = ref.dir, pattern = "bim$", full.names = TRUE)
jax.bim.files = dir(path = out.dir, pattern = "bim$", full.names = TRUE)

stopifnot(length(ref.bim.files) == length(jax.bim.files))

for(i in 1:length(ref.bim.files)) {

  spl = strsplit(ref.bim.files[i], "_")[[1]]
  print(spl[2])

  ref = read.delim(ref.bim.files[i], header = F)
  jax = read.delim(jax.bim.files[i], header = F)
  stopifnot(nrow(ref) == nrow(jax))
  stopifnot(ncol(ref) == ncol(jax))

  stopifnot(ref[,6] == jax[,6])
  jax[,5] = ref[,5]
  write.table(jax, file = jax.bim.files[i], sep = "\t", row.names = F, 
              col.names = F, quote = F)

} # for(i)

########################################################################
# Now run plink_bed2vcf.sh.
########################################################################

