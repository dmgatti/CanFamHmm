################################################################################
# Read in the SHAPEIT phased files for the JAX data and write out to *.rds files.
# Create haplotypes in a overlapping sliding window.
# Daniel Gatti
# dan.gatti@jax.org
# August 15, 2017
################################################################################
options(stringsAsFactors = FALSE)

window = 10
input.dir = "/projects/dgatti/Dog/GeneseekData/Cleaned/shapeit/"
jax.dir = "/projects/dgatti/Dog/GeneseekData/Cleaned/shapeit/"
hmm.dir = paste0("/projects/dgatti/Dog/InputData/PureBreedsAllHQ/shapeit/win",
          window ,"/")
output.dir = paste0("/projects/dgatti/Dog/GeneseekData/Cleaned/shapeit/win",
             window ,"/")

setwd(input.dir)

# Get the phased haplotype files.
hap.files = dir(pattern = "phased.rds$")

# Get the map files.
map.files = dir(pattern = "map$")

# Order the files by chromosome.
h.chr = regexpr("[0-9]", hap.files)
h.chr = substring(hap.files, h.chr, h.chr + 1)
h.chr = as.numeric(sub("_", "", h.chr))
hap.files = hap.files[order(h.chr)]

m.chr = regexpr("[0-9]", map.files)
m.chr = substring(map.files, m.chr, m.chr + 1)
m.chr = as.numeric(sub("\\.", "", m.chr))
map.files = map.files[order(m.chr)]

# We need these because SHAPEIT removed some markers with high no-call rates.
jax.files = dir(path = jax.dir, pattern = "phased.rds$", full.names = T)

# Verify that the chromosomes are the same in both set of files.
j.chr = regexpr("chr", jax.files)
h.chr = regexpr("chr", hap.files)
j.chr = substring(jax.files, j.chr, j.chr + 5)
h.chr = substring(hap.files, h.chr, h.chr + 5)
j.chr = as.numeric(gsub("[a-z]+|_", "", j.chr))
h.chr = as.numeric(gsub("[a-z]+|_", "", h.chr))
jax.files = jax.files[order(j.chr)]

num.states = 2^window
vec = 2^(0:(window - 1))

# Read in each haplotypefile, calculate the haplotypes for each sample in 
# the sliding window and write out to RDS files.
for(i in 1:length(hap.files)) {

  print(hap.files[i])

  hap = readRDS(hap.files[i])
  map = read.delim(map.files[i], sep = "", header = F)
  jax = readRDS(jax.files[i])
  hap = hap[,colnames(jax)]

  # Non-overlapping windows.
#  num.windows = floor(ncol(hap) / window)
#  window.markers = colnames(hap)[seq(floor(window / 2), ncol(hap), window)]

  # Overlapping windows.
  num.windows = ncol(hap) - window + 1
  window.markers = colnames(hap)[floor(window / 2):(ncol(hap) - ceiling(window / 2))]

  map = map[map[,2] %in% window.markers,]
  write.table(map, file = paste0(output.dir, map.files[i]), sep = " ", row.names = F, 
              col.names = F, quote = F)

  mat = matrix(0, nrow(hap), num.windows, dimnames = 
        list(rownames(hap), window.markers))

  hap = t(hap)

  for(j in 1:num.windows) {

    rng = j:(j + window - 1)
    mat[,j] = colSums(vec * hap[rng,])

  } # for(j)

  saveRDS(mat, file = paste0(output.dir, hap.files[i]))

} # for(i)


