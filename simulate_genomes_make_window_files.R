################################################################################
# Read in the simulated data and write out to *.rds files.
# Create haplotypes in a non-overlapping sliding window.
# Daniel Gatti
# dan.gatti@jax.org
# August 15, 2017
################################################################################
options(stringsAsFactors = FALSE)

window = 7
input.dir = "/projects/dgatti/Dog/SimulatedData/PureBreedHQ/phased/shapeit/multi/"
sim.dir   = "/projects/dgatti/Dog/SimulatedData/PureBreedHQ/phased/shapeit/multi/"
output.dir = paste0("/projects/dgatti/Dog/SimulatedData/PureBreedHQ/phased/shapeit/multi/win",
             window ,"/")

setwd(input.dir)

# Get the phased haplotype files.
hap.files = dir(pattern = "phased.rds$")

# We need these because SHAPEIT removed some markers with high no-call rates.
jax.files = dir(path = sim.dir, pattern = "phased.rds$", full.names = T)

num.states = 2^window
vec = 2^(0:(window - 1))

# Read in each haplotypefile, calculate the haplotypes for each sample in 
# the sliding window and write out to RDS files.
for(i in 1:length(hap.files)) {

  print(hap.files[i])

  hap = readRDS(hap.files[i])
  jax = readRDS(jax.files[i])
  hap = hap[,colnames(jax)]

  num.windows = ncol(hap) - window + 1
  window.markers = colnames(hap)[floor(window / 2):(ncol(hap) - ceiling(window / 2))]

  mat = matrix(0, nrow(hap), num.windows, dimnames = 
        list(rownames(hap), window.markers[1:num.windows]))

  hap = t(hap)

  for(j in 1:num.windows) {

    rng = j:(j + window - 1)
    mat[,j] = colSums(vec * hap[rng,])

  } # for(j)

  saveRDS(mat, file = paste0(output.dir, hap.files[i]))

} # for(i)


