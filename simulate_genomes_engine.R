################################################################################
# Simulate genomes from teh reference data.
# The caller may set the chromosome, number of samples to simulate, the chunk
# size in Mb and the reference sample directory. This will allow easy switches
# between phasing algorithms. 
# NOTE: At the moment, we only simulate a single block size across the whole 
# genome. I will add support for a distribution of sizes later.
# Daniel Gatti
# dan.gatti@jax.org
# Aub 2, 2017
################################################################################
options(stringsAsFactors = FALSE)
source("/projects/dgatti/Dog/Scripts/CanFamHmm/simulate_genome.R")

# Arg 1: is the number of samples to simulate.
# Arg 2: Size in Mb of the chunks to simulate.
# Arg 3: is TRUE if the simulated alleles should be scrambled, FALSE otherwise.
# Arg 4: Full path to phased reference directory containing the *.rds files.
# Arg 5: Full path to MAP/FAM directory.
# Arg 6: Full path to output directory.
args = commandArgs(trailingOnly = TRUE)

if(length(args) != 6) {
  stop(paste("6 arguments required: n (numeric), chunk size in Mb,",
       "phased (TRUE or FALSE), reference path, map path and output path."))
} # if(length(args) != 6)

n   = as.numeric(args[1])
chunk = as.numeric(args[2])
phased = as.logical(args[3])
ref.dir = args[4]
map.dir = args[5]
out.dir = args[6]

# n = 1000
# chunk = c(1, 2:20)
# phased = TRUE
# ref.dir = "/projects/dgatti/Dog/ReferenceData/PureBreedsHQ/shapeit/"
# map.dir = "/projects/dgatti/Dog/ReferenceData/PureBreedsHQ/"
# out.dir = "/projects/dgatti/Dog/SimulatedData/PureBreedHQ/phased/shapeit/multi/"

for(chr in 1:39) {

  print(chr)

  # Get file names for the current chromosome.
  hap.file = dir(path = ref.dir, pattern = paste0("chr", chr, "_phased.rds$"),
             full.names = T)
  map.file = dir(path = map.dir, pattern = paste0("chr", chr, ".map$"),
             full.names = T)
  fam.file = dir(path = map.dir, pattern = paste0("chr", chr, ".fam$"),
             full.names = T)

  # Simulate the genomes.
  sim = simulate_genome(n = n, chunk_sizes = chunk,  
        hap.file = hap.file, map.file = map.file, fam.file = fam.file,
        phased = phased)

  # Add a PLINK header to the observations. We're making all females right now.
  obs = sim$obs
  hdr = data.frame(rownames(obs), rownames(obs), rep(0, nrow(obs)), 
        rep(0, nrow(obs)), rep(2, nrow(obs)), rep(-9, nrow(obs)))

  # Save the results.
  if(phased) {

    # Save as RDS for HMM.
    saveRDS(obs, file = paste0(out.dir, "simulated_chr", chr, "_phased.rds"))
    write.table(hdr, file = paste0(out.dir, "simulated_chr", chr, "_phased.fam"),
                sep = "\t")
    saveRDS(sim$truth, file = paste0(out.dir, "truth_chr", chr, ".rds"))

  } else {

  # Save as PLINK text file for phasing software.
  ##########################################
  stop("Need to write out to PLINK format!")
  ##########################################

  } # else

  saveRDS(sim$breeds, file = paste0(out.dir, "breeds.rds"))

} # for(chr)

