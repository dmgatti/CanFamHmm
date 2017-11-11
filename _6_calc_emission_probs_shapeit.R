################################################################################
# Read in each SHAPEIT phased haplotype file and write out *.rds files for 
# input to the HMM and calculate emission probabilities.
# Daniel Gatti
# dan.gatti@jax.org
# June 14, 2017
################################################################################
options(stringsAsFactors = F)
library(Rcpp)
library(BiocParallel)
library(VariantAnnotation)

source.dir = "/projects/dgatti/Dog/Scripts/CanFamHmm/"
sourceCpp(paste0(source.dir, "HMM.cpp"))

window = 10
input.dir = "/projects/dgatti/Dog/ReferenceData/PureBreedsHQ/shapeit/"
output.dir = paste0("/projects/dgatti/Dog/InputData/PureBreedsAllHQ/shapeit/win",
             window,"/")
map.dir = "/projects/dgatti/Dog/ReferenceData/PureBreedsHQ/"
obs.dir = paste0("/projects/dgatti/Dog/GeneseekData/Cleaned/shapeit/")

# Command line arguments.
# Arg 1: input file directory.
# Arg 2: output file directory.
# Arg 3: map file directory.
# Arg 4: observation file directory with phased data, not window subsetted.
# Arg 5: window size.
args = commandArgs(trailingOnly = TRUE)

if(length(args) != 5) {
  stop(paste("5 arguments required: input.dir, output.dir, map.dir, obs.dir",
       "and window."))
} # if(length(args) != 5)

input.dir  = args[1]
output.dir = args[2]
map.dir = args[3]
obs.dir = args[4]
window = as.numeric(args[5])

#hap.file = paste0(input.dir, "PureBreedsHQ_chr38_phased.haps.gz")
#fam.file = paste0(input.dir, "PureBreedsHQ_chr38_phased.sample")


# Read in a VCF file and return a matrix to use in estimating the emission
# probabilities. 
# hap.file: character string containing the full path to the VCF file.
# Returns: matrix with samples in rows and markers in columns.
read_geno = function(hap.file, fam.file) {

  hap = read.table(hap.file, sep = "", header = F)
  fam = read.delim(fam.file, sep = " ", header = T)
  fam = fam[-1,]
  map = hap[,1:5]
  rownames(hap) = hap[,2]
  hap = as.matrix(hap[,-(1:5)])
  colnames(hap) = paste(rep(fam[,1], each = 2), LETTERS[1:2], sep = "_")

  # NOTE: We want to transpose the data to place haplotypes in rows and
  #       markers in columns.

  return(t(hap))

} # read_geno()



# Calculate the emission probabilities in a window of SNPs.
# Arguments:
# geno.file: String containing the full path to the genotype *.rds file.
# fam.file: String containing full path to the family information in PLINK format.
# obs.file: String containing full path to the observation files.
# window: integer that is the SNP window size.
calc_emission_probs = function(geno.file, fam.file, obs.file, window) {

  geno = readRDS(file = geno.file)
  fam  = read.delim(fam.file, header = T, sep = " ")
  fam = fam[-1,]
  obs = readRDS(obs.file)
  # Make sure that the reference genotype and the observations contain the same
  # markers.
  geno = geno[,colnames(obs)]

  # Make breed numbers. Start at 0 since C is 0-based.
  # We need to expand the breeds so that each occurs twice.
  hap.breeds = as.numeric(factor(fam[,1])) - 1
  hap.breeds = rep(hap.breeds, each = 2) # There are 2 haplotypes per sample.
  num.states = 2^window

  # Get the emission probs from C++. 
  emis = get_emission_probs(geno = geno, breeds = hap.breeds, window = window, 
                            error = min(0.001, 1.0 / num.states))

  # Set up the row and column names for the emission probs.
  # We need to repeat each breed num.states times.
  sorted.unique.breeds = levels(factor(fam[,1]))
  emis.breeds = rep(sorted.unique.breeds, each = num.states)
  emis.breeds = paste(emis.breeds, 0:(num.states - 1), sep = "_")
  window.markers = colnames(geno)[floor(window / 2):(ncol(geno) - ceiling(window / 2))]

  dimnames(emis) = list(emis.breeds, window.markers[1:ncol(emis)])

  return(emis)

} # calc_emission_probs()


# Read in the phased files and write out to *.rds files.
setwd(input.dir)

hap.files = dir(path = input.dir, pattern = "phased.haps.gz$", full.names = TRUE)
fam.files = dir(path = input.dir, pattern = "sample$", full.names = TRUE)
map.files = dir(path = map.dir, pattern = "map$", full.names = TRUE)

# Verify that they're the same length.
stopifnot(length(hap.files) == length(fam.files))
stopifnot(length(hap.files) == length(map.files))

for(i in 1:length(hap.files)) {

#  print(i)
#  geno = read_geno(hap.file = hap.files[i], fam.file = fam.files[i])
#  saveRDS(geno, file = sub("haps.gz", "rds", hap.files[i]))

} # for(i)


# Make objects to pass into calc_emission_probs().
# Note that, for SHAPEIT, we need to subset the markers to match those in
# the obseration set, which is produced in step 12.
setwd(input.dir)
geno.files = dir(pattern = "phased.rds$")
fam.files  = dir(pattern = "sample$")
obs.files  = dir(path = obs.dir, pattern = "_phased.rds", full.names = T)
chr = strsplit(geno.files, "chr")
chr = sapply(chr, "[", 2)
chr = strsplit(chr, split = "\\_")
chr = as.numeric(sapply(chr, "[", 1))

stopifnot(length(geno.files) == length(obs.files))
stopifnot(length(geno.files) == length(map.files))

for(i in 1:length(geno.files)) {

  print(chr[i])

  emis = calc_emission_probs(geno.file = geno.files[i], fam.file = fam.files[i],
                             obs.file = obs.files[i], window = window)

  saveRDS(emis, file = paste0(output.dir, "emission_probs_chr", chr[i], ".rds"))

} # for(i)

