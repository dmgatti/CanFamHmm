################################################################################
# Read in each Beagle VCF and write out *.rds files for input to the HMM and
# calculate emission probabilities.
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

# Command line arguments.
# Arg 1 is input file directory.
# Arg 2 is output file directory.
# Arg 3 is window size.
input.dir  = args[1]
output.dir = args[2]
window = as.numeric(args[3])

input.dir = "/projects/dgatti/Dog/ReferenceData/PureBreedsHQ/"
output.dir = "/projects/dgatti/Dog/InputData/PureBreedsAllHQ/"
window = 1

vcf.file = paste0(input.dir, "PureBreedsHQ_chr38.phased.vcf.gz")
fam.file = paste0(input.dir, "PureBreedsHQ_chr38.fam")
map.file = paste0(input.dir, "PureBreedsHQ_chr38.map")


# Read in a VCF file and return a matrix to use in estimating the emission
# probabilities. Tragically, the VCF sample IDs are the pasted
# vcf.file: character string containing the full path to the VCF file.
# Returns: matrix with samples in rows and markers in columns.
read_geno = function(vcf.file, fam.file, map.file) {

  vcf = readVcf(file = vcf.file, genome = "CanFam3.1")
  fam = read.delim(fam.file, sep = " ", header = F)
  map = read.delim(map.file, header = F)

  # Extract genotypes.
  vcf  = expand(vcf)
  geno = geno(vcf)$GT

  # NOTE: We want to transpose the data to place haplotypes in rows and
  #       markers in columns.

  g2 = matrix(0, nrow = 2 * ncol(geno), ncol = nrow(geno), dimnames = 
       list(paste(rep(colnames(geno), each = 2), rep(c("A", "B"), 
       ncol(geno)), sep = "_"), rownames(geno)))

  for(i in 1:nrow(geno)) {

    spl = as.numeric(unlist(strsplit(geno[i,], "\\|")))
    g2[,i] = spl

  } # for(i)

  return(g2)

} # read_geno()




# Calculate the emission probabilities in a window of SNPs.
# Arguments:
# geno.file: String containing the full path to the genotype *.rds file.
# fam.file: String containing full path to the family information in PLINK format.
# window: integer that is the SNP window size.
calc_emission_probs = function(geno.file, fam.file, window) {

  geno = readRDS(file = geno.file)
  fam  = read.delim(fam.file, header = F, sep = " ")

  # Make breed numbers. Start at 0 since C is 0-based.
  # We need to expand the breeds so that each occurs num.states times.
  breeds = as.numeric(factor(fam[,1])) - 1
  num.states = 2^window
  breeds = rep(breeds, each = num.states)

  # Get the emission probs from C++.
  emis = get_emission_probs(geno = geno, breeds = breeds, window = window, 
                            error = 0.01)

  # Set up the row and column names for the emission probs.
  # We need to repeat each breed num.states times.
  sorted.unique.breeds = levels(factor(fam[,1]))
  emis.breeds = rep(sorted.unique.breeds, each = num.states)
  emis.breeds = paste(emis.breeds, 0:(num.states - 1), sep = "_")

  dimnames(emis) = list(emis.breeds, colnames(geno))

  return(emis)

} # calc_emission_probs()


# Read in the phased VCF files and write out to *.rds files.
setwd(input.dir)

vcf.files = dir(path = input.dir, pattern = "phased.vcf.gz$", full.names = TRUE)
fam.files = dir(path = input.dir, pattern = "fam$", full.names = TRUE)
map.files = dir(path = input.dir, pattern = "map$", full.names = TRUE)

# Verify that they're the same length.
stopifnot(length(vcf.files) == length(fam.files))
stopifnot(length(vcf.files) == length(map.files))

for(i in 1:length(vcf.files)) {

  print(i)
  geno = read_geno(vcf.file = vcf.files[i], fam.file = fam.files[i],
         map.file = map.files[i])
  saveRDS(geno, file = sub("vcf.gz", "rds", vcf.files[i]))

} # for(i)


# Make objects to pass into calc_emission_probs().
setwd(input.dir)
geno.files = dir(pattern = "phased.rds$")
fam.files  = dir(pattern = "fam$")
chr = strsplit(geno.files, "chr")
chr = sapply(chr, "[", 2)
chr = strsplit(chr, split = "\\.")
chr = as.numeric(sapply(chr, "[", 1))

for(i in 1:length(geno.files)) {

  print(i)

  emis = calc_emission_probs(geno.file = geno.files[i], fam.file = fam.files[i], 
                  window = window)
  saveRDS(emis, file = paste0(output.dir, "emission_probs_chr", chr[i], ".rds"))

} # for(i)




