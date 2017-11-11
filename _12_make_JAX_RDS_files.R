################################################################################
# Read in the phased VCF files for the JAX data and write out to *.rds files.
# 
# Daniel Gatti
# dan.gatti@jax.org
# July 31, 2017
################################################################################
options(stringsAsFactors = FALSE)

library(VariantAnnotation)

input.dir = "/projects/dgatti/Dog/GeneseekData/Cleaned/"

setwd(input.dir)

# Get the phased VCF files.
vcf.files = dir(pattern = "phased.vcf.gz$")

# Read in each VCF file, transpose the genotypes (as 0 or 1) and write them
# out a *.rds file for faster I/O.
for(f in vcf.files) {

  print(f)

  vcf = readVcf(file = f, genome = "CanFam3.1")
  geno = geno(vcf)$GT
  dog.names = strsplit(colnames(geno), split = "_")
  dog.names = paste(unlist(dog.names), LETTERS[1:2], sep = "_")

  mat = matrix(0, length(dog.names), nrow(geno), dimnames = 
        list(dog.names, rownames(geno)))

  for(i in 1:ncol(geno)) {

    spl = strsplit(geno[,i], "\\|")
    spl = matrix(as.numeric(unlist(spl)), nrow = 2)
    mat[(2 * (i-1) + 1):(2*i),] = spl

  } # for(i)

  saveRDS(mat, file = sub("vcf.gz$", "rds", f))

} # for(f)


