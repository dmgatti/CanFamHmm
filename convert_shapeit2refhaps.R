################################################################################
# Convert the reference haplotype files, which are in HAP/SAMPLE  format, to 
# the SHAPEIT HAP/LEGEND/SAMPLE format required for reference haplotypes.
# See https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#reference and
#     https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#hapsample.
#
# Daniel Gatti
# dan.gatti@jax.org
# Aug 2, 2017
################################################################################
options(stringsAsFactors = F)

library(R.utils)

#ref.dir = "/projects/dgatti/Dog/ReferenceData/PureBreedsHQ/shapeit/"
#out.dir = "/projects/dgatti/Dog/ReferenceData/PureBreedsHQ/shapeit/ref/"


# Arguments:
# REF_DIR: path to the directory where the HAP/SAMPLE files are.
# OUT_DIR: path to the output directory.
args = commandArgs(trailingOnly = TRUE)

if(length(args) != 3) {
  stop(paste("3 arguments required: CHR, REF_DIR and OUT_DIR"))
}

chr     = as.numeric(args[1])
ref.dir = args[2]
out.dir = args[3]


  print(paste("CHR", chr))

  # Get the input HAP and SAMPLE filenames.
  hap.file = dir(path = ref.dir, pattern = paste0("chr", chr, "_phased.haps.gz$"),
             full.names = T)
  if(length(hap.file) != 1) {
    stop(paste("HAP file for CHR", chr, "not found in ref.dir", ref.dir))
  }

  sample.file = dir(path = ref.dir, pattern = paste0("chr", chr, "_phased.sample$"),
                full.names = T)
  if(length(sample.file) != 1) {
    stop(paste("SAMPLE file for CHR", chr, "not found in ref.dir", ref.dir))
  }

  # Read in the data.
  hap = read.delim(hap.file, header = F, sep = "")
  sample = read.delim(sample.file, header = T, sep = "")
  sample = sample[-1,]

  stopifnot(2 * nrow(sample) == ncol(hap) - 5)

  # SAMPLE file.
  print("Writing REFERENCE files...")
  new.sample = data.frame(sample = sample[,2], population = sample[,1],
               group = sample[,1], sex = sample$sex)
  write.table(new.sample, file = paste0(out.dir, "ref_chr", chr, ".sample"),
              sep = " ", row.names = F, quote = F)

  # LEGEND file.
  legend = data.frame(id = hap[,2], position = hap[,3], a0 = hap[,4], 
           a1 = hap[,5], chr = hap[,1])
  write.table(legend, file = paste0(out.dir, "ref_chr", chr, ".legend"),
              sep = " ", row.names = F, quote = F)

 # HAP file.
  hap = hap[,-(1:5)]
  new.hap.filename = paste0(out.dir, "ref_chr", chr, ".haps")
  write.table(hap, file = new.hap.filename, sep = " ", row.names = F,
              col.names = F, quote = F)
  gzip(new.hap.filename, overwrite = TRUE, remove = TRUE)

