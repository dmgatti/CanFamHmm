################################################################################
# Download the individual GEO GSM files. The series matrix files are not always
# in a consistent format. Extract the allele calls from the GSM files and 
# convert them from AA/BB format to ACGT format.
# 
# Daniel Gatti
# dan.gatti@jax.org
# May 24, 2017
################################################################################
# Illumina Canine HD array is GPL17481.
# NOTE: I tried to create one loop that downloaded the data, processed it and then
#       deleted the GSM file. It kept disconnecting from GEO, so I broke it up.
options(stringsAsFactors = F)
library(GEOquery)
library(R.utils)

dest.dir = "/hpcdata/dgatti/Dog/ExternalData/GEO"
setwd(dest.dir)

gpl = getGEO("GPL17481", destdir = dest.dir)

# Get the sample IDs.
samples = Meta(gpl)$sample_id
series  = Meta(gpl)$series_id

# Download the files and zip them.
for(i in 1:length(samples)) {

  result = tryCatch({
    f = getGEOfile(GEO = samples[i], destdir = dest.dir, AnnotGPL = FALSE, amount = "full")
  }, warning = function(w) {
    print(w)
  }, error = function(e) {
    # Wait for a while with looping code and try again.
    while(nchar(f) == 0) {
      f = getGEOfile(GEO = samples[i], destdir = dest.dir, AnnotGPL = FALSE, amount = "full")
    }
  }, finally = {
    gzip(paste0(samples[i], ".soft"), overwrite = T, remove = T)
  })

} # for(s)
