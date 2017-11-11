################################################################################
# Read in the JAX data and write out a *.rds file. 
################################################################################
options(stringsAsFactors = F)

# Read in the work functions.
source("/projects/dgatti/Dog/Scripts/CanFamHmm/jax_input_fxns.R")

#### STEP 1: Extract the JAX data from the Geneseek files.

# JAX1
# NOTE: Geneseek renamed the samples, so we have to read in the data
#       and rename the samples.
canFam3 = read.delim("/projects/dgatti/Dog/CanFam3/canineHD.cf3.1.map",
          header = F)
extract_jax_data(input.dir = "GeneseekData/JAX1/", canFam3 = canFam3)

# Read in the sample disambiguation file. (Geneseek assigned new numbers
# to each sample.
data = readRDS("GeneseekData/JAX1/jax_cleaned_data.rds")

samples = read.delim("GeneseekData/JAX1/sample_disambiguation.txt")
samples[,2] = make.names(samples[,2])
samples = samples[match(rownames(data), samples[,1]),]
stopifnot(rownames(data) == samples[,1])

# Replace sample names.
rownames(data) = samples[,2]
data[,1] = samples[,2]
data[,2] = samples[,2]
saveRDS(data, "GeneseekData/JAX1/jax_cleaned_data.rds")


# JAX2
extract_jax_data(input.dir = "GeneseekData/JAX2/", canFam3 = canFam3)


# JAX3
extract_jax_data(input.dir = "GeneseekData/JAX3/", canFam3 = canFam3)


# JAX4
extract_jax_data(input.dir = "GeneseekData/JAX4/", canFam3 = canFam3)

# Change some sample IDs.
data = readRDS("GeneseekData/JAX4/jax_cleaned_data.rds")

rownames(data)[rownames(data) == "Ben"] = "Ben2"
data[,1:2] = gsub("^Ben$", "Ben2", data[,1:2])
rownames(data)[rownames(data) == "Olaf"] = "Olaf2"
data[,1:2] = gsub("^Olaf$", "Olaf2", data[,1:2])
rownames(data)[rownames(data) == "Watson"] = "Watson2"
data[,1:2] = gsub("^Watson$", "Watson2", data[,1:2])
rownames(data)[rownames(data) == "Brady-Ouellette"] = "Brady2"
data[,1:2] = gsub("^Brady-Ouellette$", "Brady2", data[,1:2])

saveRDS(data, file = "GeneseekData/JAX4/jax_cleaned_data.rds")

# Combine the data.
setwd("/projects/dgatti/Dog/GeneseekData/")
jax1 = readRDS("JAX1/jax_cleaned_data.rds")
jax2 = readRDS("JAX2/jax_cleaned_data.rds")
jax3 = readRDS("JAX3/jax_cleaned_data.rds")
jax4 = readRDS("JAX4/jax_cleaned_data.rds")
jax.markers = readRDS("JAX1/jax_cleaned_markers.rds")

stopifnot(ncol(jax1) == ncol(jax2))
stopifnot(ncol(jax1) == ncol(jax3))
stopifnot(ncol(jax1) == ncol(jax4))

rownames(jax3) = make.names(rownames(jax3))
rownames(jax4) = make.names(rownames(jax4))

jax = rbind(jax1, jax2, jax3, jax4)

saveRDS(jax, file = "Cleaned/jax_cleaned.rds")
saveRDS(jax.markers, file = "Cleaned/jax_markers.rds")

