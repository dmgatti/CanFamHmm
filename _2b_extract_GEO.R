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
source("/hpcdata/dgatti/Dog/Scripts/analysis_pipeline_v5/GEO_fxns.R")

dest.dir = "/hpcdata/dgatti/Dog/ExternalData/GEO"
setwd("/hpcdata/dgatti/Dog/")

gpl = getGEO("GPL17481", destdir = dest.dir)

# Get the sample IDs.
samples = Meta(gpl)$sample_id
series  = Meta(gpl)$series_id

# Read in the Illumina annotation.
annot = read.delim("Annotation/Canine230K_Consortium_20004694_A1_StrandReport_FDT.txt",
        skip = 5)

setwd(dest.dir)

# Read in the first file and get the number of markers.
files = dir(pattern = "soft.gz$")

# Go through each file and get the column names.
cn = as.list(files)
names(cn) = files
nr = as.list(files)
names(nr) = files

# Set up sample info data.frame.
info = data.frame(sample = samples, breed = rep(NA, length(samples)),
       sex = rep(NA, length(samples)), series  = rep(NA, length(samples)))

for(i in 1:length(files)) {

  x = scan(files[i], what = "char", sep = "\n", nlines = 100)
  
  nr[i] = x[grep("^!Sample_data_row_count", x)]

  if(length(nr[i]) > 0) {

    bounds = grep("sample_table", x)
    cn[[i]] = strsplit(x[bounds[1] + 1], "\t")[[1]]
    cn[[i]] = gsub("^([A-Z]+_)?[0-9]+\\.", "", cn[[i]])
    cn[[i]] = gsub(" ", "_", cn[[i]])
    cn[[i]] = toupper(cn[[i]])
  } # if(length(nr[i]) > 0)

  gr = grep("^!Sample_characteristics_ch1 = breed:", x)
  if(length(gr) > 0) {
    info$breed[i] = sub("^!Sample_characteristics_ch1 = breed: ", "", x[gr])
  }

  gr = grep("^!Sample_characteristics_ch1 = gender:", x)
  if(length(gr) > 0) {
    info$sex[i] = sub("^!Sample_characteristics_ch1 = gender: ", "", x[gr])
  }

  gr = grep("^!Sample_series_id = ", x)
  if(length(gr) > 0) {
    info$series[i] = sub("^!Sample_series_id = ", "", x[gr])
  }

} # for(i)
nr = as.numeric(sub("^!Sample_data_row_count = ", "", nr))

info$breed = gsub(" ", "_", info$breed)
info$breed = gsub("'", "_", info$breed)
info$breed[info$breed == "Alaskan_Malamute (AMAL)"] = "Alaskan_Malamute"
info$breed[info$breed == "Bernese_Mountain_Dogs"]   = "Bernese_Mountain_Dog"
info$breed[info$breed == "Border_collies"]          = "Border_Collie"
info$breed[info$breed == "Brittany_Spaniels"]       = "Brittany_Spaniel"
info$breed[info$breed == "Cocker_spaniel"]          = "Cocker_Spaniel"
info$breed[info$breed == "English_bulldog"]         = "English_Bulldog"
info$breed[info$breed == "Flatcoated_Retriever"]    = "Flat_Coated_Retriever"
info$breed[info$breed == "Greenland_sledge_dog"]    = "Greenland_Sledge_Dog"
info$breed[info$breed == "New_Foundland"]           = "Newfoundland"
info$breed[info$breed == "Pembroke_Welsh_corgi_(PWC)"] = "Pembroke_Welsh_Corgi"
info$breed[info$breed == "Miniature_Poodle"]        = "Toy_Poodle"
info$breed[info$breed == "Poodles"]                 = "Standard_Poodle"
info$breed[info$breed == "Portuguese_Water_dog"]    = "Portuguese_Water_Dog"
info$breed[info$breed == "Siberian_Husky_(HUSK)"]   = "Siberian_Husky"
info$breed[info$breed == "greyhound"]               = "Greyhound"
info$breed[info$breed == "rottweiler"]              = "Rottweiler"

info$sex[info$sex == "female"] = "Female"
info$sex[info$sex == "male"]   = "Male"

# NOTE: At this point, there are some sex values missing. We fill these in
# below once we have the genotype data.
write.table(info, file= "GEO_breed_sex.txt", sep = "\t", quote = F, row.names = F)

# There are three types of column names. Note that we have to substitute
# a lot of text to get down to the three classes. When the number of rows = 0,
# there are no column names.
unique.colnames = unique(cn)
unique.colnames = unique.colnames[!is.na(unique.colnames)]

# Read in the first file and set up the storage matrices.
x = scan(files[1], what = "char", sep = "\n")
nr = x[grep("^!Sample_data_row_count", x)]
nr = strsplit(nr, " = ")[[1]]
nmarkers = as.numeric(nr[2])

bounds = grep("sample_table", x)
data = x[(bounds[1] + 2):(bounds[2] - 1)]
data = strsplit(data, "\t")
data = matrix(unlist(data), nrow = nmarkers, byrow = T)
colnames(data) = strsplit(x[bounds[1] + 1], "\t")[[1]]

geno = matrix("", nrow = nmarkers, ncol = length(samples), dimnames = list(data[,"ID_REF"], 
       samples))
gcscore = matrix("", nrow = nrow(geno), ncol = ncol(geno), dimnames = dimnames(geno))
theta = matrix("", nrow = nrow(geno), ncol = ncol(geno), dimnames = dimnames(geno))
rho = matrix("", nrow = nrow(geno), ncol = ncol(geno), dimnames = dimnames(geno))
baf = matrix("", nrow = nrow(geno), ncol = ncol(geno), dimnames = dimnames(geno))
lrr = matrix("", nrow = nrow(geno), ncol = ncol(geno), dimnames = dimnames(geno))

tmp.cn = sapply(cn, paste, collapse = "")
cn = split(cn, tmp.cn)

# Organize the files into different groups. Read in the ones in the same
# format together.
# Column IDs: ID_REF, VALUE, GC_SCORE, THETA, R, B_ALLELE_FREQ, LOG_R_RATIO
# 173662 rows
for(i in 1:length(cn[[1]])) {

  print(paste(i, "of", length(cn[[1]])))

  x = scan(names(cn[[1]])[i], what = "char", sep = "\n")
  data.column = which(colnames(geno) == sub("\\.soft\\.gz$", "", names(cn[[1]])[i]))

  # Get the data from the file.
  bounds = grep("sample_table", x)
  data = x[(bounds[1] + 2):(bounds[2] - 1)]
  data = strsplit(data, "\t")
  data = matrix(unlist(data), nrow = length(data), byrow = T)
  colnames(data) = strsplit(x[bounds[1] + 1], "\t")[[1]]
  colnames(data) = toupper(colnames(data))

  geno[,data.column]    = data[,"VALUE"]
  gcscore[,data.column] = data[,"GC_SCORE"]
  theta[,data.column]   = data[,"THETA"]
  rho[,data.column]     = data[,"R"]
  baf[,data.column]     = data[,"B_ALLELE_FREQ"]
  lrr[,data.column]     = data[,"LOG_R_RATIO"]

} # for(i)

# Column IDs: ID_REF, VALUE, SCORE, THETA, R, B_ALLELE_FREQ, LOG_R_RATIO
# 173662 rows
for(i in 1:length(cn[[2]])) {

  print(paste(i, "of", length(cn[[2]])))

  x = scan(names(cn[[2]])[i], what = "char", sep = "\n")
  data.column = which(colnames(geno) == sub("\\.soft\\.gz$", "", names(cn[[2]])[i]))

  # Get the data from the file.
  bounds = grep("sample_table", x)
  data = x[(bounds[1] + 2):(bounds[2] - 1)]
  data = strsplit(data, "\t")
  data = matrix(unlist(data), nrow = length(data), byrow = T)
  colnames(data) = strsplit(x[bounds[1] + 1], "\t")[[1]]
  colnames(data) = gsub("^([A-Z]+_)?[0-9]+\\.", "", colnames(data))
  colnames(data) = gsub(" ", "_", colnames(data))
  colnames(data) = toupper(colnames(data))

  geno[,data.column]    = data[,"VALUE"]
  gcscore[,data.column] = data[,"SCORE"]
  theta[,data.column]   = data[,"THETA"]
  rho[,data.column]     = data[,"R"]
  baf[,data.column]     = data[,"B_ALLELE_FREQ"]
  lrr[,data.column]     = data[,"LOG_R_RATIO"]

  gr = grep("^!Sample_characteristics_ch1 = breed:", x)
  if(length(gr) > 0) {
    info$breed[data.column] = sub("^!Sample_characteristics_ch1 = breed: ", "", x[gr])
  }

  gr = grep("^!Sample_characteristics_ch1 = gender:", x)
  if(length(gr) > 0) {
    info$sex[data.column] = sub("^!Sample_characteristics_ch1 = gender: ", "", x[gr])
  }

  gr = grep("^!Sample_series_id = ", x)
  if(length(gr) > 0) {
    info$series[data.column] = sub("^!Sample_series_id = ", "", x[gr])
  }

} # for(i)


# Column IDs: ID_REF, VALUE, SCORE, THETA, R, GTYPE, B_ALLELE_FREQ, LOG_R_RATIO
# 172115 rows.
for(i in 1:length(cn[[3]])) {

  print(paste(i, "of", length(cn[[3]])))

  x = scan(names(cn[[3]])[i], what = "char", sep = "\n")
  data.column = which(colnames(geno) == sub("\\.soft\\.gz$", "", names(cn[[3]])[i]))

  # Get the data from the file.
  bounds = grep("sample_table", x)
  data = x[(bounds[1] + 2):(bounds[2] - 1)]
  data = strsplit(data, "\t")

  # For some reason, some markres have only 7 values.
  wh = which(sapply(data, length) == 7)
  data[wh] = lapply(data[wh], "c", "")

  data = matrix(unlist(data), nrow = length(data), byrow = T)
  colnames(data) = strsplit(x[bounds[1] + 1], "\t")[[1]]
  colnames(data) = gsub("^([A-Z]+_)?[0-9]+\\.", "", colnames(data))
  colnames(data) = gsub(" ", "_", colnames(data))
  colnames(data) = toupper(colnames(data))

  m = match(data[,"ID_REF"], rownames(geno))

  geno[m, data.column]    = data[,"GTYPE"]
  gcscore[m, data.column] = data[,"SCORE"]
  theta[m, data.column]   = data[,"THETA"]
  rho[m, data.column]     = data[,"R"]
  baf[m, data.column]     = data[,"B_ALLELE_FREQ"]
  lrr[m, data.column]     = data[,"LOG_R_RATIO"]

  gr = grep("^!Sample_characteristics_ch1 = breed:", x)
  if(length(gr) > 0) {
    info$breed[data.column] = sub("^!Sample_characteristics_ch1 = breed: ", "", x[gr])
  }

  gr = grep("^!Sample_characteristics_ch1 = gender:", x)
  if(length(gr) > 0) {
    info$sex[data.column] = sub("^!Sample_characteristics_ch1 = gender: ", "", x[gr])
  }

  gr = grep("^!Sample_series_id = ", x)
  if(length(gr) > 0) {
    info$series[data.column] = sub("^!Sample_series_id = ", "", x[gr])
  }

} # for(i)

write.table(info, file = "GEO_breed_sex.txt", quote = F, sep = "\t")
saveRDS(geno, file = "GEO_geno.rds")
saveRDS(gcscore, file = "GEO_gc_score.rds")
saveRDS(theta, file = "GEO_theta.rds")
saveRDS(rho, file = "GEO_rho.rds")
saveRDS(baf, file = "GEO_baf.rds")
saveRDS(lrr, file = "GEO_lrr.rds")

# GSE51976 has some numbers in the genotype data.
# Get the genotype data.
geno = readRDS(file = "GEO_geno.rds")

# Get the Series Matrix file.
gse = getGEO("GSE51976", destdir = dest.dir)[[1]]
gse = gse[rownames(geno),]
stopifnot(rownames(geno) == rownames(gse))
rng = match(colnames(gse), colnames(geno))
geno[,rng] = exprs(gse)
rm(gse)

# Dreger: GSE 90441 has the data in a separate file.
# Add the genotypes to the matrix.
gse.dir = "/hpcdata/dgatti/Dog/ExternalData/Dreger_etal/GSE90441/"
data = readRDS(paste0(gse.dir, "GSE90441_geno.rds"))
colnames(data) = sub("\\.GType", "", colnames(data))

# Read in breed infofor GSE90441.
info = read.csv(paste0(gse.dir, "GSE90441_breed_codes.csv"))
stopifnot(all(colnames(data) == info$number))

# Subset the markers to match.
data = data[rownames(geno),]
all(rownames(geno) == rownames(data))

# Replace the NA data with the GSE90441 genotypes.
rng = match(info$sample, colnames(geno))
geno[,rng] = data

saveRDS(geno, file = "GEO_geno.rds")

# Set the sex in the breed info file.
# Read in the annotation data.
annot = read.delim("../../Annotation/Canine230K_Consortium_20004694_A1_StrandReport_FDT.txt",
        skip = 5)
annot = annot[match(rownames(geno), annot$SNP_Name),]
stopifnot(nrow(annot) == nrow(geno))

info = read.delim("GEO_breed_sex.txt")

chrx = which(annot$Chr == "X")
chrx = chrx[which(rowMeans(geno[chrx,] == "NC") < 0.05)]
chrx = chrx[rowMeans(geno[chrx,] == "AB") < 0.25]
chry = grep("^Y", annot$SNP_Name)

sex = rep("F", nrow(geno))

hetx = colSums(geno[chrx,] == "AB")
ncy = colSums(geno[chry,] == "NC")
plot(hetx, ncy)
plot(density(hetx, adjust = 0.5))
abline(v = 220)

sex[hetx < 220] = "M"
cbind(info$sex, sex)
info$sex = sex
write.table(info, file = "GEO_breed_sex.txt", quote = F, sep = "\t")


##################################
# Remove the downloaded GSM files.
#setwd(dest.dir)
#files = dir(pattern = "soft.gz$")

#for(f in files) {

#  file.remove(f)

#} # for(f)


