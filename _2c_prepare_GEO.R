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

#####
# Convert data to ACGT format.
options(stringsAsFactors = F)

setwd("/data/dgatti/Dog/")

geo.dir = "/hpcdata/dgatti/Dog/ExternalData/GEO/"
annot.dir = "Annotation/"

# Read in the annotation data.
annot = read.delim(paste0(annot.dir,
        "Canine230K_Consortium_20004694_A1_StrandReport_FDT.txt"),
        skip = 5)

# Oddly, there are no CT or GT SNPs.
table(paste0(annot$Top_AlleleA, annot$Top_AlleleB))

# Read in the GEO data.
geo = readRDS("/hpcdata/dgatti/Dog/ExternalData/GEO/GEO_geno.rds")

# Subset the annotation to match the GEO data.
annot = annot[order(annot$Chr, annot$Coord),]
annot = annot[annot$SNP_Name %in% rownames(geo),]
geo = geo[annot$SNP_Name,]

stopifnot(nrow(geo) == nrow(annot))

geo[geo == "NC"] = "00"
geo[geo == ""] = "00"
new.alleles = paste0(annot$Top_AlleleA, annot$Top_AlleleB)
for(i in 1:nrow(geo)) {

  geo[i,] = chartr("AB", new.alleles[i], geo[i,])

} # for(i)

saveRDS(geo, file = paste0(geo.dir, "GEO_geno_ACGT.rds"))

# Keep only the Series with enough dogs.
info = read.delim(paste0(geo.dir, "GEO_breed_sex.txt"))
series_to_keep = c("GSE52147", "GSE51734", "GSE51973", "GSE51976", "GSE52221", 
                   "GSE55134", "GSE66677", "GSE70454", "GSE80315", "GSE80735",
                   "GSE83151", "GSE83154", "GSE83160", "GSE90441", "GSE96736")
info = info[info$series %in% series_to_keep,]
geo = geo[,info$sample]

# Format in PLINK PED/MAP format.
map = read.delim("/hpcdata/dgatti/Dog/ExternalData/Shannon_etal/CanFam3.1/Shannon_cf3.1.map",
                 header = F, sep = "")
dim(map)

# NOTE: We're losing a lot of markers. There may have been marker name changes
#       between the Shannon data and the new array. 
# TBD: Try to retain more markers and impute them? The missing data is in the 
#      Shannon data set, which is larger than the GEO data.

# Subset to keep only the markers in Shannon.
geo = t(geo)
sum(map[,2] %in% colnames(geo))

map = map[map[,2] %in% colnames(geo),]
geo = geo[,colnames(geo) %in% map[,2]]
geo = geo[,map[,2]]
stopifnot(all(colnames(geo) == map[,2]))

# Write out the GEO map file (156198 rows).
write.table(map, file = "/hpcdata/dgatti/Dog/ExternalData/GEO/CanFam3.1/GEO_cf3.1.map",
            row.names = F, col.names = F, quote = F, sep = " ")

# Convert the GEO genotypes to space-separated values.
allele = unique(as.vector(geo))
repl = strsplit(allele, "")
repl = sapply(repl, paste, collapse = " ")
for(i in 1:length(allele)) {

  geo[geo == allele[i]] = repl[i]

} # for(i)

# Create a FAM header for the data.
fam = matrix("0", nrow(geo), 6)
# Family
fam[,1] = info$breed
# Sample
fam[,2] = rownames(geo)
# Sex
sex = factor(info$sex, levels = c("Male", "Female"))
# Lots of NAs in sex. Determine sex from X & Y Chr data.
chrx = which(map[,1] == "39")
chrx = chrx[rowMeans(apply(geo[,chrx], 1, "%in%", c("A C", "A G", "A T", "C G"))) < 0.25]
chry = which(map[,1] == "40")

hetx = colMeans(apply(geo[,chrx], 1, "%in%", c("A C", "A G", "A T", "C G")))
ncy  = rowMeans(geo[,chry] == "0 0")
sex.guess = rep("Female", nrow(info))
sex.guess[hetx < 0.05] = "Male"

na.rows = which(is.na(sex))
sex[na.rows] = sex.guess[na.rows]
sum(is.na(sex))
fam[,5] = sex

# Phenotype
fam[,6] = "-9"

geo = cbind(fam, geo)
write.table(geo, file = "/hpcdata/dgatti/Dog/ExternalData/GEO/CanFam3.1/GEO_cf3.1.ped",
            row.names = F, col.names = F, quote = F, sep = "\t")
write.table(map, file = "/hpcdata/dgatti/Dog/ExternalData/GEO/CanFam3.1/GEO_cf3.1.map",
            row.names = F, col.names = F, quote = F, sep = "\t")

saveRDS(geo, file = "/hpcdata/dgatti/Dog/ExternalData/GEO/CanFam3.1/GEO_cf3.1.rds")
saveRDS(map, file = "/hpcdata/dgatti/Dog/ExternalData/GEO/CanFam3.1/GEO_cf3.1_map.rds")






