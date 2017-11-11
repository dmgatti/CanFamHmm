################################################################################
# Prepare Shannon data. Add breeds, order markers by CanFam3.1 and save as
# RDS and PED/MAP files.
# 
# Daniel Gatti
# dan.gatti@jax.org
# May 24, 2017
################################################################################
# After downloading the Shannon data from Dryad, run PLINK to convert the 
# BED and BIM files to PED and MAP files.
# module load plink/1.9
# plink --bfile ShannonBoyko_All --recode --out --tab ShannonNew

options(stringsAsFactors = FALSE)

library(data.table)

setwd("/hpcdata/dgatti/Dog/")

# Read in the combined data PLINK file.
data = fread("ExternalData/Shannon_etal/ShannonBoyko_New.ped")
data = as.data.frame(data)
data = as.matrix(data)

# Read in the breed info and merge into data file.
breeds = read.delim("ExternalData/Shannon_etal/ShannonBoyko_All.nopheno.txt",
         header = F)
stopifnot(nrow(breeds) == nrow(data))
m = match(data[,2], breeds[,1])
data[,1] = gsub(" ", "_", breeds[m,2])

# Sort by breed.
data = data[order(data[,1]),]

# Read in the old marker map.
map = read.delim("ExternalData/Shannon_etal/ShannonBoyko_New.map",
              header = F, sep = "\t")
dimnames(map) = list(map[,2], c("chr", "marker", "cM", "bp"))

# Add the marker names to the data.
colnames(data) = c("family", "sample", "father", "mother", "sex", "phenotype",
                 map[,2])

saveRDS(data, file = "ExternalData/Shannon_etal/ShannonBoyko_New.rds")


######
# START HERE ONCE THE ABOVE CODE HAS BEEN RUN.

data = readRDS(file = "ExternalData/Shannon_etal/ShannonBoyko_New.rds")

# Split up the data.
hdr  = data[,1:6]
data = data[,-(1:6)]

# Read in the new CanFam3.1 marker locations.
annot = read.delim("/hpcdata/dgatti/Dog/Annotation/Canine230K_Consortium_20004694_A1_StrandReport_FDT.txt",
          sep = "", comment.char = "#")

sum(colnames(data) %in% annot$SNP_Name)

data  = data[,colnames(data) %in% annot$SNP_Name]
annot = annot[annot$SNP_Name %in% colnames(data),]

# Sort the markers by chromosome and add Y and M as 40 and 41.
tmp.chr = annot$Chr
names(tmp.chr) = annot$SNP_Name
tmp.chr[tmp.chr == "X"] = "39"
tmp.chr[grep("^Y", names(tmp.chr))] = "40"
tmp.chr[grep("^chrM", names(tmp.chr))] = "41"
annot$Chr = tmp.chr

annot = annot[order(as.numeric(annot$Chr), annot$Coord),]

stopifnot(annot$SNP_Name %in% colnames(data))

data = data[,annot$SNP_Name]

# Make a PLINK map file with the CanFam3.1 coordinates.
markers = annot[,c(5, 2, 1, 6)]
rownames(markers) = markers[,2]
map = map[markers[,2],]
stopifnot(map[,2] == markers[,2])
markers[,3] = map[,3]

# Save the data.
data = cbind(hdr, data)
output.dir = "/hpcdata/dgatti/Dog/ExternalData/Shannon_etal/CanFam3.1/"
saveRDS(data, file = paste0(output.dir, "Shannon_cf3.1.rds"))
saveRDS(markers, file = paste0(output.dir, "Shannon_cf3.1_map.rds"))

write.table(data, file = paste0(output.dir, "Shannon_cf3.1.ped"), quote = F, 
            row.names = F, col.names = F)
# 159,070 markers.
write.table(markers, file = paste0(output.dir, "Shannon_cf3.1.map"), quote = F, 
            row.names = F, col.names = F)


