##########
# Cluster each breed individually and look at the clustering to 
# evaluate the breed assignment.
options(stringsAsFactors = F)
library(BiocParallel)

setwd("/hpcdata/dgatti/Dog/")

ref.dir = "/hpcdata/dgatti/Dog/ReferenceData/"
fig.dir = "/hpcdata/dgatti/Dog/Figures/"
shannon.dir = "/hpcdata/dgatti/Dog/ExternalData/Shannon_etal/CanFam3.1/"

# Load in Shannon data.
#data    = readRDS(file = paste0(shannon.dir, "Shannon_cf3.1.rds"))
#markers = readRDS(file = paste0(shannon.dir, "Shannon_cf3.1_map.rds"))
#fam  = data[,1:6]
#data = as.matrix(data[,-(1:6)])

# Load in GEO data.
geo.dir = "/hpcdata/dgatti/Dog/ExternalData/GEO/CanFam3.1/"
#geo         = readRDS(file = paste0(geo.dir, "GEO_cf3.1.rds"))
#geo.markers = readRDS(file = paste0(geo.dir, "GEO_cf3.1_map.rds"))
#geo.fam  = geo[,1:6]
#geo = as.matrix(geo[,-(1:6)])

# Intersect the markers.
#marker.set = intersect(markers[,2], geo.markers[,2])

#print(paste("Retaining", length(marker.set), "markers."))

# Subset the data.
#keep = which(markers[,2] %in% marker.set)
#data    = data[,keep]
#markers = markers[keep,]

#keep = which(geo.markers[,2] %in% marker.set)
#geo         = geo[,keep]
#geo.markers = markers[keep,]

#stopifnot(ncol(geo) == ncol(data))
#stopifnot(ncol(geo.markers) == ncol(markers))

# Merge data sets.
#data = rbind(data, geo)
#fam  = rbind(fam, geo.fam)

#save(data, fam, markers, file = paste0(ref.dir, "ref_data_tmp.Rdata"))

#rm(geo, geo.markers, geo.fam)


load(file = paste0(ref.dir, "ref_data_tmp.Rdata"))

# Only keep samples with no-call rates less than
keep = which(rowMeans(data == "0 0") < 0.2)
data = data[keep,]
fam  = fam[keep,]

# Synch up alleles.
data[data == "C A"] = "A C"
data[data == "G A"] = "A G"
data[data == "T A"] = "A T"
data[data == "G C"] = "C G"

# Fix names that are different between data sets.
unique.breeds = sort(unique(fam[,1]))
fam[,1] = gsub(" ", "_", fam[,1])
fam[,1] = gsub("-", "_", fam[,1])
fam[,1] = gsub("'", "_", fam[,1])
fam[,1] = sub("Dogs$", "Dog", fam[,1])
fam[,1] = sub("dog$", "Dog", fam[,1])
fam[,1] = sub("_\\(AMAL\\)$", "", fam[,1])
fam[,1] = sub("collies$", "Collie", fam[,1])
fam[,1] = sub("Spaniels$", "Spaniel", fam[,1])
fam[,1] = sub("spaniel", "Spaniel", fam[,1])
fam[,1] = sub("^Brittany$", "Brittany_Spaniel", fam[,1])
fam[,1] = sub("^Chinese_Shar_pei$", "Shar_Pei", fam[,1])
fam[,1] = sub("Dogue_De_Bordeaux", "Dogue_de_Bordeaux", fam[,1])
fam[,1] = sub("bull", "Bull", fam[,1])
fam[,1] = sub("Eurasian", "Eurasier", fam[,1])
fam[,1] = sub("Flatcoated_Retriever", "Flat_Coated_Retriever", fam[,1])
fam[,1] = sub("sledge", "Sledge", fam[,1])
fam[,1] = sub("^New_Foundland$", "Newfoundland", fam[,1])
fam[,1] = sub("^Poodles$", "Poodle", fam[,1])
fam[,1] = sub("^Wolf$", "Gray_Wolf", fam[,1])
fam[grep("Village_Dog_DRC", fam[,1]),1] = "Village_Dog_DRC"
fam[grep("Village_Dog_Egypt", fam[,1]),1] = "Village_Dog_Egypt"
fam[grep("Village_Dog_Egypt", fam[,1]),1] = "Village_Dog_Egypt"
fam[grep("Village_Dog_Fiji", fam[,1]),1] = "Village_Dog_Fiji"
fam[grep("Village_Dog_French_Polynesia", fam[,1]),1] = "Village_Dog_French_Polynesia"
fam[grep("Village_Dog_India", fam[,1]),1] = "Village_Dog_India"
fam[grep("Village_Dog_Indonesia", fam[,1]),1] = "Village_Dog_Indonesia"
fam[grep("Village_Dog_Lebanon", fam[,1]),1] = "Village_Dog_Lebanon"
fam[grep("Village_Dog_Liberia", fam[,1]),1] = "Village_Dog_Liberia"
fam[grep("Village_Dog_Mexico", fam[,1]),1] = "Village_Dog_Mexico"
fam[grep("Village_Dog_Namibia", fam[,1]),1] = "Village_Dog_Namibia"
fam[grep("Village_Dog_Papua_New_Guinea", fam[,1]),1] = "Village_Dog_Papua_New_Guinea"
fam[grep("Village_Dog_Peru", fam[,1]),1] = "Village_Dog_Peru"
fam[grep("Village_Dog_Solomon_Islands", fam[,1]),1] = "Village_Dog_Solomon_Islands"
fam[grep("Village_Dog_Turkey", fam[,1]),1] = "Village_Dog_Turkey"
fam[grep("Village_Dog_Vietnam", fam[,1]),1] = "Village_Dog_Vietnam"
unique.breeds = sort(unique(fam[,1]))

pdf(paste0(fig.dir, "breed_cor.pdf"), width = 12, height = 8)
for(i in 1:length(unique.breeds)) {

  print(paste(i, "of", length(unique.breeds)))

  wh = which(fam[,1] == unique.breeds[i])
  if(length(wh) > 2) {
    tmp = data[wh,]
    tmp = data.frame(tmp, stringsAsFactors = T)
    tmp = lapply(tmp, as.numeric)
    tmp = matrix(unlist(tmp), nrow = length(tmp[[1]]), ncol = length(tmp),
          dimnames = list(paste(fam[wh,1], fam[wh,2], sep = "_"), NULL))
    cl = hclust(as.dist(1.0 - cor(t(tmp))), method = "average")
    plot(cl, main = unique.breeds[i])
  }

} # for(i)
dev.off()


