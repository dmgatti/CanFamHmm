################################################################################
# Merge the Shannon and GEO data, cluster by breed and remove samples that
# don't cluster with their breed.
# Daniel Gatti
# dan.gatti@jax.org
# June 3, 2017
################################################################################
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

# Only keep samples with no-call rates less than 20% 
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

# Compute kinship and cluster the data.
cor_fxn = function(index) {

  colMeans(data[,index] == data)

} # cor_fxn()

data = t(data)

param = MulticoreParam(workers = 10, progressbar = TRUE, log = TRUE)

kin = bplapply(X = 1:ncol(data), FUN = cor_fxn, BPPARAM = param)

K = matrix(0, ncol(data), ncol(data), dimnames = list(colnames(data), 
    colnames(data)))
for(i in 1:length(kin)) {
  K[,i] = kin[[i]]
}

cl = hclust(as.dist(1.0 - K), method = "average")
saveRDS(cl, file = "global_cluster.rds")

png(paste0(fig.dir, "global_kinship_clust.png"), width = 25000, height = 2000)
plot(cl, cex.axis = 0.7)
dev.off()

data = data[,cl$order]
fam = fam[cl$order,]
K = K[cl$order, cl$order]

dimnames(K) = list(fam[,1], fam[,1])

saveRDS(cbind(fam, t(data)), file = paste0(ref.dir, "ref_genotypes.rds"))
saveRDS(markers, file = paste0(ref.dir, "ref_markers.rds"))

data    = readRDS(file = paste0(ref.dir, "ref_genotypes.rds"))
markers = readRDS(file = paste0(ref.dir, "ref_markers.rds"))

saveRDS(K, file = paste0(ref.dir, "ref_kinship.rds"))
rm(kin)

K = readRDS(file = paste0(ref.dir, "ref_kinship.rds"))

png(paste0(fig.dir, "global_kinship.png"), width = 20000, height = 20000)
image(1:nrow(K), 1:ncol(K), K, axes = F, ann = F)
at = seq(1, nrow(K), 2)
mtext(side = 1, line = 0.2, at = at, text = rownames(K)[at])
mtext(side = 2, line = 0.2, at = at, text = rownames(K)[at], las = 1)
mtext(side = 3, line = 0.2, at = at, text = rownames(K)[at])
mtext(side = 4, line = 0.2, at = at, text = rownames(K)[at], las = 1)
dev.off()


pdf(paste0(fig.dir, "global_kinship_breed_order.pdf"))

  unique.breeds = sort(unique(data[,1]))
  for(br in unique.breeds) {
    wh = which(data[,1] == br)
    plot(wh, pch = 16, main = br, col = 0)
    text(1:length(wh), wh, labels = data[wh,1])
  } # for(i)

dev.off()


# Remove or change strain names based on clustering.
# Make a vector of samples that were changed.
changed = NULL

# Airedale Terrier
gr = which(data[,1] == "Airedale_Terrier")
# "PFZ35H02" is a Staffordshire Terrier
wh = which(data[,2] == "PFZ35H02")
new.name = "American_Staffordshire_Terrier"
changed = rbind(changed, c(data[wh,1:2], new.name, "Change"))

# American Cocker Spaniel
# Remove two samples that don't cluster: GSM1330057, GSM1330056.
gr = which(data[,1] == "Cocker_Spaniel")
wh = gr[gr< 6000]
changed = rbind(changed, cbind(data[wh,1:2], rep(NA, length(wh)),
          rep("Remove", length(wh))))
# Change the "Cocker_Spaniel" above index 6530 to "American_Cocker_Spaniel".
wh = which(data[,1] == "Cocker_Spaniel")
new.name = "American_Cocker_Spaniel"
changed = rbind(changed, cbind(data[wh,1:2], rep(new.name, length(wh)),
          rep("Change", length(wh))))

# American Pit Bull Terrier.
# American Staffordshire Terrier is the American breed.
# Staffordshire Bull Terrier is the English breed.
# In the data, American Pit Bull Terrier, American Staffordshire Bull Terrier,
# American Staffordshire Terrier and Staffordshire Bull Terrier all
# appear. There is a cluster near 440 and a larger cluster near 1640.
# Remove the ones near 440 because they cluster with Glen of Imaal Terriers.
gr = grep("Staffordshire|Pit_Bull", data[,1])
wh = gr[gr < 500]
changed = rbind(changed, cbind(data[wh,1:2], rep(NA, length(wh)),
          rep("Remove", length(wh))))
# Rename the rest "American_Staffordshire_Terrier"
wh = grep("Staffordshire|Pit_Bull", data[,1])
new.name = "American_Staffordshire_Terrier"
changed = rbind(changed, cbind(data[wh,1:2], rep(new.name, length(wh)),
          rep("Change", length(wh))))

# Bearded Collie.
gr = grep("Bearded_Collie", data[,1])
# Remove sample "PFZ35E05" because it clusters with other samples.
wh = gr[gr > 4000]
changed = rbind(changed, c(data[wh,1:2], NA, "Remove"))

# Belgian Tervuren
gr = grep("Belgian_Tervuren", data[,1])
# One of teh Belgian Tervurens is a Belgian Sheep Dog.
wh = gr[gr < 2075]
new.name = "Belgian_SheepDog"
changed = rbind(changed, c(data[wh,1:2], new.name, "Change"))

# Boxer
# There are 2 clusters of Boxers, one near 120 and one near 7540.
# The cluster near 7540 is part of a large cluster that is separate from
# all other samples.
gr = which(data[,1] == "Boxer")
# Remove samples > 7000.
wh = gr[gr > 7000]
changed = rbind(changed, cbind(data[wh,1:2], rep(NA, length(wh)),
          rep("Remove", length(wh))))

# BullDog
gr = which(data[,1] == "BullDog")
# All of the BullDogs cluster with English_BullDogs.
wh = gr
new.name = "English_BullDog"
changed = rbind(changed, cbind(data[wh,1:2], rep(new.name, length(wh)),
          rep("Change", length(wh))))

# Bull Terrier
gr = which(data[,1] == "Bull_Terrier")
# The sample near 1600 clusters with American Staffordshire Terrier.
wh = gr[gr > 1600]
new.name = "American_Staffordshire_Terrier"
changed = rbind(changed, c(data[wh,1:2], new.name, "Change"))

# Cairn_Terrier
# There are 2 groups of Cairn Terriers, one near 250 and one head 2975.
# The samples near 250 cluster with Yorkshire_Terrier and the samples near 
# 2975 cluster near Scottish_Terrier and West_Highland_White_Terrier.
# The two groups split between one with dog names and others. Maybe it's
# two different sets of Cairn Terriers.
# Rename the samples near 250 as Cairn_Terrier_1 and the samples near 2975
# Cairn_Terrier_2.
gr = which(data[,1] == "Cairn_Terrier")
wh = gr[which(gr < 2000)]
new.name = "Cairn_Terrier_1"
changed = rbind(changed, cbind(data[wh,1:2], rep(new.name, length(wh)),
          rep("Change", length(wh))))
wh = gr[which(gr > 2000)]
new.name = "Cairn_Terrier_2"
changed = rbind(changed, cbind(data[wh,1:2], rep(new.name, length(wh)),
          rep("Change", length(wh))))

# Canaan_Dog
# One dog doesn't cluser with the others. Remove "Janey".
gr = which(data[,1] == "Canaan_Dog")
wh = gr[gr < 1050]
changed = rbind(changed, c(data[wh,1:2], NA, "Remove"))

# Cane Corso
# There are two groups that cluster separately, but near Neopolitan Mastiff.
# Leaving them unchanged because they cluster near a related breed.

# Carolina Dog.
# On Carolina Dog clusters with Huskies and Alaska Dogs.
# Remove PFZ34D08.
gr = which(data[,1] == "Carolina_Dog")
wh = gr[gr < 2000]
changed = rbind(changed, c(data[wh,1:2], NA, "Remove"))

# Catahoula_Leopard_Dog
# Sample Li86 clusters with Retrievers. 
gr = which(data[,1] == "Catahoula_Leopard_Dog")
# Remove Li86.
wh = gr[gr > 5500]
changed = rbind(changed, c(data[wh,1:2], NA, "Remove"))

# Chihuahua
# There are two groups of Chihuahuas, one near 4770 and one near 7160.
# They come from different data sets and it's not clear which one is correct.
# The cluster near 4770 clusters with Chinese Crested.
# The cluster near 7160 clusters with Briards and Kuvasz's. 
# Remove the cluster near 7160.
gr = which(data[,1] == "Chihuahua")
wh = gr[gr > 7000]
changed = rbind(changed, cbind(data[wh,1:2], rep(NA, length(wh)),
          rep("Remove", length(wh))))

# Chinese_Shar_pei
gr = which(data[,1] == "Chinese_Shar_pei")
# Remove PFZ9H06 because it clusters with the Chows.
wh = gr[gr > 1350]
changed = rbind(changed, c(data[wh,1:2], NA, "Remove"))

# English_BullDog
gr = which(data[,1] == "English_BullDog")
# Remove PFZ11A12 because it clusters with Village Dogs.
wh = gr[gr < 1650]
changed = rbind(changed, c(data[wh,1:2], NA, "Remove"))

# English_Cocker_Spaniel
gr = which(data[,1] == "English_Cocker_Spaniel")
# Two samples, "CoS_GT103" and "CoS_GT104" cluster with Spaniels.
# Remove "CoS_GT103" and "CoS_GT104".
wh = gr[gr < 2000]
changed = rbind(changed, cbind(data[wh,1:2], rep(NA, 2), rep("Remove", 2)))

# English_Springer_Spaniel
gr = which(data[,1] == "English_Springer_Spaniel")
# Two samples, GSM2400856 and PFZ43C03 cluster with other breeds.
# Remove GSM2400856 and PFZ43C03.
wh = gr[gr < 2000]
changed = rbind(changed, cbind(data[wh,1:2], rep(NA, length(wh)),
          rep("Remove", length(wh))))

# Flat_Coated_Retriever
gr = which(data[,1] == "Flat_Coated_Retriever")
# Remove PFZ15B05 because it clusters with other Retrievers.
wh = gr[gr > 5400]
changed = rbind(changed, c(data[wh,1:2], NA, "Remove"))

# German_Shorthaired_Pointer
gr = which(data[,1] == "German_Shorthaired_Pointer")
# Remove PFZ10D10 because it clusters with other Retrievers.
wh = gr[gr > 4000]
changed = rbind(changed, c(data[wh,1:2], NA, "Remove"))

# Golden_Retriever
gr = which(data[,1] == "Golden_Retriever")
# Remove GSM1330188 because it clusters with Rhodesian_Ridgebacks.
wh = gr[gr < 4000]
changed = rbind(changed, c(data[wh,1:2], NA, "Remove"))

# Gordon_Setter
# Gordon_Setters have two clusters, one near 795 and one near 6825.
# The cluster near 795 clusters with Spitz dogs.
# The cluster near 6825 clusters near other Setters.
gr = which(data[,1] == "Gordon_Setter")
# Remove dogs in the 795 cluster.
wh = gr[gr < 6000]
changed = rbind(changed, cbind(data[wh,1:2], rep(NA, length(wh)), 
          rep("Remove", length(wh))))

# Greyhound
gr = which(data[,1] == "Greyhound")
# There are two clusers of Greyhouds, one near 3640 and one near 7193.
# There are 278 dogs in the cluster at 7193. But it's part of an outlier
# cluster. There are only 43 dogs in teh cluster at 3640, but it clusters 
# with Whippets.
# All of the dogs near 7193 are from GSE52147.
# Remove dogs in the 7193 cluster.
wh = gr[gr > 7000]
changed = rbind(changed, cbind(data[wh,1:2], rep(NA, length(wh)), 
          rep("Remove", length(wh))))

# Havanese
gr = which(data[,1] == "Havanese")
# Remove PFZ13C03 because it clusters with Bichon Frises.
wh = gr[gr > 4700]
changed = rbind(changed, c(data[wh,1:2], NA, "Remove"))

# Jack_Russell_Terrier
gr = which(data[,1] == "Jack_Russell_Terrier")
# Remove PFZ16C10 because it clusters with South American Village Dogs.
wh = gr[gr > 4000]
changed = rbind(changed, c(data[wh,1:2], NA, "Remove"))

# Labrador_Retriever
gr = which(data[,1] == "Labrador_Retriever")
# There are two Labrador clusters, one near 544 and one near 5908.
# There are also a few outliers in between.
# Remove the samples in between.
wh = gr[gr > 1000 & gr < 5300]
changed = rbind(changed, cbind(data[wh,1:2], rep(NA, length(wh)), 
          rep("Remove", length(wh))))
# Rename the each cluster as 1 and 2. 
wh = gr[which(gr < 1000)]
new.name = "Labrador_Retriever_1"
changed = rbind(changed, cbind(data[wh,1:2], rep(new.name, length(wh)), 
          rep("Change", length(wh))))
wh = gr[which(gr > 5200)]
new.name = "Labrador_Retriever_2"
changed = rbind(changed, cbind(data[wh,1:2], rep(new.name, length(wh)), 
          rep("Change", length(wh))))
 
# Maltese.
gr = which(data[,1] == "Maltese")
# Remove two samples that don't cluster with Maltese: PFZ29F05 & PFZ37C11.
wh = gr[gr < 4600 | gr > 4750]
changed = rbind(changed, cbind(data[wh,1:2], rep(NA, length(wh)), 
          rep("Remove", length(wh))))

# Pembroke_Welsh_Corgi
gr = which(data[,1] == "Pembroke_Welsh_Corgi")
# There are two clusters of Corgis, one near 2500 and one near 7520.
# The cluster near 2500 cluster with the Cardigan Welsh Corgis.
# The cluster near 7520 is part of a large outlier cluster with Greyhounds and 
# from GSE80735.
# Remove samples in the cluster near 7520.
wh = gr[gr > 7000]
changed = rbind(changed, cbind(data[wh,1:2], rep(NA, length(wh)), 
          rep("Remove", length(wh))))

# Peruvian_Inca_Orchid
gr = which(data[,1] == "Peruvian_Inca_Orchid")
# Remove two samples that don't cluster with the rest: PFZ2F05 & PFZ2C05.
wh = gr[gr < 3000 | gr > 6000]
changed = rbind(changed, cbind(data[wh,1:2], rep(NA, length(wh)), 
          rep("Remove", length(wh))))

# Poodle
gr = which(data[,1] == "Poodle")
# There are 2 clusters of Poodles, one near 680 and one near 7095.
# The cluster near 680 clusters with Standard_Poodles.
# The cluster near 7095 clusters with Toy_Poodle.
# Change the Poodle in the low cluster to Standard_Poodle.
wh = gr[gr < 1000]
new.name = "Standard_Poodle"
changed = rbind(changed, cbind(data[wh,1:2], rep(new.name, length(wh)), 
          rep("Change", length(wh))))
# Change the Poodle in the high cluster to Toy_Poodle.
wh = gr[gr > 7000]
new.name = "Toy_Poodle"
changed = rbind(changed, cbind(data[wh,1:2], rep(new.name, length(wh)), 
          rep("Change", length(wh))))

# Rat_Terrier
gr = which(data[,1] == "Rat_Terrier")
# It's not clear if there is one breed here, but two samples are far from 
# the others.
# Remove two samples: PFZ35A11 & PFZ33E05
wh = gr[gr > 4000]
changed = rbind(changed, cbind(data[wh,1:2], rep(NA, length(wh)), 
          rep("Remove", length(wh))))

# Schipperke
gr = which(data[,1] == "Schipperke")
# One sample doesn't cluster with the rest.
# Remove one sample: "PFZ43E07
wh = gr[gr > 6000]
changed = rbind(changed, c(data[wh,1:2], NA, "Remove"))

# Siberian_Husky
gr = which(data[,1] == "Siberian_Husky")
# Remove two huskies that don't cluster with the others: 
wh = gr[gr > 1500]
changed = rbind(changed, cbind(data[wh,1:2], rep(NA, length(wh)), 
          rep("Remove", length(wh))))

# Sloughi
gr = which(data[,1] == "Sloughi")
# It looks like one Sloughi is a Saluki: GSM2195606
# Change GSM2195606 to Saluki.
wh = gr[gr < 1000]
new.name = "Saluki"
changed = rbind(changed, c(data[wh,1:2], new.name, "Change"))

# Standard_Schnauzer
gr = which(data[,1] == "Standard_Schnauzer")
# Change 2 Standard Schnauzer because they cluster with Mineature Schnauzers:
# PFZ9A09 & PFZ18A01
wh = gr[gr > 6000]
new.name = "Miniature_Schnauzer"
changed = rbind(changed, cbind(data[wh,1:2], rep(new.name, length(wh)), 
          rep("Change", length(wh))))

# Terrier_Yorkshire
gr = which(data[,1] == "Terrier_Yorkshire")
# These all cluster together near Yorkshire_Terrier.
# Change these to "Yorkshire_Terrier_2".
wh = gr
new.name = "Yorkshire_Terrier_2"
changed = rbind(changed, cbind(data[wh,1:2], rep(new.name, length(wh)), 
          rep("Change", length(wh))))

# Tibetan_Mastiff
gr = which(data[,1] == "Tibetan_Mastiff")
# There are two clusters of Tibetan Mastiffs: One near 1290 and one near 1528.
# Each one clusters with a different east-Asian dogs.
# Make 2 clusters called Tibetan_Mastiff_1 and Tibetan_Mastiff_2.
wh = gr[gr < 1400]
new.name = "Tibetan_Mastiff_1"
changed = rbind(changed, cbind(data[wh,1:2], rep(new.name, length(wh)), 
          rep("Change", length(wh))))
wh = gr[gr > 1400]
new.name = "Tibetan_Mastiff_2"
changed = rbind(changed, cbind(data[wh,1:2], rep(new.name, length(wh)), 
          rep("Change", length(wh))))

# Village_Dog_Afghanistan
gr = which(data[,1] == "Village_Dog_Afghanistan")
# Remove 2 dogs that don't cluster with the others: PFZ25B01 & 8317
wh = gr[gr < 1050]
changed = rbind(changed, cbind(data[wh,1:2], rep(NA, length(wh)), 
          rep("Remove", length(wh))))

# Village_Dog_Brazil
gr = which(data[,1] == "Village_Dog_Brazil")
# Remove 3 dogs that don't cluster with the others: PFZ3H03, PFZ3H02, PFZ1B09.
wh = gr[gr < 3000 | gr > 6000]
changed = rbind(changed, cbind(data[wh,1:2], rep(NA, length(wh)), 
          rep("Remove", length(wh))))

# Village_Dog_Colombia
gr = which(data[,1] == "Village_Dog_Colombia")
# Remove 2 dogs that don't cluster with the others: PFZ3G01, PFZ3F02.
wh = gr[gr < 2000 | gr > 5000]
changed = rbind(changed, cbind(data[wh,1:2], rep(NA, length(wh)), 
          rep("Remove", length(wh))))

# Village_Dog_DRC
gr = which(data[,1] == "Village_Dog_DRC")
# There are 2 clusters here, one near 1000 and one near 5178.
# Make 2 clusters called Village_Dog_DRC_1 and Village_Dog_DRCf_2.
wh = gr[gr < 2000]
new.name = "Village_Dog_DRC_1"
changed = rbind(changed, cbind(data[wh,1:2], rep(new.name, length(wh)), 
          rep("Change", length(wh))))
wh = gr[gr > 4000]
new.name = "Village_Dog_DRC_2"
changed = rbind(changed, cbind(data[wh,1:2], rep(new.name, length(wh)), 
          rep("Change", length(wh))))

# Village_Dog_Dominican_Republic
gr = which(data[,1] == "Village_Dog_Dominican_Republic")
# Remove 1 dog that doesn't cluster with the others: PFZ2C06.
wh = gr[gr < 4000]
changed = rbind(changed, c(data[wh,1:2], NA, "Remove"))

# Village_Dog_Egypt
gr = which(data[,1] == "Village_Dog_Egypt")
# Remove 1 dog that doesn't cluster with the others: 1664.
wh = gr[gr < 500]
changed = rbind(changed, c(data[wh,1:2], NA, "Remove"))

# Village_Dog_Fiji
gr = which(data[,1] == "Village_Dog_Fiji")
# Remove 2 dogs that don't cluster with the others: FJ37 & FJ10.
wh = gr[gr < 700 | gr > 4600]
changed = rbind(changed, cbind(data[wh,1:2], rep(NA, length(wh)), 
          rep("Remove", length(wh))))

# Village_Dog_French_Polynesia
gr = which(data[,1] == "Village_Dog_French_Polynesia")
# Remove the dogs in the clusters less than index = 4000.
wh = gr[gr < 4000]
changed = rbind(changed, cbind(data[wh,1:2], rep(NA, length(wh)), 
          rep("Remove", length(wh))))

# Village_Dog_India
gr = which(data[,1] == "_Dog_India")
# Remove the dogs in the low clusters less than index = 500.
wh = gr[gr < 500]
changed = rbind(changed, cbind(data[wh,1:2], rep(NA, length(wh)), 
          rep("Remove", length(wh))))

# Village_Dog_Indonesia
gr = which(data[,1] == "Village_Dog_Indonesia")
# Remove the dogs in the clusters greater than index = 1500.
wh = gr[gr > 1500]
changed = rbind(changed, cbind(data[wh,1:2], rep(NA, length(wh)), 
          rep("Remove", length(wh))))

# Village_Dog_Namibia
gr = which(data[,1] == "Village_Dog_Namibia")
# Remove the dogs in the clusters greater than index = 4000.
wh = gr[gr > 4000]
changed = rbind(changed, c(data[wh,1:2], NA, "Remove"))

# Village_Dog_Nepal
gr = which(data[,1] == "Village_Dog_Nepal")
# Remove the dogs in the clusters less than index = 1100.
wh = gr[gr < 1100]
changed = rbind(changed, cbind(data[wh,1:2], rep(NA, length(wh)), 
          rep("Remove", length(wh))))

# Village_Dog_Peru
gr = which(data[,1] == "Village_Dog_Peru")
# Remove the dogs in the clusters not near 4500
wh = gr[gr < 3000 | gr > 5500]
changed = rbind(changed, cbind(data[wh,1:2], rep(NA, length(wh)), 
          rep("Remove", length(wh))))

# Village_Dog_Portugal
gr = which(data[,1] == "Village_Dog_Portugal")
# Remove the dogs in the clusters not near 4500
wh = gr[gr < 3000 | gr > 5500]
changed = rbind(changed, cbind(data[wh,1:2], rep(NA, length(wh)), 
          rep("Remove", length(wh))))

# Village_Dog_Puerto_Rico
gr = which(data[,1] == "Village_Dog_Puerto_Rico")
# Remove the dogs in the clusters not near 4500.
wh = gr[gr > 4600]
changed = rbind(changed, cbind(data[wh,1:2], rep(NA, length(wh)), 
          rep("Remove", length(wh))))

# Village_Dog_Solomon_Islands
gr = which(data[,1] == "Village_Dog_Solomon_Islands")
# Remove the dogs in the clusters not near 4500.
wh = gr[gr < 2000]
changed = rbind(changed, c(data[wh,1:2], NA, "Remove"))

# Xoloitzcuintli
gr = which(data[,1] == "Xoloitzcuintli")
# Remove the dogs in the clusters not near 4500.
wh = gr[gr > 5500]
changed = rbind(changed, c(data[wh,1:2], NA, "Remove"))

# Yorkshire_Terrier
gr = which(data[,1] == "Yorkshire_Terrier")
# There are two clusters of Yorkshire_Terriers, one near 300 and one near 4300.
# Make 2 clusters called Yorkshire_Terrier_1 and Yorkshire_Terrier_2.
wh = gr[gr < 3000]
new.name = "Yorkshire_Terrier_1"
changed = rbind(changed, cbind(data[wh,1:2], rep(new.name, length(wh)), 
          rep("Change", length(wh))))
wh = gr[gr > 3000]
new.name = "Yorkshire_Terrier_2"
changed = rbind(changed, cbind(data[wh,1:2], rep(new.name, length(wh)), 
          rep("Change", length(wh))))


# Change sample names.
change = grep("Change", changed[,4])
new.names = changed[change,3]
m = match(changed[change, "sample"], data[,2])
data[m,1] = new.names
rownames(K)[m] = new.names
colnames(K)[m] = new.names

# Remove samples.
remove = grep("Remove", changed[,4])
m = match(changed[remove,"sample"], data[,2])
data = data[-m,]
K = K[-m, -m]

# Save the data.
saveRDS(data, file = paste0(ref.dir, "ref_genotypes_cleaned.rds"))
saveRDS(K,    file = paste0(ref.dir, "ref_kinship_cleaned.rds"))

# Save the file containing the changes that we made.
write.table(changed, file = paste0(ref.dir, "ref_changed_samples.txt"), 
            sep = "\t", quote = F)

cl = hclust(as.dist(1.0 - K), method = "average")

png(paste0(fig.dir, "global_kinship_clust_after_changes.png"), width = 25000, 
    height = 2000)
plot(cl, cex.axis = 0.7)
dev.off()

pdf(paste0(fig.dir, "global_kinship_breed_order_after_changes.pdf"))

  unique.breeds = sort(unique(data[,1]))
  for(br in unique.breeds) {
    wh = which(data[,1] == br)
    plot(wh, pch = 16, main = br, col = 0)
    text(1:length(wh), wh, labels = data[wh,1])
  } # for(i)

dev.off()

# Interpolate marker cM.
# Some of the cM data is not ordered by physical position. Find these markers
# and linearly interpolate between the nearest markers that are in order.
# Arguments: map: a 4 column data.frame for one chromosome.
interpolate_markers = function(map) {

  # Get markers where the difference between the current marker and the next is negative.
  wh = which(diff(map[,3]) < 0) + 1

  # If everything is sorted, return the original data.
  if(length(wh) == 0) {
    return(map)
  }
  
  # Handle boundary cases.
  if(wh[1] == 1) {

    map[1,3] = map[2,3] - 0.001
    wh = wh[-1]

  } # if(wh[1] == 1)

  if(wh[length(wh)] == nrow(map)) {

    mean.cm.per.bp = mean(diff(map[,3])) / mean(diff(map[,4]))
    map[nrow(map), 3] = map[nrow(map)-1, 3] + diff(map[(nrow(map)-1):nrow(map), 4]) * mean.cm.per.bp
    wh = wh[-length(wh)]

  } # if(wh[length(wh)] == nrow(map))

  for(i in wh) {

    prev = i - 1
    post = i + 1  # Note: 'next' is a reserved keyword, so I can't use it.

    # This only runs when there are more than one marker out of order in a row.
    if(map[post,3] < map[prev,3]) {

      while(map[post,3] < map[prev,3]) {

        post = post + 1

      } # while(map[post,3] < map[prev,3])

      nmarkers = post - prev - 1
      spacing = (map[post,3] - map[prev,3]) / (nmarkers + 1)
      map[i:(post-1),3] = map[prev,3] + spacing * 1:nmarkers

    } # if(map[post,3] < map[prev,3])

    map[i,3] = map[prev,3] + (map[post,3] - map[prev,3]) * (map[i,4] - map[prev,4]) / (map[post,4] - map[prev,4])

  } # for(i)

  wh = which(diff(map[,3]) < 0)
  if(length(wh) > 0) {

    stop("Number of markers is out of order after sorting.")

  } # if(length(wh) > 0)

  return(map)

} # interpolate_markers()

for(i in 1:39) {

  rng = which(markers[,1] == i)
  markers[rng,] = interpolate_markers(map = markers[rng,])

} # for(i)


# Write out data with different subsets.
# 1. AllLQ: All data, all markers and samples.
# 2. AllHQ: All data, < 5% no-call.
# 3. PureBreedsLQ: Pure breeds, all markers and samples.
# 4. PureBreedsHQ: Pure breeds, < 5% no-call.

write_plink = function(data, markers, out.dir) {

  # Write out the data by chromosome in PLINK format.
  unique.chr = unique(markers[,1])

  print(paste("Num. samples =", nrow(data)))
  print(paste("Num. markers =", nrow(markers)))

  for(i in 1:length(unique.chr)) {

    print(paste("CHR", unique.chr[i]))

    ss = which(markers[,1] == unique.chr[i])
    write.table(markers[ss,], file = paste0(out.dir, "_chr", unique.chr[i], ".map"), 
                sep = "\t", row.names = F, col.names = F, quote = F)

    ss.data = c(1:6, ss + 6)
    write.table(data[,ss.data], file = paste0(out.dir, "_chr", unique.chr[i], ".ped"), 
                sep = "\t", row.names = F, col.names = F, quote = F)

  } # for(i)

} # write_plink()

# 1. AllLQ: All data, all markers and samples.
write_plink(data = data, markers = markers,
            out.dir = "/hpcdata/dgatti/Dog/ReferenceData/AllLQ/AllLQ")

# 2. AllHQ: All data, < 5% no-call.
sel.mkr = which(colMeans(data == "0 0") < 0.05)
sel.sam = which(rowMeans(data == "0 0") < 0.05)
# Note: data has 6 more columns than markers.
write_plink(data = data[sel.sam, sel.mkr], markers = markers[sel.mkr[-(1:6)] - 6,],
            out.dir = "/hpcdata/dgatti/Dog/ReferenceData/AllHQ/AllHQ")

# 3. PureBreedsLQ: Pure breeds, all markers and samples.
data = data[-grep("Mix|Village|Unknown", data[,1]),]
write_plink(data = data, markers = markers,
            out.dir = "/hpcdata/dgatti/Dog/ReferenceData/PureBreedsLQ/PureBreedsLQ")

# 4. PureBreedsHQ: Pure breeds, < 5% no-call. Breeds with >= 10 samples
sel.mkr = which(colMeans(data == "0 0") < 0.05)
sel.sam = which(rowMeans(data == "0 0") < 0.05)
tbl = table(data[sel.sam,1])
tbl = tbl[tbl >= 10]
sel.sam = intersect(sel.sam, which(data[,1] %in% names(tbl)))
stopifnot(min(table(data[sel.sam,1])) == 10)
# Keep the wolves.
sel.sam = c(grep("Gray_Wolf",data[,1]), sel.sam)
# Keep the American Eskimo Dog.
sel.sam = c(grep("American_Eskimo_Dog",data[,1]), sel.sam)

write_plink(data = data[sel.sam, sel.mkr], markers = markers[sel.mkr[-(1:6)] - 6,],
            out.dir = "/hpcdata/dgatti/Dog/ReferenceData/PureBreedsHQ/PureBreedsHQ")




