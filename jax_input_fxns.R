################################################################################
# Read in the JAX data from Geneseek, remove samples with low call rates,
# move the markers over to CanFam3 and write out to an RDS file.
# Daniel Gatti
# dan.gatti@jax.org
# Apr. 4, 2017
################################################################################
options(stringsAsFactors = FALSE)
setwd("/hpcdata/dgatti/Dog/")

# Create a function that reads in Geneseek data, cleans it, and outputs one
# PLINK file.
# Arguments: input.dir: Path to a directory containing a single Geneseek 
#                       data set.
#            output.dir: Path to the location for output files. Default =
#                        input.dir.
#            canFam3: data.frame containing the CanFam3.1 marker locations.
#                     Columns contain chr, marker name, cM & bp.
extract_jax_data = function(input.dir, output.dir = input.dir, 
                   canFam3) {

  # Sort the canFam3 markers.
  canFam3 = canFam3[order(canFam3[,1], canFam3[,4]),]

  # Load in the Illumina data.
  print("Reading FinalReport file...")
  final.report.file = dir(path = input.dir, pattern = "_FinalReport.txt$",
                      full.names = TRUE)
  jax = scan(final.report.file, what = "char", sep = "\n", skip = 9)
  jax = strsplit(jax, split = "\t")
  jax = matrix(unlist(jax), nrow = length(jax), ncol = length(jax[[1]]),
        dimnames = list(NULL, jax[[1]]), byrow = T)
  jax = jax[-1,]

  # Create a matrix like the PLINK PED matrix.
  jax = split(data.frame(jax), jax[,2])

  # Place data into a matrix.
  plink.colnames = c("family", "sample", "father", "mother", "sex", "phenotype")
  data = matrix("", nrow = length(jax), ncol = 6 + nrow(jax[[1]]),
         dimnames = list(names(jax), c(plink.colnames, jax[[1]][,1])))
  x = matrix(0, nrow(data), ncol(data), dimnames = dimnames(data))[,-(1:6)]
  y = matrix(0, nrow(data), ncol(data), dimnames = dimnames(data))[,-(1:6)]

  data[,1] = names(jax)
  data[,2] = names(jax)
  data[,3] = 0
  data[,4] = 0
  data[,5] = 0
  data[,6] = -9

  print("Extracting Alleles ...")
  for(i in 1:length(jax)) {

    data[i,-(1:6)] = paste(jax[[i]][,"Allele1...Top"], jax[[i]][,"Allele2...Top"], 
                   sep = " ")
    x[i,] = as.numeric(jax[[i]][,"X"])
    y[i,] = as.numeric(jax[[i]][,"Y"])

  } # for(i)

  # Remove samples with no-call rates > 20%.
  call.rate = rowMeans(data == "- -")
  write.table(call.rate, file = paste0(output.dir, "call_rate.txt"),
              sep = "\t")
  remove = which(call.rate > 0.2)
  if(length(remove > 0)) {

    print(paste("Removing samples", paste(data[remove,1], collapse = ","),
          "due to high no-call rates."))
    data = data[-remove,]
    x = x[-remove,]
    y = y[-remove,]

  } # if(length(remove > 0)  

  # Order the markers by genome location.
  hdr = data[,1:6]
  data = data[,-(1:6)]
  data = data[,match(canFam3[,2], colnames(data))]
  stopifnot(colnames(data) == canFam3[,2])
  data = cbind(hdr, data)

  x = x[,canFam3[,2]]
  y = y[,canFam3[,2]]

  # Get the sample sex. NOTE: PLINK uses 1 = male and 2 = female.
  sex = rep("F", nrow(data))
  names(sex) = rownames(data)
  chrx = which(canFam3[,1] == 39)  # Chr 39 is distal Chr X
  chry = which(canFam3[,1] == 40)  # Chr 40 is Chr Y
  rmx = rowMeans(x[,chrx] + y[,chrx], na.rm = T)
  rmy = rowMeans(x[,chry] + y[,chry], na.rm = T)
  plot(rmx, rmy)
  sex[rmy > 0.6] = "M"
  sex = factor(sex, levels = c("M", "F"))
  data[,5] = as.numeric(sex)

  # Change no-call symbol from "-" to 0 for PLINK.
  data[data == "- -"] = "0 0"

  # Write out data.
  print(paste("Write out data to", output.dir))
  saveRDS(x, file = paste0(output.dir, "jax_cleaned_x.rds"))
  saveRDS(y, file = paste0(output.dir, "jax_cleaned_y.rds"))
  saveRDS(data,    file = paste0(output.dir, "jax_cleaned_data.rds"))
  saveRDS(canFam3, file = paste0(output.dir, "jax_cleaned_markers.rds"))

} # extract_jax_data()


# Read in the RDS files in the input directory and write them out in 
# plink PED and MAP format.
# Bonus: convert to BED/BIM/FAM format.
convert2plink = function(input.dir, camFam3) {

} # convert2plink()