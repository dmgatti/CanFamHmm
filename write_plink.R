################################################################################
# Given a family data.frame, a data matrix with space-delimited allele calls
# and a marker map, write out a set of PLINK PED and MAP files.
# Daniel Gatti
# dan.gatti@jax.org
# July 16, 2017
################################################################################
# Arguments:
# fam: data.frame containing PLINK family data. num.samples x 6.
# data: matrix containing space-delimited allele calls. num.samples x num.markers.
# map: data.frame containing marker information in PLINK format. num.markers x 4.
# out.prefix: string containing the prefix to append to the PLINK file
#             suffixes.
write_plink = function(fam, data, map, file.prefix) {

  # Error checks.
  num.samples = nrow(fam)
  num.markers = nrow(map)
  stopifnot(num.samples == nrow(data))
  stopifnot(num.markers == ncol(data))
  stopifnot(ncol(fam) == 6)
  stopifnot(ncol(map) == 4)
  
  # Write out PED file.
  data = cbind(fam, data)
  write.table(data, file = paste0(file.prefix, ".ped"), sep = "\t", quote = F, 
              row.names = F, col.names = F)

  # Write out MAP file.
  write.table(map, file = paste0(file.prefix, ".map"), sep = "\t", quote = F, 
              row.names = F, col.names = F)  

} # write_plink()