################################################################################
# Simulate one dog genome using the reference data.
# Daniel Gatti
# dan.gatti@jax.org
# June 27, 2017
# Contains: select_chunk_bounds(), sample_breeds(), read_vcf_files(), 
#           read_map_files(), simulate_genome().
################################################################################
options(stringsAsFactors = F)
library(VariantAnnotation)

################################################################################
# Select the chunk positions based on the given chunk sizes and probabilities.
# Arguments:
# chunk: numeric vector containing the discrete chunk sizes to sample.
# probs: numeric vector containing the probability of drawing each chunk size.
# map: data.frame containing the markers positions in PLINK format.
select_chunk_bounds = function(chunk, probs, map) {

#  chrlen = floor(map[nrow(map), 4])
  chrlen = map[nrow(map), 4]

  # We have to handle a single chunk size differently becuase of the way
  # that sample() works.
  if(length(chunk) == 1) {

    n.chunks = ceiling(chrlen / chunk)
    bounds = seq(0, chrlen, chunk)
    if(n.chunks == 1) {

      bounds = c(0, chrlen)

    } # if(n.chunks == 1)

  } else {

    n.chunks = ceiling(map[nrow(map), 4] / min(chunk))
    sizes = sample(x = chunk, size = n.chunks, replace = TRUE, prob = probs)
    sizes = c(0, sizes[cumsum(sizes) <= chrlen])
    bounds = cumsum(sizes)

  } # else

  # Make sure the the end of the chromsome lines up.
  bounds[length(bounds)] = chrlen

  return(bounds)

} # select_chunk_bounds()


################################################################################
# Given the reference genotypes and the location of the chunk, select a
# sample from a breed with uniform probability.
# We sample chunks from each chromatid and then assemble them.
# hap: Numeric matrix containing the phased reference haplotypes. Haplotypes
#      in rows and markers in columns.
# fam: Data.frame containing sample information n PLINK FAM format.
# chunk: numeric vector containing the discrete chunk sizes to sample.
# probs: numeric vector containing the probability of drawing each chunk size.
# map: data.frame containing the markers positions in PLINK format.
sample_breeds = function(hap, fam, chunk, probs, map) {

  breeds = factor(fam[,1], levels = sort(unique(fam[,1])))

  # SHAPEIT uses the hap file format. Get the haplotype breeds by stripping
  # off the trailing _(A|B)
  hap.breeds = sub("_(A|B)$", "", rownames(hap))

  obs   = matrix(0, nrow(map), ncol = 2, dimnames = list(map[,2], 1:2))
  truth = matrix(0, nrow(map), ncol = 2, dimnames = list(map[,2], 1:2))

  for(i in 1:2) {

    # Note: mulitply bounds by 1e6 to be in bp.
    bounds = select_chunk_bounds(chunk = chunk, probs = probs, map = map) * 1e6
    n.chunks = length(bounds) - 1
    br = sample(x = breeds, size = n.chunks, replace = T)

    # Multiply the map Mb by 1e6 to get bp.
    map.bp = map[,4] * 1e6

#    rng = 1
    for(j in 1:n.chunks) {

      # Select a dog from the current breed.
      wh = which(hap.breeds == br[j])
      wh = sample(wh, 1)

      mkr.rng = which(map.bp > bounds[j] & map.bp <= bounds[j+1])

      if(length(mkr.rng) > 0) {

        chunk.gt = hap[wh,mkr.rng]
 #       rng = rng:(rng + length(chunk.gt) - 1)

        truth[mkr.rng, i] = rep(br[j], length(mkr.rng))
        obs[mkr.rng, i]   = chunk.gt
#        rng = rng[length(rng)] + 1

      } # if(length(mkr.rng) > 0)

    } # for(j)

  } # for(i)

  retval = list(truth = truth, obs = obs)
  attr(retval, "breeds") = breeds

  return(retval)

} # sample_breeds()


################################################################################
# Read in the VCF files for all chromosomes.
# Argument:
# ref.dir: full path to the directory containing the phased reference VCF files.
# Return a list containing vcf objects.
read_vcf_files = function(ref.dir) {

  vcf.files = dir(path = ref_dir, pattern = paste0("chr[0-9]+.phased.vcf.gz$"),
              full.names = T)

  spl = sub(".phased.vcf.gz$", "", vcf.files)
  spl = strsplit(spl, "chr")
  chrs = as.numeric(sapply(spl, "[", 2))
  vcf.files = vcf.files[order(chrs)]
  chrs = sort(chrs)

  vcf.list = as.list(vcf.files)
  names(vcf.list) = chrs

  for(chr in 1:length(vcf.files)) {

    print(paste("CHR", chr ))

    vcf.list[[chr]] = readVcf(vcf.files[chr], genome = "CanFam3.1")

  } # for(chr)

  return(vcf.list)

} # read_vcf_files()


################################################################################
# Read in the MAP files for all chromosomes.
# ref.dir: full path to the directory containing the MAP files in PLINK format.
# Return a list containing data.frames with MAPs. NOTE: positions converted
# to Mb.
read_map_files = function(ref.dir) {

  map.files = dir(path = ref_dir, pattern = paste0("map$"), full.names = T)

  spl = sub(".map$", "", map.files)
  spl = strsplit(spl, "chr")
  chrs = as.numeric(sapply(spl, "[", 2))
  map.files = map.files[order(chrs)]
  chrs = sort(chrs)

  map.list = lapply(map.files, read.delim, header = F)

  for(i in 1:length(map.list)) {

    # Convert marker positions from bp to Mb.
    map.list[[i]][,4] = map.list[[i]][,4] * 1e-6

  } # for(i)

  return(map.list)

} # read_map_files()



################################################################################
# Simulate a given number of dog genomes on a single chromosome.
# Arguments:
# n: integer that is the number of samples to simulate. Default = 1.
# chunk_sizes: numeric vector of chunk sizes in Mb to sample from the reference 
#              data. Default = 20.
# chunk_probs: numeric vector of probabilities of same length as chunk_sizes
#              that represents the probability of sampling the corresponding
#              chunk. Must sum to 1. Default = 1.
# hap.file: full path to the phased reference *.rds file.
# map.file: full path to the phased reference MAP file.
# fam.file: full path to the phased reference FAM file.
# phased: logical variable that is TRUE if the returned genotypes should remain
#         phased and FALSE if they should be scrambled. Default = TRUE.
# Returns: data.frame in PLINK format with the phased genotypes for one dog.
simulate_genome = function(n = 1, chunk_sizes = 20,
                  chunk_probs = rep(1/length(chunk_sizes), length(chunk_sizes)),
                  hap.file, map.file, fam.file, phased = TRUE) {

  if(length(chunk_sizes) != length(chunk_probs)) {
    stop(paste("length(chunk_sizes) != length(chunk_probs)\nPlease make",
         "sure that these are both of equal length."))
  } # if(length(chunk_sizes) != length(chunk_probs))

  if(abs(sum(chunk_probs) - 1.0) > 1e-8) {
    stop(paste("sum(chunk_probs) != 1.0\nPlease make",
         "sure that chunk_probs sums to 1."))
  } # if(abs(sum(chunk_probs) - 1.0) > 1e-8)

  # Read in the haplotype, MAP and FAM files.
  print("Reading files...")
  hap = readRDS(hap.file)
  map = read.delim(map.file, header = F, sep = "")
  map[,4] = map[,4] * 1e-6

  # Read in FAM file to get the breed info.
  fam = read.delim(fam.file, sep = " ", header = F, stringsAsFactors = F)
  breeds = sort(unique(fam[,1]))

  stopifnot(ncol(hap) == nrow(map))
  stopifnot(nrow(hap) == 2 * nrow(fam))

  # We select from breeds with uniform probability. Then select a sample from
  # within a breed with uniform probability.
  geno  = matrix(0, nrow = 2 * n, ncol = nrow(map), dimnames = 
          list(paste(rep(paste0("sim", 1:n), each = 2), c("A", "B"), sep = "_"), map[,2]))
  truth = matrix(0, nrow = nrow(geno), ncol = ncol(geno), dimnames = dimnames(geno))

  for(i in 1:n) {

    print(paste(i, "of", n))

    # Select chunk boundaries.
    tmp = sample_breeds(hap = hap, fam = fam, chunk = chunk_sizes,
          probs = chunk_probs, map = map)

    if(phased == FALSE) {

      swap = sample(1:nrow(tmp[[2]]), nrow(tmp[[2]]) / 2)
      tmp[[2]][swap,] = tmp[[2]][swap,2:1]

    } # if(phased == FALSE)

    sim.smaple.rng = (2 * (i-1) + 1):(2 * i)
    truth[sim.smaple.rng,] = t(tmp[[1]])
    geno[sim.smaple.rng,]  = t(tmp[[2]])

  } # for(i)

  return(list(truth = truth, obs = geno, breeds = breeds))

} # simulate_genome()

