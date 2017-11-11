################################################################################
# Compare the number of recombinations per dog and the number of breeds needed
# to account fo 80% of the genome between window sizes.
# Daniel Gatti
# dan.gatti@jax.org
# Aug. 16, 2017
################################################################################
options(stringsAsFactors = F)

base.dir = "/projects/dgatti/Dog/HMM/HMM_v7_SHAPEIT/"
fig.dir = "/projects/dgatti/Dog/Figures/window_size/"

# Count the number of recombinations in a viterbi object.
count_recomb = function(viterbi) {

  viterbi = lapply(viterbi, function(z) { diff(z[,1]) })
  viterbi = sapply(viterbi, function(z) { z[z != 0] })

  return(length(unlist(viterbi)))

} # count_recomb()

# Get the mean Viterbi probability across the genome.
mean_viterbi_probs = function(viterbi) {

  viterbi = lapply(viterbi, function(z) { z[,2] })
  return(mean(unlist(viterbi)))

} # mean_viterbi_probs()

# Count the number of breeds needed to account for 90% of the
# genome.
count_breeds = function(viterbi, prop) {

  viterbi = unlist(lapply(viterbi, function(z) { z[,1] }))
  vit.tbl = sort(table(viterbi), decreasing = T)
  vit.tbl = vit.tbl / sum(vit.tbl)
  return(min(which(cumsum(vit.tbl) > prop)))

} # count_breeds()


# For each directory, read in the viterbi files and count the number of
# recombinations and the percentage composition of the top 4 breeds.
summarize_dogs = function(path) {

  print(path)

  # Get the viterbi files.
  vit.files = dir(path = path, pattern = "viterbi.rds$", full.names = T)

  # Get the dog names.
  dog.names = unique(gsub("_(A|B)_HMM_viterbi.rds$", "", vit.files))

  # Make summary data structures.
  num.recomb = vector("numeric", length(dog.names))
  names(num.recomb) = dog.names
  prop80 = vector("numeric", length(dog.names))
  names(prop80) = dog.names
  mean.prob = vector("numeric", length(dog.names))
  names(mean.prob) = dog.names

  for(d in dog.names) {

    # Read in the viterbi files.
    dog.vit.files = vit.files[grep(paste0("^", d, "_"), vit.files)]
    stopifnot(length(dog.vit.files) == 2)

    vit = c(readRDS(dog.vit.files[1]), readRDS(dog.vit.files[2]))
    num.recomb[d] = count_recomb(vit)
    prop80[d] = count_breeds(vit, 0.8)
    mean.prob[d] = mean_viterbi_probs(vit)

  } # for(d)

  return(list(num.recomb, prop80, mean.prob))

} # summarize_dogs()

# Get the directories for each window size.
win.dirs = list.dirs(path = base.dir)
win.dirs = win.dirs[grep("win[0-9]", win.dirs)]

win1 = summarize_dogs(win.dirs[1])
win7 = summarize_dogs(win.dirs[2])

# Make figures of the results.
png(paste0(fig.dir, "win1_win7_comparison.png"), width = 1000, height = 1000,
    res = 128)
layout(matrix(1:4, 2, 2))
lim = range(c(win1[[1]], win7[[1]]))
plot(win1[[1]], win7[[1]], xlab = "win1", ylab = "win7", xlim = lim, ylim = lim,
     main = "Number of Recombinations")
abline(0, 1, col = 2)

lim = range(c(win1[[2]], win7[[2]]))
plot(win1[[2]], win7[[2]], xlab = "win1", ylab = "win7", xlim = lim, ylim = lim, 
     main = "Number of Breeds to Explain 80%")
abline(0, 1, col = 2)

lim = range(c(win1[[3]], win7[[3]]))
plot(win1[[3]], win7[[3]], xlab = "win1", ylab = "win7",  xlim = lim, ylim = lim,
     main = "Mean Viterbi Probability")
abline(0, 1, col = 2)
dev.off()



