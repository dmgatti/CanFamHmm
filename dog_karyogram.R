################################################################################
# Dog kayrogram.
options(stringsAsFactors = F)
library(RColorBrewer)

# Arguments: map: marker annotation file in PLINK format.
# NOTE: We only plot Chr 1 through 39 (X).
chr.skeletons = function(map,  ...) {
  
  chrlen = sapply(map, max)
  names(chrlen)[length(chrlen)] = "X"

  # Plot a coordinate system to draw on.
  par(font = 2, font.axis = 2, font.lab = 2, las = 1, plt = c(0.07, 0.99, 
      0.08, 0.99))
  plot(-1, -1, col = 0, xlim = c(0, max(chrlen) + 10), 
       ylim = c(0.5, length(chrlen) + 0.5), cex.axis = 1.25, yaxs = "i",
       xaxs = "i", yaxt = "n", xlab = "", ylab = "", ...)

  # Draw chr skeletons
  for(i in 1:length(chrlen)) {
    lines(c(0, chrlen[i]), rep(length(chrlen) - i + 1, 2), col = 0)
  } # for(i)

  # Draw light lines at 10 Mb and heavier lines at 50 Mb.
  abline(v = 0:30 * 10, col = "grey80")
  abline(v = 0:5  * 50, col = "grey50", lwd = 1.2)
  mtext(side = 2, line = 0.25, at = length(chrlen):1, cex = 1.25, 
        text = names(chrlen))
  mtext(side = 1, line = 2.5, cex = 1.5, text = "Megabases")
  mtext(side = 2, line = 2, font = 2, cex = 1.5, text = "Chromosome", las = 0)

} # chr.skeletons()


# Set up breed colors. This will create the same colors for all breeds.
get_breed_colors = function(breed.names) {

  breed.colors = colors()
  breed.colors = breed.colors[-grep("gr[ae]y|ivory|white|beige", breed.colors)]
  breed.colors = breed.colors[seq(2, length(breed.colors), length.out = 
                 length(breed.names))]

  names(breed.colors) = breed.names

  return(breed.colors)
 
} # get_breed_colors()


# This will just return colors for up to 12 breeds.
get_breed_colors2 = function(breed.names) {

  if(length(breed.names) > 12) {
    stop("get_breed_colors2: can't handle more than 12 colors.")
  }

#  col = c("#F68073", "#EBAC0F", "#128B21", "#527491", "#B96107", "#262E75", 
#          "#DC3F12", "#48BE98", "#F0AF4D", "#61445D", "#AFE374", "#2F1510")

  col = c("#33a02c", "#ff7f00", "#1f78b4", "#e31a1c", "#6a3d9a", "#b15928",
          "#a6cee3", "#fb9a99", "#fdbf6f", "#b2df8a", "#cab2d6", "#efcf33")

  breed.colors = col[1:length(breed.names)]
  names(breed.colors) = breed.names

  return(breed.colors)
 
} # get_breed_colors2()


################################################################################
# Pseudophase the two chromosomes by minimizing the number of recombinations.
# Arguments: haps: numeric matrix with breed indices for each haplotype. Named
#                  markers in rows and 2 haplotypes in columns. Breed IDs are
#                  1-based and must match the breeds in "breeds".
#            breeds: Vector of breed names, ordered in the same way that the
#                    breeds are numbered in haps.
################################################################################
pseudophase = function(haps, breeds) {

  # Create a matrix with counts of the nubmer of alleles for each breed at
  # each marker.
  mat = matrix(0, nrow = nrow(haps), ncol = length(breeds), dimnames =
        list(names(haps[,1]), breeds))
  idx.mat = cbind(1:nrow(haps), haps[,1])
  mat[idx.mat] = 1
  idx.mat = cbind(1:nrow(haps), haps[,2])
  mat[idx.mat] = mat[idx.mat] + 1

  # Set haplotype A.
  idx = 1
  hap1 = rep(0, nrow(mat))
  while(idx < nrow(mat)) {

    br = min(which(mat[idx,] > 0))
    start = idx
    end = min(which(diff(c(mat[idx:nrow(mat),br], 0)) < 0)) + idx - 1
    hap1[start:end] = br
    mat[start:end,br] = mat[start:end,br] - 1
    idx = end + 1

  } # while(idx < nrow(mat))

  # Set haplotype B.
  hap2 = rep(0, nrow(mat))
  idx = 1
  while(idx < nrow(mat)) {

    br = min(which(mat[idx,] > 0))
    start = idx
    end = min(which(diff(c(mat[idx:nrow(mat),br], 0)) < 0)) + idx - 1
    hap2[start:end] = br
    mat[start:end,br] = mat[start:end,br] - 1
    idx = end + 1

  } # while(idx < nrow(mat))
 
  return(cbind(hap1, hap2))

} # pseudophase()


# Fill in the haplotypes for a single vector of haplotypes on one chromosome.
fill_hap = function(hap) {

  brks = diff(is.na(hap))
  st = which(brks > 0)
  en = which(brks < 0)

  if(length(c(st,en)) > 0) {

    # Handle boundary conditions at the beginning and end.
    if(is.na(hap[1])) {
      st = c(1, st)
    } # if(is.na(hap[1]))

    if(is.na(hap[length(hap)])) {
      en = c(en, length(hap))
    } # if(is.na(hap[1]))

    stopifnot(length(st) == length(en))
    stopifnot(st < en)

    # For each NA block, fill in the block with the surrounding breeds.
    for(j in 1:length(st)) {

      # Proximal breed. If it's 0, then this block extends to the beginning
      # of the haplotype. Use the distal breed.
      prox = hap[st[j] - 1]      
      if(st[j] == 1) {
        prox = hap[en[j] + 1]
      } 

      # Distal breed. If it's 0, then this block extends to the end
      # of the haplotype. Use the proximal breed.
      dist = hap[en[j] + 1]
      if(en[j] == length(hap)) {
        dist = hap[st[j] - 1]
      }

      # Fill in with surrounding breeds from start to mid and mid to end.
      mid = round((en[j] + st[j]) / 2)
      hap[st[j]:mid] = prox
      hap[mid:en[j]] = dist

    } # for(j)

  } # if(length(st) > 0)

  stopifnot(all(!is.na(hap)))

  return(hap)

} # fill_hap()



################################################################################
# Identify short regions from breeds that comprise a small fraction of the 
# genome and write over these regions with the proximal and distal breeds.
# Arguments: hapA: list containing the Viterbi path for haplotype A.
#                  One per chromosome, in order. Probs already removed.
#            hapB: list containing the Viterbi path for haplotype B.
#                  One per chromosome, in order. Probs already removed.
################################################################################
simplify_haps = function(hapA, hapB) {

  # Get the overall breed composition.
  haps = unlist(hapA, hapB)
  tbl = table(haps)
  tbl = sort(tbl / sum(tbl), decreasing = T)

  breeds2change = as.numeric(names(tbl)[tbl <= 0.01])

  hapA = lapply(hapA, function(z) { z[z %in% breeds2change] = NA; z })
  hapB = lapply(hapB, function(z) { z[z %in% breeds2change] = NA; z })

  # If any chromosomes are completely NA, set them to the top 3 breeds.
  for(i in 1:length(hapA)) {

    if(all(is.na(hapA[[i]]))) {

      top.breeds = as.numeric(names(tbl)[1:3])
      len = length(hapA[[i]])
      hapA[[i]] = c(top.breeds[1], rep(top.breeds, round((tbl[1:3] / 
                    sum(tbl[1:3])) * len)))[1:len]

    } # if(all(is.na(hapA[[i]])))

    if(all(is.na(hapB[[i]]))) {

      top.breeds = as.numeric(names(tbl)[1:3])
      len = length(hapB[[i]])
      hapB[[i]] = c(top.breeds[1], rep(top.breeds, round((tbl[1:3] / 
                    sum(tbl[1:3])) * len)))[1:len]

    } # if(all(is.na(hapB[[i]])))

  } # for(i)

  # Fill in haps.
  hapA = lapply(hapA, fill_hap)
  hapB = lapply(hapB, fill_hap)

  return(list(hapA = hapA, hapB = hapB))

} # simplify_haps()


get_breed_composition = function(hapA, probA, hapB, probB, breeds, confidence.thr,
                        composition.thr) {

  # Zero out breeds with low probabilities.
  tmpA = mapply(function(x, y) { x[y < confidence.thr] = NA; x }, x = hapA, y = probA)
  tmpB = mapply(function(x, y) { x[y < confidence.thr] = NA; x }, x = hapB, y = probB)
  breed = c(unlist(tmpA), unlist(tmpB))
  breed = table(breed)
  names(breed) = breeds[as.numeric(names(breed))]
  breed = breed / sum(breed)
  breed = sort(breed, decreasing = TRUE)
  num.breeds = min(which(cumsum(breed) / sum(breed) >= composition.thr))
#  num.breeds = sum(breed > 0.03)
  num.breeds = min(num.breeds, 12)

  return(list(breed.comp = breed, num.breeds = num.breeds))

} # get_breed_composition()


################################################################################
# Arguments: hapA: list containing the Viterbi path for haplotype A.
#                      One per chromosome, in order.
#            hapB: list containing the Viterbi path for haplotype B.
#                      One per chromosome, in order.
#            map: list containing the marker positions for the haplotypes
#                      on each chromosome. In Mb.
#            breeds: character vector containing the breed names in the same order
#                    as they appear in the viterbi files.
#            title: character string to place at top of plot.
#            sex: either M or F. Changes how Chr X is drawn.
#            legend: boolean variable indicating whether to plot the legend.
#            colors: either the value "global", in which case colors are selected
#                    for all breeds or a set of colors to use, or a single
#                    integer that is the number of breeds to show.
dog_karyogram = function(hapA, hapB, map, breeds, title, sex, legend = T,
                colors = "global") {

  stopifnot(length(hapA) == length(hapB))
  stopifnot(length(hapA) == length(map))

  confidence.thr  = 0.0
  composition.thr = 1.0

  # For now, strip out the probs in column 2, we may use them later.
  probA = lapply(hapA, function(z) { z[,2] })
  probB = lapply(hapB, function(z) { z[,2] })

  hapA = lapply(hapA, function(z) { z[,1] })
  hapB = lapply(hapB, function(z) { z[,1] })

  stopifnot(sapply(hapA, length) == sapply(hapB, length))

  # Get overall breed composition.
  tmp = get_breed_composition(hapA, probA, hapB, probB, breeds, confidence.thr,
        composition.thr)
  breed = tmp$breed.comp
  num.breeds = tmp$num.breeds

  # Pseudophase the haplotypes to place breeds on the same chromosome.
  haps = vector("list", length(hapA))
  for(i in 1:length(hapA)) {

    haps[[i]] = cbind(hapA[[i]], hapB[[i]])

  } # for(i)

  haps = lapply(haps, pseudophase, breeds = breeds)

  # Convert back to hapA & hapB.
  hapA = lapply(haps, function(z) { z[,1] })
  hapB = lapply(haps, function(z) { z[,2] })
  names(hapA) = 1:length(hapA)
  names(hapB) = 1:length(hapA)

  breed.colors = get_breed_colors(breeds)
  if(is.numeric(colors)) {

    breed.colors = get_breed_colors2(names(breed)[1:num.breeds])
    tmp = rep(0, length(breeds))
    names(tmp) = breeds
    tmp[names(breed.colors)] = breed.colors
    breed.colors = tmp

  } # if(is.numeric(colors))

  chr.skeletons(map)
  mtext(side = 3, line = 0.5, font = 2, cex = 2, text = title)

  offset = 0.01
  width  = 0.35

  for(i in 1:length(hapA)) {

    mkrs = map[[i]]

    # Haplotype A
    brks = c(1, which(diff(hapA[[i]]) != 0), length(mkrs))

    rect.bounds = matrix(0, length(brks) - 1, 4, dimnames = 
                  list(NULL, c("xleft", "ybottom", "xright", "ytop")))
    rect.bounds[,"ybottom"] = length(map) - i + 1 + offset
    rect.bounds[,"ytop"] = length(map) - i + 1 + offset + width
    rect.colors = rep(0, nrow(rect.bounds) - 1)
    rect.probs  = rep(0, nrow(rect.bounds) - 1)
    gt = rep(0, nrow(rect.bounds) - 1)
 
    for(j in 2:length(brks)) {

      gt[j-1] = hapA[[i]][brks[j]]
      rect.bounds[j-1, "xleft"]  = mkrs[brks[j-1]]
      rect.bounds[j-1, "xright"] = mkrs[brks[j]]
      rect.colors[j-1] = breed.colors[gt[j-1]]
      rect.probs[j-1] = median(probA[[i]][brks[j-1]:brks[j]])

    } # for(j)

    # Add transparency to the blocks to reflect confidence.
#    rect.probs[rect.probs < confidence.thr]  = 0
#    rect.probs[rect.probs >= confidence.thr] = 1
#    rect.probs[breeds[gt] %in% names(breed)[1:num.breeds]] = 1
#    rect.colors = rgb(t(col2rgb(rect.colors)), alpha = rect.probs * 255,
#                  maxColorValue = 255)
#    rect.colors[rect.probs < confidence.thr] = 0

    rect(xleft = rect.bounds[,1], ybottom = rect.bounds[,2], 
         xright = rect.bounds[,3], ytop = rect.bounds[,4], density = NA,
         col = rect.colors, border = 1)

    if(names(hapA)[i] != "39" | sex == "F") {

      # Haplotype B
      brks = c(1, which(diff(hapB[[i]]) != 0), length(mkrs))

      rect.bounds = matrix(0, length(brks) - 1, 4, dimnames = 
                    list(NULL, c("xleft", "ybottom", "xright", "ytop")))
      rect.bounds[,"ybottom"] = length(map) - i + 1 - offset - width
      rect.bounds[,"ytop"] = length(map) - i + 1 - offset
      rect.colors = rep(0, nrow(rect.bounds) - 1)
      rect.probs  = rep(0, nrow(rect.bounds) - 1)
      gt = rep(0, nrow(rect.bounds) - 1)
 
      for(j in 2:length(brks)) {

        gt[j-1] = hapB[[i]][brks[j]]
        rect.bounds[j-1, "xleft"]  = mkrs[brks[j-1]]
        rect.bounds[j-1, "xright"] = mkrs[brks[j]]
        rect.colors[j-1] = breed.colors[gt[j-1]]
        rect.probs[j-1] = median(probB[[i]][brks[j-1]:brks[j]])

      } # for(j)

      # Add transparency to the blocks to reflect confidence.
#      rect.probs[rect.probs < confidence.thr]  = 0
#      rect.probs[rect.probs >= confidence.thr] = 1
#      rect.probs[breeds[gt] %in% names(breed)[1:num.breeds]] = 1
#      rect.colors = rgb(t(col2rgb(rect.colors)), alpha = rect.probs * 255,
#                    maxColorValue = 255)
      rect.colors[rect.probs < confidence.thr] = 0

      rect(xleft = rect.bounds[,1], ybottom = rect.bounds[,2], 
           xright = rect.bounds[,3], ytop = rect.bounds[,4], density = NA,
           col = rect.colors, border = 1)

    } #  if(i < 39 | sex == "F")

  } # for(i)

  if(legend) {
    usr = par("usr")
    # Note: Using strwidth() to fix the error sizing the box.
    legend(x = usr[1] + 0.6 * diff(usr[1:2]), 
           y = usr[3] + 0.6 * diff(usr[3:4]), 
           legend = gsub("_", " ", names(breed)[1:num.breeds]),
           bg = "white", cex = 1.25, text.width = 44,
           fill = breed.colors[names(breed)[1:num.breeds]], pt.cex = 2,
           title = title)
  } # if(legend)

} # dog_karyogram()



################################################################################
# Plot the genome in two side-by-side plots; chr 1:20 and chr 21:39.
# Arguments: hapA: list containing the Viterbi path for haplotype A.
#                      One per chromosome, in order.
#            hapB: list containing the Viterbi path for haplotype B.
#                      One per chromosome, in order.
#            map: list containing the marker positions for the haplotypes
#                      on each chromosome. In Mb.
#            breeds: character vector containing the breed names in the same order
#                    as they appear in the viterbi files.
#            title: character string to place at top of plot.
#            sex: either M or F. Changes how Chr X is drawn.
dog_karyogram2 = function(hapA, hapB, map, breeds, title, sex) {

  layout(matrix(1:2, 1, 2))
  dog_karyogram(hapA[1:20], hapB[1:20], map[1:20], breeds, title, sex, legend = FALSE)
  dog_karyogram(hapA[21:39], hapB[21:39], map[21:39], breeds, title, sex, 
                legend = TRUE)

} # dog_karyogram2()



###############
dog_breed_barplot = function(hapA, hapB, breeds, title, colors = "global") {

  # For now, strip out the probs in column 2, we may use them later.
  probA = lapply(hapA, function(z) { z[,2,drop = FALSE] })
  probB = lapply(hapB, function(z) { z[,2,drop = FALSE] })
  hapA = lapply(hapA, function(z) { z[,1,drop = FALSE] })
  hapB = lapply(hapB, function(z) { z[,1,drop = FALSE] })

  confidence.thr  = 0.8
  composition.thr = 0.9

  # Get overall breed composition.
  tmp = get_breed_composition(hapA, probA, hapB, probB, breeds, confidence.thr,
        composition.thr)
  breed = tmp$breed.comp
  num.breeds = tmp$num.breeds

  breed.colors = get_breed_colors(breeds)
  if(is.numeric(colors)) {

    breed.colors = get_breed_colors2(names(breed)[1:num.breeds])

  } # if(is.numeric(colors))

  names(breed.colors) = gsub("_", " ", names(breed.colors))
  names(breed) = gsub("_", " ", names(breed))
  breed = breed * 100

  if(num.breeds == 1) {
    par(plt = c(0.25, 0.99, 0.5, 0.92))
  } else {
    par(plt = c(0.1, 0.99, 0.5, 0.92))
  } # else

  x = barplot(breed[1:num.breeds], col = breed.colors[names(breed)[1:num.breeds]],
          names.arg = NA, las = 2, xlab = "", ylab = "Breed Percentage",
          cex.lab = 1.25, ylim = c(0, max(breed)))
  mtext(side = 3, line = 0, cex = 1.5, text = title)
  usr = par("usr")
  y.offset = diff(usr[3:4]) * 0.02
  text(x = x, y = -y.offset, labels = names(breed)[1:num.breeds], srt = 45,
       adj = c(1,1), xpd = T, cex = 1.5)

} # dog_breed_barplot()

