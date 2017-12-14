#CanFamHMM

## Background 

CanFamHmm is an effort to perform haplotype reconstruction in domestic dogs using publicly available data. Several researchers have published genotype data for pure bred dogs and these can serve as refrence data to reconstruct mixed breed dog genomes in terms up pure breed haplotypes. 

The algorithm uses a hidden Markov model to sweep across the genome and call the most probably breed at each locus. I use the reference data to estimate the probability of observing each allele within each breed at each marker. Then, I tune the probability of seeing a recombination between markers. With these two sets of probabilities, I sweep across the genome using the Viterbi algorithm to estimate the most probable breed.

