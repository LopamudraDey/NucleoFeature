
#sequence <- "ATGCATGCATGCGCGCGTGTGTAAAATTTTT"


curvature_value <- function(sequence_input, k = 2) {
  sequence_input <- toupper(as.character(sequence_input))
  nucleotides <- unlist(strsplit(sequence_input, ""))

  total_kmers <- length(nucleotides) - k + 1
  if (total_kmers < 1) return(NULL)

  kmers <- sapply(1:total_kmers, function(i) paste(nucleotides[i:(i + k - 1)], collapse = ""))


  counts <- table(kmers)

  # curvature values from Olson et al. (1998)
  curvature_values <- c("AA" = 0.127, "AC" = -0.054, "AG" = -0.089, "AT" = 0.183,
                        "CA" = -0.077, "CC" = -0.152, "CG" = 0.058, "CT" = -0.211,
                        "GA" = -0.150, "GC" = 0.129, "GG" = 0.054, "GT" = 0.096,
                        "TA" = -0.100, "TC" = 0.118, "TG" = 0.053, "TT" = -0.097)


  possible_kmers <- expand.grid(rep(list(c("A", "C", "G", "T")), k))
  all_kmers <- apply(possible_kmers, 1, paste, collapse = "")

  full_counts <- setNames(rep(0, length(all_kmers)), all_kmers)
  full_counts[names(counts)] <- counts

  # Normalize
  kmer_freqs <- full_counts / sum(full_counts)

  curvature_score <- sum(kmer_freqs[names(curvature_values)] * curvature_values, na.rm = TRUE)

  return(curvature_score)
}
#sequence <- "ATGCATGCATGCGCATTTACGTGTGTTTGCGCGGGGGG"
#curvature_value <- curvature_value(sequence, k = 2)

#print(paste("DNA Curvature Score:", curvature_value))

