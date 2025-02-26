kmer_frequency <- function(sequence_input, k) {
  # Ensure input is character
  sequence_input <- toupper(as.character(sequence_input))
  nucleotides <- unlist(strsplit(sequence_input, ""))

  # Generate k-mers
  total_kmers <- length(nucleotides) - k + 1
  if (total_kmers < 1) return(NULL)

  kmers <- sapply(1:total_kmers, function(i) paste(nucleotides[i:(i + k - 1)], collapse = ""))


  counts <- table(kmers)

  # Get all possible k-mers
  possible_kmers <- expand.grid(rep(list(c("A", "C", "G", "T")), k))
  all_kmers <- apply(possible_kmers, 1, paste, collapse = "")

  # Fill missing k-mers with zero
  full_counts <- setNames(rep(0, length(all_kmers)), all_kmers)
  full_counts[names(counts)] <- counts

  # Normalize
  kmer_freqs <- full_counts / sum(full_counts) * 100
  kmer_df <- data.frame(Kmer = names(kmer_freqs), Frequency = as.numeric(kmer_freqs))

 return(kmer_df)

}




