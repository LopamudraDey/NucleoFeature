pseknc_composition <- function(sequence_input, k, lambda = 2, w = 0.05) {
  # Ensure input is uppercase
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
  
  # Normalize frequencies
  kmer_freqs <- full_counts / sum(full_counts)
  
  # Compute sequence-order correlation
  tau_values <- numeric(lambda)
  for (i in 1:lambda) {
    paired_kmers <- sapply(1:(total_kmers - i), function(j) {
      paste(nucleotides[j:(j + k - 1)], collapse = "")
    })
    paired_kmers_next <- sapply((1 + i):(total_kmers), function(j) {
      paste(nucleotides[j:(j + k - 1)], collapse = "")
    })
    
    tau_values[i] <- sum(kmer_freqs[paired_kmers] * kmer_freqs[paired_kmers_next])
  }
  
  # Compute PseKNC vector
  pseknc_vector <- c(kmer_freqs, w * tau_values)
  pseknc_df <- data.frame(Feature = c(names(kmer_freqs), paste0("Tau_", 1:lambda)), 
                          Value = as.numeric(pseknc_vector))
  
  return(pseknc_df)
}
