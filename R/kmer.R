#kmer_frequency <- function(sequence_input, k) {
  # Ensure input is character
  #sequence_input <- toupper(as.character(sequence_input))
 # nucleotides <- unlist(strsplit(sequence_input, ""))

  # Generate k-mers
 # total_kmers <- length(nucleotides) - k + 1
 # if (total_kmers < 1) return(NULL)

#  kmers <- sapply(1:total_kmers, function(i) paste(nucleotides[i:(i + k - 1)], collapse = ""))

  # Count occurrences
#  counts <- table(kmers)

  # Get all possible k-mers
#  possible_kmers <- expand.grid(rep(list(c("A", "C", "G", "T")), k))
#  all_kmers <- apply(possible_kmers, 1, paste, collapse = "")

  # Fill missing k-mers with zero
 # full_counts <- setNames(rep(0, length(all_kmers)), all_kmers)
#  full_counts[names(counts)] <- counts

  # Normalize
 # freq <- full_counts / sum(full_counts) * 100
#  kmer_df <- data.frame(Kmer = names(kmer_freqs), Frequency = as.numeric(kmer_freqs))

#  return(kmer_df)
#  return(freq)
#}

kmer_frequency <- function(seq, k) {
    # Convert to uppercase and remove unwanted characters
  seq <- toupper(as.character(seq))
    seq <- gsub("[^ACGT]", "", seq)  # Remove invalid characters

    # Get all k-mers from the sequence
    kmers <- substring(seq, 1:(nchar(seq) - k + 1), k:nchar(seq))

    # Generate ALL possible k-mers of length k
    possible_kmers <- unique(outer(c("A", "C", "G", "T"), rep("", k - 1), paste0))
    for (i in 2:k) {
        possible_kmers <- unique(outer(possible_kmers, c("A", "C", "G", "T"), paste0))
    }

    # Count occurrences of each k-mer
    kmer_counts <- table(factor(kmers, levels = possible_kmers))

    # Normalize frequencies (percentage)
    kmer_freqs <- (kmer_counts / sum(kmer_counts)) * 100

    # Convert to a data frame with column names
    kmer_df <- data.frame(Kmer = names(kmer_freqs), Frequency = as.numeric(kmer_freqs))

    return(kmer_df)
}



