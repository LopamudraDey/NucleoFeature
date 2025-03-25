reverse_seqcomplement <- function(sequence) {
  # we are replacing A with T, T with A, C with G, G with C
  complement <- c(A = "T", T = "A", C = "G", G = "C")
  sequence <- unlist(strsplit(sequence, split = ""))
  rev_complement_sequence <- complement[sequence]
  return(paste(rev(rev_complement_sequence), collapse = ""))

}
kmer_compositionreverse <- function(sequence, k) {
  kmers <- substring(sequence, 1:(nchar(sequence) - k + 1), k:nchar(sequence))
  kmer_table <- table(kmers)  # Counts the frequency of each k-mer
  return(kmer_table)
}

sequence_input <- "AGCTAGCTAG"

# Compute reverse complement of the sequence
reverse_complement_sequence <- reverse_complement(sequence_input)
cat("Reverse Complement Sequence:", reverse_complement_sequence, "\n")
