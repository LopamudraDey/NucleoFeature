shannon_entropy <- function(sequence_in) {
  sequence_in <- as.character(sequence_in)
  sequence_in <- toupper(unlist(strsplit(sequence_in, "")))
  freq <- table(sequence_in) / length(sequence_in)  # Probability of each nucleotide
  entropy <- -sum(freq * log2(freq))  # Shannon entropy formula
  return(entropy)
}


