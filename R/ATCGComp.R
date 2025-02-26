
ACTG_composition <- function(sequence_input) {
  # Check if input is a character
  if (!is.character(sequence_input)) {
    stop("Error: Input must be a character string")
  }

  # Convert to uppercase and split into individual nucleotides
  sequence_input <- toupper(sequence_input)
  nucleotides <- unlist(strsplit(sequence_input, ""))

  # Calculate frequency of each nucleotide
  counts <- table(nucleotides)
  composition <- counts / sum(counts) * 100  # Convert to percentage

  return(composition)
}

