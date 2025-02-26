gc_content <- function(sequence_in) {
  sequence_in <- as.character(sequence_in)
  sequence_in <- toupper(unlist(strsplit(sequence_in, "")))
  gc <- sum(sequence_in == "G" | sequence_in == "C") / length(sequence_in) * 100
  return(gc)
}


