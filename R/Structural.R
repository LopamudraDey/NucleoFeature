# compute Melting Temperature (Tm) function 
compute_tm <- function(sequence) {
  sequence <- toupper(sequence)  
  counts <- table(strsplit(sequence, "")[[1]])  # Count nucleotides

  A <- ifelse("A" %in% names(counts), counts["A"], 0)
  T <- ifelse("T" %in% names(counts), counts["T"], 0)
  G <- ifelse("G" %in% names(counts), counts["G"], 0)
  C <- ifelse("C" %in% names(counts), counts["C"], 0)
  total <- A + T + G + C

  # Compute Melting Temperature (Tm) using the empirical formula
  if (total > 0) {
    tm_empirical <- 64.9 + 41 * ((G + C) - 16.4) / total
  } else {
    tm_empirical <- NA
  }


  return(list(Tm_Empirical = tm_empirical))
}
}

# compute DNA Stability Index function. Depends on GC content. A high GC content correlates with higher stability.
compute_stability <- function(sequence) {
  sequence <- toupper(sequence)
  counts <- table(strsplit(sequence, "")[[1]])

  G <- ifelse("G" %in% names(counts), counts["G"], 0)
  C <- ifelse("C" %in% names(counts), counts["C"], 0)
  total <- sum(counts)

  stability_index <- ifelse(total > 0, (G + C) / total, NA)  # GC content as proxy for stability
  return(stability_index)
}

# Function to estimate Z-DNA Probability. I have used 4 Z-dna motifs
compute_z_dna_probability <- function(sequence) {
  sequence <- toupper(sequence)  # Convert to uppercase

  # Known Z-DNA motifs (alternating purines and pyrimidines)
  z_dna_motifs <- c("CGCGCG", "GCGCGC", "TATA", "ATAT")

  # Count occurrences of each motif
  motif_counts <- sapply(z_dna_motifs, function(motif) {
    matches <- gregexpr(motif, sequence)[[1]]
    sum(matches > 0)  # Count valid matches
  })

  # Compute probability as a fraction of sequence length
  total_motif_count <- sum(motif_counts)
  z_dna_prob <- total_motif_count / nchar(sequence)

  return(z_dna_prob)
}

