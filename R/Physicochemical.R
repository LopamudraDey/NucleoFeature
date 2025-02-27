# Define physicochemical property values for dinucleotides
dinucleotide_properties <- data.frame(
  Dinucleotide = c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT",
                   "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"),
  Hydrophobicity = c(1.0, 0.8, 1.2, 0.6, 0.8, 0.5, 1.4, 0.7,
                     1.2, 1.5, 0.9, 0.7, 0.6, 0.7, 0.9, 1.0),
  Stacking_Energy = c(-0.5, -0.8, -1.2, -0.3, -0.8, -1.1, -2.0, -0.7,
                      -1.2, -1.9, -1.4, -0.6, -0.3, -0.7, -0.9, -0.5),
  Bendability = c(0.8, 0.7, 0.9, 0.6, 0.7, 0.5, 1.3, 0.6,
                  0.9, 1.2, 0.8, 0.6, 0.6, 0.6, 0.7, 0.8),
  Flexibility = c(0.7, 0.6, 0.8, 0.5, 0.6, 0.4, 1.1, 0.5,
                  0.8, 1.0, 0.7, 0.5, 0.5, 0.5, 0.6, 0.7)
)

# Function to calculate physicochemical properties of a DNA sequence
compute_physicochemical_properties <- function(sequence) {
  sequence <- toupper(sequence)
  seq_length <- nchar(sequence)
  dinucleotides <- substring(sequence, 1:(seq_length-1), 2:seq_length)

  # Initialize property sums
  hydrophobicity_sum <- 0
  stacking_energy_sum <- 0
  bendability_sum <- 0
  flexibility_sum <- 0

  # Compute values for each dinucleotide
  for (di in dinucleotides) {
    row <- subset(dinucleotide_properties, Dinucleotide == di)
    if (nrow(row) > 0) {
      hydrophobicity_sum <- hydrophobicity_sum + row$Hydrophobicity
      stacking_energy_sum <- stacking_energy_sum + row$Stacking_Energy
      bendability_sum <- bendability_sum + row$Bendability
      flexibility_sum <- flexibility_sum + row$Flexibility
    }
  }

  # Normalize values per dinucleotide
  dinucleotide_count <- length(dinucleotides)
  if (dinucleotide_count > 0) {
    hydrophobicity_avg <- hydrophobicity_sum / dinucleotide_count
    stacking_energy_avg <- stacking_energy_sum / dinucleotide_count
    bendability_avg <- bendability_sum / dinucleotide_count
    flexibility_avg <- flexibility_sum / dinucleotide_count
  } else {
    hydrophobicity_avg <- NA
    stacking_energy_avg <- NA
    bendability_avg <- NA
    flexibility_avg <- NA
  }

  return(list(
    Hydrophobicity = hydrophobicity_avg,
    Stacking_Energy = stacking_energy_avg,
    Bendability = bendability_avg,
    Flexibility = flexibility_avg
  ))
}

# Function to estimate nucleosome positioning signals
compute_nucleosome_positioning <- function(sequence) {
  sequence <- toupper(sequence)
  nucleosome_motif <- "WWSSWW"  # W = A/T, S = G/C
  count <- sum(grepl(nucleosome_motif, sequence))
  return(count / nchar(sequence))  # Normalize by sequence length
}

# Function to estimate DNA bending stiffness (simplified as GC content)
compute_bending_stiffness <- function(sequence) {
  sequence <- toupper(sequence)
  g_count <- sum(strsplit(sequence, "")[[1]] == "G")
  c_count <- sum(strsplit(sequence, "")[[1]] == "C")
  return((g_count + c_count) / nchar(sequence))  # GC content as stiffness measure
}
