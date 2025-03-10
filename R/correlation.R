compute_autocorrelation <- function(dna_seq, max_lag = 3) {
  # Define nucleotide EIIP values
  nucleotide_properties <- list(
    "A" = 0.1260, "T" = 0.1335,
    "G" = 0.0806, "C" = 0.1340
  )


  nucleotides <- unlist(strsplit(dna_seq, ""))
  numerical_seq <- sapply(nucleotides, function(n) ifelse(!is.null(nucleotide_properties[[n]]), nucleotide_properties[[n]], NA))

  numerical_seq <- na.omit(numerical_seq)
  n <- length(numerical_seq)

  ac_features <- numeric(max_lag)
  morans_features <- numeric(max_lag)
  gearys_features <- numeric(max_lag)

  # Compute autocorrelation features
  mean_val <- mean(numerical_seq, na.rm = TRUE)
  total_variance <- sum((numerical_seq - mean_val)^2, na.rm = TRUE)

  for (lag in 1:max_lag) {
    if (n > lag) {
      ac_features[lag] <- sum((numerical_seq[1:(n - lag)] - mean_val) * (numerical_seq[(1 + lag):n] - mean_val)) / (n - lag)
      morans_features[lag] <- ((n - lag) * sum((numerical_seq[1:(n - lag)] - mean_val) * (numerical_seq[(1 + lag):n] - mean_val))) / (total_variance * (n - 1))
      gearys_features[lag] <- sum((numerical_seq[1:(n - lag)] - numerical_seq[(1 + lag):n])^2) / (2 * total_variance)
    } else {
      ac_features[lag] <- NA
      morans_features[lag] <- NA
      gearys_features[lag] <- NA
    }
  }

  return(list(
    Auto_Covariance = ac_features,
    Moran_I = morans_features,
    Geary_C = gearys_features
  ))
}


