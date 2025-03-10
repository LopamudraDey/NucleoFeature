
helical_twist <- function(dna_sequence) {

  bp_count <- nchar(dna_sequence)

  # assuming B-DNA has 10.5 base pairs per turn
  bp_per_turn <- 10.5

  turns <- bp_count / bp_per_turn

  twist_angle <- 360 * turns

  return(twist_angle)
}


#dna_sequence <- "ATGCATGCATGCATGCATGCATGCATGCATGC"


# helical_twist <- helical_twist(dna_sequence)
# print(paste("Helical Twist (degrees):", helical_twist))
