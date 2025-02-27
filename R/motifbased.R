detect_tfbs <- function(sequence, motifs = c(
  "TATAAA", "CACGTG", "GATA", "GGGCGG", "TGGGGA", "TTTAAA", "CCGCCC", "AGGAGG",
  "ATGCAAAT", "CCATGG", "CGCGCG", "AGCT", "CTCF", "GGGACTTTCC", "TGAGTCA",
  "AATTAA", "GGAGGA", "CGTACG", "GAAGGA", "TTAGGG","GTCCCCAGGGGA","GTCCCCTGGGGA", "AAACCACAA",
  "AAACCACAC", # RUNX1 binding variant (R = A, M = C)
  "AAACCACGA", # RUNX1 binding variant (R = G, M = A)
  "AAACCACGC"),
  tf_names = c(
    "TATA Box", "E-box (Myc)", "GATA Factors", "SP1", "NF-kB", "Octamer", "AP-2",
    "PU.1", "PAX6", "CCAAT-box (NF-Y)", "CpG Island", "AP-1", "CTCF", "NFAT",
    "Jun-Fos (AP-1)", "Polyadenylation", "ETS Family", "Forkhead (FOXO)",
    "RUNX", "Telomeric Repeat (TRF)","Ebf1","Ebf1","RUNX1","RUNX1","RUNX1","RUNX1" ))
{
  matches <- sapply(motifs, function(motif) {
    length(gregexpr(motif, sequence, ignore.case = TRUE)[[1]]) - 1
  })

  # Create a data frame with TF names and frequencies
  result_df <- data.frame(
    BindingSite = tf_names,
    Motif = motifs,
    Frequency = matches,
    stringsAsFactors = FALSE
  )

  return(result_df)
}


# Function to identify CpG Islands
detect_cpg_island <- function(seq, min_length = 200, gc_threshold = 0.5, oe_threshold = 0.6) {
  seq <- toupper(seq)  # Convert to uppercase to ensure consistency

  # Compute GC content
  gc_count <- sum(strsplit(seq, "")[[1]] %in% c("G", "C"))
  gc_content <- gc_count / nchar(seq)  # Calculate GC content as a fraction

  # Count CpG dinucleotides
  cpg_count <- sum(gregexpr("CG", seq, fixed = TRUE)[[1]] > 0)

  # Count occurrences of 'G' and 'C' separately
  g_count <- sum(strsplit(seq, "")[[1]] == "G")
  c_count <- sum(strsplit(seq, "")[[1]] == "C")

  # Compute observed/expected CpG ratio
  expected_cpg <- (g_count * c_count) / nchar(seq)
  oe_ratio <- ifelse(expected_cpg > 0, cpg_count / expected_cpg, 0)

  # Check if the sequence qualifies as a CpG island
  if (nchar(seq) >= min_length && gc_content >= gc_threshold && oe_ratio >= oe_threshold) {
    return(list(CpG_Island = TRUE, CpG_Count = cpg_count, GC_Content = gc_content, OE_Ratio = oe_ratio))
  } else {
    return(list(CpG_Island = FALSE, CpG_Count = cpg_count, GC_Content = gc_content, OE_Ratio = oe_ratio))
  }
}

# 3️⃣ Function to approximate conserved elements using sequence similarity
detect_conserved_elements <- function(sequence, reference_sequence) {
  alignment <- pairwiseAlignment(sequence, reference_sequence, type = "global")
  score <- pid(alignment)  # Percentage identity score
  return(score)
}
