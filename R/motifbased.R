detect_tfbs <- function(sequence, motifs = c(
  "TATAAA",   # TATA box (Core promoter)
  "CACGTG",   # E-box (Myc binding site)
  "GATA",     # GATA transcription factors
  "GGGCGG",   # SP1 transcription factor
  "TGGGGA",   # NF-kB binding site
  "TTTAAA",   # Octamer binding site
  "CCGCCC",   # AP-2 binding site
  "AGGAGG",   # PU.1 binding site
  "ATGCAAAT", # PAX6 binding site
  "CCATGG",   # CCAAT-box (NF-Y binding)
  "CGCGCG",   # CpG island recognition
  "AGCT",     # AP-1 binding site
  "CTCF",     # CTCF insulator binding
  "GGGACTTTCC", # NFAT binding site
  "TGAGTCA",  # Jun-Fos (AP-1) binding
  "AATTAA",   # Polyadenylation signal
  "GGAGGA",   # ETS family binding
  "CGTACG",   # Forkhead (FOXO) binding
  "GAAGGA",   # RUNX binding site
  "TTAGGG"    # Telomeric repeat (TRF binding)
)) {
  matches <- sapply(motifs, function(motif) {
    length(gregexpr(motif, sequence, ignore.case = TRUE)[[1]]) - 1
  })
  return(matches)
}

# 2️⃣ Function to identify CpG Islands
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
