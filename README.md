##### Nucleotide Feature Generation###############
The package contains 4 functions ACTG_composition(), frequency(), gc_content() and shannon_entropy().
# How to use this file?
open R
install.packages("remotes")
# Install your package from GitHub
remotes::install_github("LopamudraDey/NucleoFeature")
library(NucleoFeature)
# Example 
result <- NucleoFeature::ACTG_composition("AGCTAGCTAGCCTTTT")
print(result)

nucleotides
    A     C     G     T 
18.75 25.00 18.75 37.50

# Example 2
frequency()calculates K-mers.
It takes 2 inputs: sequence and size of K
result <- NucleoFeature::frequency("AGCTAGCTAGCCTTTT", 2)
   Kmer Frequency
1    AA  0.000000
2    CA  0.000000
3    GA  0.000000
4    TA 13.333333
5    AC  0.000000
6    CC  6.666667
7    GC 20.000000
8    TC  0.000000
9    AG 20.000000
10   CG  0.000000
11   GG  0.000000
12   TG  0.000000
13   AT  0.000000
14   CT 20.000000
15   GT  0.000000
16   TT 20.000000


# Example 3
gc_contect() calculates % of G+C in the sequence 
seq="ACTCTCAGCT"
gc_content(seq)
50
