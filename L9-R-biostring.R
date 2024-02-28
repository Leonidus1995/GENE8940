# load package into your current environment
library(Biostrings)

# read sequence
dna <- readDNAStringSet("/home/fg69001/GENE8940/L9-R-test.fasta")

# sequence transformation
dna_complement <- complement(dna)
dna_reversed = reverse(dna)
dna_reverse_complement = reverseComplement(dna)

# write new sequence as fasta file
writeXStringSet(dna_reverse_complement, "/work/gene8940/fg69001/test_reverse_complement.fasta")
