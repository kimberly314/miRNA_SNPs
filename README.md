miRNA_SNPs Repository

Computational Genetics Script(s) for Analyzing Variation in microRNA

OVERVIEW:

mir_sequences.py takes in a data file with microRNA SNPs (variations) and a file with microRNA sequences and outputs both microRNA sequence variations for folding analysis.

mir_2_sequences is a continuation of mir_sequences with more functionality, such as outputting reverse complements.

BACKGROUND:

This was one of my first Python projects. It was created as part of research funded by my NASA Space Grant and Merck Fellowship at the University of Washington in the Akey Lab. This portion of our research was focused on determining the effect of variations in some of the smallest functional units of the human genome, microRNA.

Our earlier research found a lower number of variations in microRNA as compared to adjacent regions. This indicates these genes are highly conserved. Folding analysis can show us how variations affect the molecular shape and function.

USAGE:

Healthy microRNA are a hairpin shape. By using mir_sequences.py to processes data from our earlier analysis, we create a file containing dual sequences for each of the microRNA of interest. One sequence contains the first allele and one contains the second. These sequences can then be put into modeling programs to observe the effect of alleles.
