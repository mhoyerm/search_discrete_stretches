Search discrete stretches

Searches in a fasta file of ORFs or proteins, stretches of nucleotides or peptides that were selected in lists and returns the protein or genome names with the respective found sequence.

Usage:

python search_discrete_stretches.py -F ORFs.fasta -N nucleotides.txt -E peptides.txt -O output.csv --ignoreframe

python search_discrete_stretches.py -F proteome.fasta -E peptides.txt -O output.csv --protein

-h [optional]: shows helpful information

-F [required]: fasta archive containing ORFs or proteins

--protein [optional]: use this command if the fasta archive is a protein file

--keepframe [optional]: use this command if you want to respect the reading frame

-E [optional]: list in txt or csv format containing wanted peptides

-N [optional]: list in txt or csv format containing wanted nucletide sequence

-O [required]: output file in csv format

