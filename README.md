# Viral-motif-module
Training module to identify viral motifs within the human genome

## PART I: Which dataset should we be using?
See README.md contained within the Comparing-seqs-by-ID folder for a detailed guide.

## PART II: Creating the motif using your annotated reference sequence as a guide
1. Obtain a reference sequence for your set of data (Recommendation: Use the GenPept full file, which includes mature peptide annotations)
2. Append your reference fasta to your fasta file within Geneious
3. Align the sequences using mafft/other multiple sequence alignment tool
4. Extract 8mers at each site confirmed in literature (4 AAs upstream and 4 AAs downstream)

More TBD next week