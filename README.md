# Viral-motif-module
Training module to identify viral motifs within the human genome

## PART I: Which dataset should we be using?
First step is to download the sequences of interest from either genbank or vipr. I'll add more notes on this later, but the critical example details include using the protein databases of either genbank or vipr and filtering by polyproteins (size ~1800-4000 AAs, or 3x as much for the CDS, remove dupes)

Just by looking at GenBank and ViPR fasta files, we see that all ViPR sequences have an associated GenBank ID.
We could assume that ViPR is completely encompassed by GenBank, but let's actually do the comparison.

In your desired github folder:

```git clone https://github.com/bvtsu/Viral-motif-module```

```cd Viral-motif-module/Comparing-seqs-by-ID/```

Let's first do this with the bash script, ```find-diff.sh```.

```cd Bash/```

```chmod +x find-diff.sh```

```./find-diff.sh /path/to/genbank-fasta-file /path/to/vipr-fasta-file```

```cd ../Examples/output```

```wc -l compared*```

You may find that there are actually some ViPR sequences that don't overlap with the GenBank list. If you individual search up these sequences, you'll likely find that they are still being classified and may be broadly associated with the viral family, rather than the genus we're interested in.

Now let's do this with python. Although addressing the question takes a bit more code, storing the fasta information is way more straightforward. Assuming you're in the output folder:

```cd ../../Python```
 
```python find-diffs.py /path/to/genbank-fasta-file /path/to/vipr-fasta-file```

This reports in your terminal the exact same values you saw in the ```wc -l compared*``` line.
As a bonus, you can also request to create a new fasta file containing the merged fasta list from both datasets by adding the "all" tag at the end of the command.

```python find-diffs.py /path/to/genbank-fasta-file /path/to/vipr-fasta-file all```

As an added bonus, if you have streamlit installed:

```streamlit run id_upsetplot.py outputs/genbank_only.txt outputs/vipr_only.txt outputs/shared.txt "Title of your plot"```

<img src="https://github.com/daugherty-lab/Viral-motif-module/blob/main/img/genbank_vs_vipr.png?raw=true" width="50%">

## PART II: Creating the motif using your annotated reference sequence as a guide

TBD