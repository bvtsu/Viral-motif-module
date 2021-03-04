# Viral-motif-module
Training module to identify viral motifs within the human genome

## Which dataset should we be using?
#### Just by looking at GenBank and ViPR fasta files, we see that all ViPR sequences have an associated GenBank ID.
#### We could assume that ViPR is completely encompassed by GenBank, but let's actually do the comparison.

In your desired github folder:

```git clone https://github.com/bvtsu/Viral-motif-module```

```cd Viral-motif-module```

Let's first do this with the bash script, ```find-diff.sh```.

```chmod +x find-diff.sh```

```./find-diff.sh /path/to/genbank-fasta-file /path/to/vipr-fasta-file```

```cd Examples/output```

```wc -l compared*```

You may find that there are actually some ViPR sequences that don't overlap with the GenBank list. If you individual search up these sequences, you'll likely find that they are still being classified and may be broadly associated with the viral family, rather than the genus we're interested in.

Now let's do this with python...