#!/bin/bash
#start with a shebang line to indicate that Bash from this absolute path will be used to interpret the following code below

#This script compares a genbank list of fasta seqs with a vipr list of fasta seqs
#Optional: One could create flags to specify -genbank or -vipr to make the script re-usable for different combinations of format types

#$1 = genbank 
#$2 = vipr
#dirname=`dirname $1`
#dirname2=`dirname $2`
fname=`basename $1`
fname2=`basename $2`

if [ -d "Examples/output/" ] 
then
    echo "Output directory exists." 
else
    echo "Creating output directory."
    mkdir Examples/output/
fi
grep ">" ${1} | cut -d ">" -f2 | cut -d "." -f1 > Examples/output/${fname}.IDs.out

grep ">" ${2} | cut -d "|" -f1 | cut -d ":" -f2 > Examples/output/${fname2}.IDs.out

comm -12 <(sort Examples/output/${fname}.IDs.out) <(sort Examples/output/${fname2}.IDs.out) | sort > Examples/output/compared-${fname}-${fname2}-shared.IDs.out
comm -23 <(sort Examples/output/${fname}.IDs.out) <(sort Examples/output/${fname2}.IDs.out) | sort > Examples/output/compared-${fname}-only.IDs.out
comm -13 <(sort Examples/output/${fname}.IDs.out) <(sort Examples/output/${fname2}.IDs.out) | sort > Examples/output/compared-${fname2}-only.IDs.out