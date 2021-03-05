#!/bin/bash
#start with a shebang line to indicate that Bash from this absolute path will be used to interpret the following code below

#This script compares a genbank list of fasta seqs with a vipr list of fasta seqs
#Optional: One could create flags to specify -genbank or -vipr to make the script re-usable for different combinations of format types

#$1 = genbank 
#$2 = vipr
fname=`basename $1` #Get only the file name from the first provided term, the genbank file path
fname2=`basename $2` #Get only the file name from the second provided term, the vipr file path
seq_dir=`dirname $1` #Get directory of where your seqs are stored, assuming they're in the same folder

if [ -d ${seq_dir}/output ] #Check if directory exists
then
    echo "Output directory exists." #If it does, print this text
else
    echo "Creating output directory." #If it doesn't, print this text
    mkdir ${seq_dir}/output #and make the directory
fi

#grep the genbank file for only header lines using the character ">",
#split every line by the ">" and take the second segment,
#split every line by the "." and take the first segment
#save the resulting ID list output
grep ">" ${1} | cut -d ">" -f2 | cut -d "." -f1 > ${seq_dir}/output/${fname}.IDs.out 

#grep the vipr file for only header lines using the character ">",
#split every line by the "|" and take the first segment,
#split every line by the ":" and take the second segment
#save the resulting ID list output
grep ">" ${2} | cut -d "|" -f1 | cut -d ":" -f2 > ${seq_dir}/output/${fname2}.IDs.out

#comm -12 produces a list of shared elements between the genbank and vipr ID lists
#comm -23 produces a list of elements UNIQUE to the genbank file
#comm -13 produces a list of elements UNIQUE to the vipr file
comm -12 <(sort ${seq_dir}/output/${fname}.IDs.out) <(sort ${seq_dir}/output/${fname2}.IDs.out) | sort > ${seq_dir}/output/compared-${fname}-${fname2}-shared.IDs.out
comm -23 <(sort ${seq_dir}/output/${fname}.IDs.out) <(sort ${seq_dir}/output/${fname2}.IDs.out) | sort > ${seq_dir}/output/compared-${fname}-only.IDs.out
comm -13 <(sort ${seq_dir}/output/${fname}.IDs.out) <(sort ${seq_dir}/output/${fname2}.IDs.out) | sort > ${seq_dir}/output/compared-${fname2}-only.IDs.out