#!/usr/bin/env python3
import os
import sys #We need this to use sys.argv
from typing import Dict
from collections import defaultdict

def fasta_to_dict(filename: str, fasta_type: str) -> Dict[str, str]:
    fasta = {}
    with open(filename) as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if fasta_type=="genbank":
                    seq_header = line[1:].split('.')[0]
                elif fasta_type=="vipr":
                    seq_header = line[4:].split('|')[0]
                if seq_header not in fasta:
                    fasta[seq_header] = []
                continue
            seq = line
            fasta[seq_header].append(seq)
    for key, value in fasta.items():
        fasta[key]=''.join(value)
    return(fasta)

def return_fasta_concat(dict1,dict2,dict3) -> Dict[str, str]:
    combined_dict=[]
    combined_dict=[dict1, dict2, dict3]
    final_dict = defaultdict(set) # defaultdict(set) to avoid duplicates
    for d in combined_dict:
        for key, value in d.items():  # use d.iteritems() in python 2
            final_dict[key].add(value)
    return(final_dict)

try:
    os.mkdir("outputs/")
except OSError:
    print("Creation of the directory %s failed, already exists")
else:
    print("Successfully created the directory.")

#Set file paths using the first and second terms following the python script name
genbank_fasta=sys.argv[1]
vipr_fasta=sys.argv[2]

genbank_dict=fasta_to_dict(genbank_fasta, "genbank")
vipr_dict=fasta_to_dict(vipr_fasta, "vipr")

genbank_only = { key : genbank_dict[key] for key in set(genbank_dict) - set(vipr_dict) }
vipr_only = { key : vipr_dict[key] for key in set(vipr_dict) - set(genbank_dict) }
vipr_only = { key : vipr_dict[key] for key in set(vipr_dict) - set(genbank_dict) }
shared = { key : vipr_dict[key] for key in set(vipr_dict) & set(genbank_dict) }
print("# of unique entries for genbank: ",len(genbank_only))
print("# of unique entries for vipr: ",len(vipr_only))
print("# of shared entries: ",len(shared))

if sys.argv[3]=="all":
    all_dicts = return_fasta_concat(genbank_only,vipr_only,shared)
    print("# of combined entries returned".format(),len(all_dicts))
    with open("outputs/combined.fasta", "w") as output_file:
        for key, value in all_dicts.items():
            #print('{0} corresponds to {1}'.format(key, list(value)[0]))
            output_file.write(">" + key + "\n" +list(value)[0] + "\n\n")

#do not forget to close it
