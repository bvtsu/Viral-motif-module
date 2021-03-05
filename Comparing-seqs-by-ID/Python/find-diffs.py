#!/usr/bin/env python3
import sys #We need this to use sys.argv
from typing import Dict

def fasta_to_dict(filename: str) -> Dict[str, str]:
    fasta = {}
    with open(filename) as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                seq_header = line[1:].split('.')[0]
                if seq_header not in fasta:
                    fasta[seq_header] = []
                continue
            seq = line
            fasta[seq_header].append(seq)
    for key, value in fasta.items():
        fasta[key]=''.join(value)
    return(fasta)

#Set file paths using the first and second terms following the python script name
genbank_fasta=sys.argv[1]
vipr_fasta=sys.argv[2]

genbank_dict=fasta_to_dict(genbank_fasta)
vipr_dict=fasta_to_dict(vipr_fasta)