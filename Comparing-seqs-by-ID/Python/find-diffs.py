#!/usr/bin/env python3
import os #for os.mkdir
import sys #for sys.argv
from typing import Dict #to give type specs
from collections import defaultdict #autofills values when values aren't found for a given key

#Define function that will read fasta files and store header, sequence pairs in dict format
def fasta_to_dict(filename: str, fasta_type: str) -> Dict[str, str]:
    fasta = {}
    with open(filename) as file: #open the fasta formatted file
        for line in file: #go line by line
            line = line.strip() #take away the whitespace
            if not line: #if it's a blank line,
                continue #go to the next line
            if line.startswith(">"): #these ">" lines are the header info
                if fasta_type=="genbank": #if this is genbank file
                    seq_header = line[1:].split('.')[0] #start from the second character, and save up to the first .
                elif fasta_type=="vipr": #if this is a vipr file
                    seq_header = line[4:].split('|')[0] #start from the 5th character, and save til the first |
                if seq_header not in fasta:
                    fasta[seq_header] = []
                continue #go to the next line, which should be the fasta seq
            seq = line #store the entire line as a string variable
            fasta[seq_header].append(seq) #associate the header (key) with the seq (value)
    for key, value in fasta.items(): #take the key, value pair
        fasta[key]=''.join(value) #overwrite the existing list of broken up fasta seqs with a concatenated form, easier to write to file later
    return(fasta) #Return the fasta dict

#Probably didn't need to make a function for this, since we're only using this once
def return_dict_concat(dict1,dict2,dict3) -> Dict[str, str]:
    combined_dict=[]
    combined_dict=[dict1, dict2, dict3] #Store all the dicts as a list
    final_dict = defaultdict(set) #defaultdict(set) to prepare a dict without duplicates
    for d in combined_dict: #iterate through each of the dicts
        for key, value in d.items(): #take the key, value pair
            final_dict[key].add(value) #relocate key, value pairs from each dict to the new, combined dict
    return(final_dict)

#Make output folder, if it doesn't already exist
try:
    os.mkdir("outputs/")
except OSError:
    print("Creation of the directory %s failed, already exists")
else:
    print("Successfully created the output directory.")

#Set file paths using the first and second terms following the python script name
genbank_fasta=sys.argv[1] #Take the first term (should be genbank) that follows the python script name
vipr_fasta=sys.argv[2] #Take the second term (should be vipr) that follows the python script name

genbank_dict=fasta_to_dict(genbank_fasta, "genbank") #Use the function defined above to store the genbank fasta as a dict
vipr_dict=fasta_to_dict(vipr_fasta, "vipr") #Use the function defined above to store the vipr fasta as a dict

genbank_only = { key : genbank_dict[key] for key in set(genbank_dict) - set(vipr_dict) } #subtract vipr pairs from genbank pairs, then use the key from the pairs to regenerate a dict of unique genbank seqs
vipr_only = { key : vipr_dict[key] for key in set(vipr_dict) - set(genbank_dict) } #subtract genbank pairs from vipr pairs, then use the key from the pairs to regenerate a dict of unique vipr seqs
shared = { key : vipr_dict[key] for key in set(vipr_dict) & set(genbank_dict) } #take only the pairs that exist in both dicts, use the key from the pairs to regenerate a dict of shared seqs
print("# of unique entries for genbank: ",len(genbank_only)) #Print total uniques for genbank
print("# of unique entries for vipr: ",len(vipr_only)) #Print total uniques for vipr
print("# of shared entries: ",len(shared)) #Print total shared

try:
    if sys.argv[3]=="all": #if you add "all" as a third term following the python script name
        all_dicts = return_dict_concat(genbank_only,vipr_only,shared) #use the function defined above to concat all the dicts
        print("# of combined entries returned".format(),len(all_dicts)) #Print total combined (differences from either end and shared)
        with open("outputs/combined.fasta", "w") as output_file: #create a new fasta file in the outputs folder
            for key, value in all_dicts.items(): #Take each key, value pair
                #print('{0} corresponds to {1}'.format(key, list(value)[0]))
                output_file.write(">" + key + "\n" +list(value)[0] + "\n\n") #Store them into the new fasta file, separated by newlines
except: #if you don't add "all"
    print("No combined fasta requested.") #Report that no new fasta file was created