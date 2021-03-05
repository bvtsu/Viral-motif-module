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

def ID_keys_to_text(in_dict: Dict[str, str], output_name: str):
    with open("outputs/{}".format(output_name), "w") as f:
        for key in in_dict:
            f.write("{}\n".format(key))

def return_dict_concat_fasta(in_dict1: Dict[str, str],in_dict2: Dict[str, str],in_dict3: Dict[str, str]) -> Dict[str, str]:
    combined_dict=[]
    combined_dict=[in_dict1, in_dict2, in_dict3] #Store all the dicts as a list
    final_dict = defaultdict(set) #defaultdict(set) to prepare a dict without duplicates
    with open("outputs/combined.fasta", "w") as output_file: #create a new fasta file in the outputs folder
        for d in combined_dict: #iterate through each of the dicts
            for key, value in d.items(): #take the key, value pair
                final_dict[key].add(value) #relocate key, value pairs from each dict to the new, combined dict
                output_file.write(">{0}\n{1}\n\n".format(key,value)) #Store them into the new fasta file, separated by newlines
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

#I should've made the first four lines of this section into a function
genbank_only = { key : genbank_dict[key] for key in set(genbank_dict) - set(vipr_dict) } #subtract vipr pairs from genbank pairs, then use the key from the pairs to regenerate a dict of unique genbank seqs
vipr_only = { key : vipr_dict[key] for key in set(vipr_dict) - set(genbank_dict) } #subtract genbank pairs from vipr pairs, then use the key from the pairs to regenerate a dict of unique vipr seqs
shared = { key : vipr_dict[key] for key in set(vipr_dict) & set(genbank_dict) } #take only the pairs that exist in both dicts, use the key from the pairs to regenerate a dict of shared seqs
shared_by_fasta = {key: vipr_dict[key] for key in vipr_dict if key in genbank_dict and vipr_dict[key] == genbank_dict[key]} #Didn't finish all of this, but one could also compare fasta values to remove dupes/see similarities between lists
print("# of unique entries for genbank: ",len(genbank_only)) #Print total uniques for genbank
print("# of unique entries for vipr: ",len(vipr_only)) #Print total uniques for vipr
print("# of shared entries: ",len(shared)) #Print total shared
print("# of shared entries by fasta: ",len(shared_by_fasta)) #Print total shared
ID_keys_to_text(genbank_only,"genbank_only.txt")
ID_keys_to_text(vipr_only,"vipr_only.txt")
ID_keys_to_text(shared,"shared.txt")


try:
    if sys.argv[3]=="all": #if you add "all" as a third term following the python script name
        all_dicts = return_dict_concat_fasta(genbank_only,vipr_only,shared) #use the function defined above to concat all the dicts
        print("# of combined entries returned".format(),len(all_dicts)) #Print total combined (differences from either end and shared)
except: #if you don't add "all"
    print("No combined fasta requested.") #Report that no new fasta file was created