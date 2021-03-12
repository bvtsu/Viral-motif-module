#!/usr/bin/env python3

import os #for os.mkdir
import sys #for sys.argv
from typing import Dict #to give type specs

#Define function that will read fasta files and store header, sequence pairs in dict format
def fasta_to_dict(filename: str) -> Dict[str, str]:
    fasta = {}
    with open(filename) as file: #open the fasta formatted file
        for line in file: #go line by line
            line = line.strip() #take away the whitespace
            if not line: #if it's a blank line,
                continue #go to the next line
            if line.startswith(">"): #these ">" lines are the header info
                seq_header = line[1:] #start from the second character, and save up to the first .
                if seq_header not in fasta:
                    fasta[seq_header] = []
                continue #go to the next line, which should be the fasta seq
            seq = line #store the entire line as a string variable
            fasta[seq_header].append(seq) #associate the header (key) with the seq (value)
    for key, value in fasta.items(): #take the key, value pair
        fasta[key]=''.join(value) #overwrite the existing list of broken up fasta seqs with a concatenated form, easier to write to file later
    return(fasta) #Return the fasta dict

#Make output folder, if it doesn't already exist
try:
    os.mkdir("outputs/")
except OSError:
    print("Creation of the directory %s failed, already exists")
else:
    print("Successfully created the output directory.")

#Set file paths using the first and second terms following the python script name
aligned_fasta=sys.argv[1] #Take the first term (should be your mafft fasta) that follows the python script name

aligned_dict=fasta_to_dict(aligned_fasta) #Use the function defined above to store the fasta as a dict
#print(aligned_dict['NP_047200'])

sites_in_raw_fasta=[170,540,763,1016,1152,1317,1652,1747,1774,1964]
list_of_cleavage_sites=['LQRQGNSV','LAPQHWKT','LTSQTLTE','VIKQGAAS','IRRQGLLT','LEPQGLKD','IRRQGNRV','QEPQAAYS','IQRQGISP','TTQQSLIV']
initial_gap_count=0
for x in range(0,len(aligned_dict['NP_047200'])):
    if aligned_dict['NP_047200'][x]=='-':
        initial_gap_count=initial_gap_count+1
    else:
        break
start_sites_in_aln=[x+initial_gap_count-1 for x in sites_in_raw_fasta]
#print(start_sites_in_aln)

start_iter=iter(start_sites_in_aln)
n_slider=next(start_iter)
add_gaps=0
site_map={}
while len(list_of_cleavage_sites) > 0:
    if list_of_cleavage_sites[0] in aligned_dict['NP_047200'][n_slider+add_gaps:n_slider+60+add_gaps].replace('-',''): 
        print("Found",list_of_cleavage_sites[0],'at',n_slider+add_gaps,':',n_slider+60+add_gaps)
        site_slider=-1
        for rev_match in range(60,-1,-1):
            if aligned_dict['NP_047200'][n_slider+add_gaps:n_slider+60+add_gaps][rev_match-1]==list_of_cleavage_sites[0][site_slider]:
                #print(list_of_cleavage_sites[0][site_slider],'at',n_slider+add_gaps+rev_match-1)
                if list_of_cleavage_sites[0] in site_map:
                    site_map[list_of_cleavage_sites[0]].insert(0,n_slider+add_gaps+rev_match-1)
                else:
                    site_map[list_of_cleavage_sites[0]] =[n_slider+add_gaps+rev_match-1]
                if len(list_of_cleavage_sites[0])+site_slider>0:
                    site_slider=site_slider-1
                else:
                    print(aligned_dict['NP_047200'][n_slider+add_gaps:n_slider+60+add_gaps],'\n')
                    break
        list_of_cleavage_sites.pop(0)
        if len(list_of_cleavage_sites) > 1:
            n_slider=next(start_iter)
        add_gaps=0
    else:
        add_gaps=add_gaps+1
print("Positions of cleavage sites in the reference sequence")
print(site_map,'\n')

print('Printing aligned 8mers for:','APY16382')
for sets in site_map.values():
    print(''.join([aligned_dict['APY16382'][i] for i in sets])) 
    #If - occurs to the left of P1, regex the original fasta and pull 1 more position to the left
    #If - occurs to the right of P1', regex the original fasta and pull 1 more position