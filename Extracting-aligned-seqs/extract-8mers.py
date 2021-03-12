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

#Temp, manually set location of cleavage sites (Reported by UniProt aka predicted + Literature)
#sites_in_raw_fasta=[170,540,763,1016,1152,1317,1652,1747,1774,1964]
#list_of_cleavage_sites=['LQRQGNSV','LAPQHWKT','LTSQTLTE','VIKQGAAS','IRRQGLLT','LEPQGLKD','IRRQGNRV','QEPQAAYS','IQRQGISP','TTQQSLIV']
#Our reference is NP_047200
sites_in_raw_fasta=[170,540,763,1152,1317,1652,1747,1774,1964] #Apparently VIKQGAAS is cut by 3CD
list_of_cleavage_sites=['LQRQGNSV','LAPQHWKT','LTSQTLTE','IRRQGLLT','LEPQGLKD','IRRQGNRV','QEPQAAYS','IQRQGISP','TTQQSLIV'] #Apparently VIKQGAAS is cut by 3CD

#Add an initial gap count to the raw fasta positions of the sites to pre-emptively speed up sliding window match
initial_gap_count=0
for x in range(0,len(aligned_dict['NP_047200'])):
    if aligned_dict['NP_047200'][x]=='-':
        initial_gap_count=initial_gap_count+1
    else:
        break
start_sites_in_aln=[x+initial_gap_count-1 for x in sites_in_raw_fasta]
#print(start_sites_in_aln)

start_iter=iter(start_sites_in_aln) #make the list iterable with next()
n_slider=next(start_iter) #use next() to store the first element
add_gaps=0 #value will increase with each encountered gap (-), to accurately store site location in the alignment
site_map={}
while len(list_of_cleavage_sites) > 0: #Run this loop as long as the cleavage list size is greater than zero
    #As cleavage sites are found within a sliding window of 60 amino acids...store site info and delete sites from cleavage list
    if list_of_cleavage_sites[0] in aligned_dict['NP_047200'][n_slider+add_gaps:n_slider+60+add_gaps].replace('-',''): 
        print("Found",list_of_cleavage_sites[0],'at',n_slider+add_gaps,':',n_slider+60+add_gaps)
        site_slider=-1 #Define slider that will act as a quering index starting at the end of the seq, to move in reverse
        for rev_match in range(60,-1,-1): #for-loop going in reverse order
            #both the sliding 60AA window and the cleavage site in question will be traced for matched sequences starting from the end of each seq
            if aligned_dict['NP_047200'][n_slider+add_gaps:n_slider+60+add_gaps][rev_match-1]==list_of_cleavage_sites[0][site_slider]:
                #print(list_of_cleavage_sites[0][site_slider],'at',n_slider+add_gaps+rev_match-1)
                if list_of_cleavage_sites[0] in site_map: #if this ref cleavage site already has a dict entry
                    site_map[list_of_cleavage_sites[0]].insert(0,n_slider+add_gaps+rev_match-1) #append the matched AA index to the front of the values list
                else: #if the ref cleavage site does NOT have a dict entry yet
                    site_map[list_of_cleavage_sites[0]] =[n_slider+add_gaps+rev_match-1] #make one and store the first value
                if len(list_of_cleavage_sites[0])+site_slider>0: #If length of cleavage site and slider are greater than zero (end of seq not yet reached)
                    site_slider=site_slider-1 #Make slider more negative, to continue tracing across the cleavage site in reverse order
                else:
                    print(aligned_dict['NP_047200'][n_slider+add_gaps:n_slider+60+add_gaps],'\n')
                    break #break for-loop and return to the while loop once all indices have been recorded
        list_of_cleavage_sites.pop(0) #This is where we delete the currently queried cleavage site from the cleavage list
        if len(list_of_cleavage_sites) > 1: #As long as there's still more than 1 seq remaining for the while loop
            n_slider=next(start_iter) #use next() to move to the next element in the iterable cleavage start position list
        add_gaps=0 #reset gap counter
    else:
        add_gaps=add_gaps+1 #keep adding to gap counter, until sliding window captures all of the cleavage site positions
print("Positions of cleavage sites in the reference sequence")
print(site_map,'\n')

#Testing for one case, but there's still a gap captured, meaning the alignment isn't ideal
#this means we need to fill in the gaps before we proceed with identity matching...
#print('Printing aligned 8mers for:','APY16382')
#for sets in site_map.values():
#    print(''.join([aligned_dict['APY16382'][i] for i in sets])) 
    #If - occurs to the left of P1, regex the original fasta and pull 1 more position to the left
    #If - occurs to the right of P1', regex the original fasta and pull 1 more position

#Thank goodness we defined this function to generate fasta dicts
raw_fasta=sys.argv[2]
raw_fasta_dict=fasta_to_dict(raw_fasta)

#Here we make a dictionary where each fasta ID (ID_key) is associated with its set of aligned clv seqs (the values)
#The strategy here is to fill in the gaps with neighboring residues to complete the 8mers
#Typically we want to avoid nested for-loops, but this one is pretty fast so we can look the other way here
#HUGE ASSUMPTION: P1 and P1' are aligned across all sequences, if they are meaningful alignments
separated_cleavage_dict={}
for ID_key in raw_fasta_dict: #iterate through IDs in fasta dictionary, doesn't matter which (aligned or raw)
    for cleavage_positions in site_map.values(): #iterate through cleavage positions we found from the ref sequence
        aligned_string=''.join([aligned_dict[ID_key][i] for i in cleavage_positions]) #combine the list of AAs into a single string of AAs
        gaps_right_of_p1prime=aligned_string[4:].count('-') #count how many gaps exist to the right of the P1-P1' clv site
        gaps_left_of_p1=aligned_string[:3].count('-') #count how many gaps exist to the left of the P1-P1' clv site
        ungapped_aligned_string=aligned_string.replace('-','') #Remove gaps to prepare string search
        raw_string_loc=raw_fasta_dict[ID_key].find(ungapped_aligned_string) #String search to find where the match exists
        full_site=raw_fasta_dict[ID_key][raw_string_loc-gaps_left_of_p1:raw_string_loc+len(ungapped_aligned_string)+gaps_right_of_p1prime] #Recover the full string, accounting for loss in positions due to gaps
        if ID_key in separated_cleavage_dict: #if this ID_key already has a dict entry
            separated_cleavage_dict[ID_key].append(full_site) #append the matched clv site to the front of the values list
        else: #if the ID_key does NOT have a dict entry yet
            separated_cleavage_dict[ID_key] =[full_site] #make one and store the first value aka cleavage site

#Begin cleanup
print(len(separated_cleavage_dict),"IDs before cleanup")

#The idea is to select polyprotein cleavage sets that are similarly represented across all sequences
#This means if even one site is missing, or uncertain (X), the entire set should be removed
#Likewise, we want to remove duplicates so that we don't mask signals in weaker motif positions (broader motif is acceptable for us)
separated_cleavage_dict_noblanks={k: v for k, v in separated_cleavage_dict.items() if '' not in v} #Dict comprehension, remove entries with missing clv sites
print(len(separated_cleavage_dict_noblanks),"IDs after removing missing cleavage sites")
separated_cleavage_dict_noblanks_x={k: v for k, v in separated_cleavage_dict_noblanks.items() if 'X' not in ''.join(v)} #Dict comprehension, remove entries with uncertain residues
print(len(separated_cleavage_dict_noblanks_x),"IDs after removing cleavage sites with uncertainty (X residues)")
# Remove duplicate values in dictionary 
temp = [] #This list compiles the first occurrence of a unique value
separated_cleavage_dict_noblanks_x_nodupes = {} #This is where we'll store our not dupes
for ID_key, cleavage_set in separated_cleavage_dict_noblanks_x.items(): #iterate through each dictionary pair
    if cleavage_set not in temp: #if this is the first instance of the value we've come across
        temp.append(cleavage_set) #put it in the list
        separated_cleavage_dict_noblanks_x_nodupes[ID_key] = cleavage_set #and store the ID_key, cleavage set pair
print(len(separated_cleavage_dict_noblanks_x_nodupes),"IDs after removing duplicate cleavage site sets")
print(list(separated_cleavage_dict_noblanks_x_nodupes.items())[-1]) #Just print the last entry of the dictionary, to make sure the filters worked


####Extra test garbage
#concat all the values in the site_map dict and search rapidly through each seq, then separate by 8mers
#concat_site_map_vals=sum(site_map.values(), [])
#aligned_dict_extraction={}
#for ID_keys in aligned_dict:
#    aligned_dict_extraction[ID_keys]=''.join([aligned_dict[ID_keys][i] for i in concat_site_map_vals])

#count_entries_with_gaps=0
#for vals in aligned_dict_extraction.values():
    #print(vals)
#    if '-' in vals:
#        count_entries_with_gaps=count_entries_with_gaps+1

#print(count_entries_with_gaps)
#print(len(aligned_dict_extraction))


