#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 15:08:00 2025
Takes a fasta file and taxonomy table, outputs fasta file with taxonomy added to headers in SINTAX format.
The fasta file must be formatted with the BIN/OTU name (i.e. whatever name is used to look up the taxonomy
table) as the right half of the header, separated from other text by a "|" character.
The same names must be in the first column of the taxonomy table.
Other column names in the taxonomy table must be the taxonomic ranks. Any rank names which are not supported
by SINTAX will be excluded from the output.
@author: robinfloyd
"""

import sys
from Bio import SeqIO

def import_taxonomy(taxonomyfile):
    taxonomy_dict = {}    # create empty dictionary for results
    tax_rank_initials = []  # create empty list for taxonomic rank initials
    sintax_valid_ranks = ['kingdom','domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']

    with open(taxonomyfile, 'r') as tax_file:           
        for line in tax_file:            
            columns = line.strip().split('\t')  # Assuming tab-separated columns  
            
            if not(tax_rank_initials): # If the list of initials is not yet defined, we are on the first (header) row                     

                start_col_for_tax = min(i for i,v in enumerate(columns) if v.lower() in sintax_valid_ranks)
                # Find the first column name matching one of the recognized taxonomic ranks

                tax_rank_names = columns[start_col_for_tax:]
                print('Taxonomy table includes the following ranks:',', '.join(tax_rank_names))

                columns_to_exclude = list(a for a,b in enumerate(tax_rank_names) if b not in sintax_valid_ranks)
                # Any columns with non-valid ranks are excluded

                columns_to_exclude.sort(reverse=True)

                for x in columns_to_exclude:
                    tax_rank_names.pop(x)
                            
                print('Output will include the following SINTAX-supported ranks:',', '.join(tax_rank_names))
    
                for rank in tax_rank_names:    # Get the initial letters of the taxonomic ranks from the header row
                   tax_rank_initials.append(rank[0].lower())
                              
            else:                    
                bin_name = columns[0]  # First column is the BIN or sequence name          
                taxonomic_info = columns[start_col_for_tax:]    # Take all columns after the 1st, these contain the taxonomy

                for y in sorted(columns_to_exclude, reverse=True):
                    if y < len(taxonomic_info):
                        taxonomic_info.pop(y)

                # Missing values are denoted by 'NA'. If only the right-most values, fine; they will simply be
                # left out of the SINTAX output. But if any values are missing in the middle of the taxonomic
                # chain, it will create errors, so any such rows need to be excluded from the reference set and
                # output to a separate "problem" file for checking (likewise if a row lacks any taxonomic info).
    
                # The next loop checks whether the last item in the list is 'NA' and deletes it if it is. 
                # Keep doing this until the last item is not 'NA' or the entire list is deleted.
                while (len(taxonomic_info) > 0) and (taxonomic_info[-1] == 'NA'):
                    taxonomic_info = taxonomic_info[:-1]
    
                # Next, check either if we have deleted the entire list (row had no taxonomy at all), or if any 
                # values remain as 'NA' (not as the last item). In either case, this row needs to be output to
                # the "problems" file and not added to the reference set.
                if (len(taxonomic_info) == 0) or ('NA' in taxonomic_info):
                    print(line.strip(), file=open('problemtaxa.tsv','a'))    
                    continue # Exit the 'for' loop and continue to the next line in the BIN taxonomy table
                                  
                # If the current line made it through the previous checks, start a new string for the taxonomy.
                new_taxonomy_string = ''
    
                # Loop through each item in the taxonomy list, adding the appropriate initial for the rank.
                # Concatenate to the same taxonomy string each time through the loop.

                for init, name in zip(tax_rank_initials,taxonomic_info):
                    new_taxonomy_string = new_taxonomy_string + init + ':' + name + ','
                
                # Once all the names have been added, final comma needs to be deleted
                new_taxonomy_string = new_taxonomy_string[:-1]
                
                # Store the results in a dictionary by BIN name
                taxonomy_dict[bin_name] = new_taxonomy_string
        
    return taxonomy_dict

# Take the names of the fasta and taxonomy files as the first two arguments when called
fastafile = sys.argv[1]
taxonomyfile = sys.argv[2]

# Count sequences in the input fasta file
num = len([1 for line in open(fastafile) if line.startswith(">")])
print(num, 'sequences found in input fasta file.')

# Create the dictionary of SINTAX-formatted taxonomy for each BIN by calling the 'import_taxonomy' function above
taxonomy_dict = import_taxonomy(taxonomyfile)

# Create the default name for the output file by taking all text before the last '.' in the input
# filename, then adding "_sintax.fasta".
dataset_name = fastafile.rsplit('.',1)[0]
sintaxfile = dataset_name + '_sintax.fasta'

# Counter variable to track the number of items processed
count = 0

# Loop through each sequence in the input fasta file
for seq_record in SeqIO.parse(fastafile, "fasta"):
    # Capture the fasta header and BIN name
    oldheader = seq_record.id
    bin_name = oldheader.rsplit('|',1)[-1]

    # Check that the BIN name is in the dictionary of SINTAX-formatted taxonomy created earlier.
    # If no match is found using the second half of the header, try the first half.
    # Non-BIN sequences such as bacteria use the ProcessID as the identifier, which is the 1st part of the line.
    if not (bin_name in taxonomy_dict.keys()):
        print(bin_name, file=open('problemseqs.tsv','a'))
        continue # Exit the 'for' loop and continue to the next sequence in the input file

    # For sequences that are in the dictionary, get the taxonomy string
    new_taxonomy_string = taxonomy_dict[bin_name]

    # Create the new header by adding the SINTAX-formatted taxonomy string separated by ";"
    newheader = oldheader + ';tax=' + new_taxonomy_string + ';'

    # Write results to the output file.    
    print('>' + newheader +'\n' + seq_record.seq, file=open(sintaxfile,'a'))
    count += 1

print(count, 'sequences written to SINTAX output file.')
