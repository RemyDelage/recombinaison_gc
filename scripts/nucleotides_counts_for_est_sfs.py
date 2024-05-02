# -*- coding:utf-8 -*-

import gzip
from collections import Counter
import glob
import os
import argparse


def nucleotides_counts(bed_file: str):
    '''
    This function allows to count the occurences of nucleotides for each sites 
    and for each species (the reference species and two outgroups) for a 
    BED file (compressed).

    Parameters
    ----------
        bed_file:str
            The compressed BED file which contain the sequences.
    
    Returns
    -------
        bases_counts : list
            A list of nucleotide counts for each site and species 
    '''
    bases_counts=[] # Initialize the output list
    with gzip.open(bed_file, "rt") as bed: # Unzip and read the BED file
        for line in bed:
            # Access to the sequences
            fields = line.strip().split("\t")
            sequences = fields[7].split(",")
            # Transform the transform lower-case characters in upper-case character to make counting easier
            sequences = [seq.upper() for seq in sequences]
            # Define each sequences for each species
            seq_ref_species = sequences[0]
            seq_outgroup1 = sequences[1]
            seq_outgroup2 = sequences[2]

            # Count the nucleotides occurences for each species (ignore hyphens and interrogation marks)
            count_ref_species = Counter(seq_ref_species)
            count_ref_species.pop("?", None)
            count_ref_species.pop("-", None)
            
            count_outgroup1 = Counter(seq_outgroup1)
            count_outgroup1.pop("?", None)
            count_outgroup1.pop("-", None)

            count_outgroup2 = Counter(seq_outgroup2)
            count_outgroup2.pop("?", None)
            count_outgroup2.pop("-", None)

            # Create a count list for each site
            bases_count_list=[]
            # Get all the nucleotides counts for each site (not only those present)
            bases_count_list.append({'A': count_ref_species['A'], 'T': count_ref_species['T'], 'G': count_ref_species['G'], 'C': count_ref_species['C']})
            bases_count_list.append({'A': count_outgroup1['A'], 'T': count_outgroup1['T'], 'G': count_outgroup1['G'], 'C': count_outgroup1['C']})
            bases_count_list.append({'A': count_outgroup2['A'], 'T': count_outgroup2['T'], 'G': count_outgroup2['G'], 'C': count_outgroup2['C']})
            
            # Add all the lists (one for each site) into the output list
            bases_counts.append(bases_count_list)
    return bases_counts
            

def format_counts(bases_counts: list):
    '''
    This function reformats the output list of the nucleotides_counts function 
    so that it can be used in EST-SFS software (Keightley et al., 2016).

    Parameters
    ----------
        bases_counts: list
            The list of counts produced by the nucleotides counts function.

    Returns
    -------
        formatted_counts: list
            The counts formated nucleotides counts for using EST-SFS :
            The reference species and the first outgroup will be separated by 
            tabulation and the outgroup species will be separated with a space. 
    '''
    # Initialize tge ouput list
    formatted_counts = []
    for site in bases_counts:
        formatted_site = []
        # Update the format for each site
        for i, species_count in enumerate(site):
            formatted_species_counts = ','.join(str(species_count[base]) for base in "ATGC")
            if i == 0:
                formatted_site.append(formatted_species_counts)
            elif i == 1:
                formatted_site.append('\t' + formatted_species_counts)
            else:
                formatted_site.append(' ' + formatted_species_counts)
        # Add the updated list into the output list
        formatted_counts.append(''.join(formatted_site))
    return formatted_counts


def write_output(bases_counts: list, output_file: str):
    '''
    This function allows to write the output file which will be used in EST-SFS
    software (Keightley et al., 2016).

    Parameters
    ----------
        bases_counts: list
            The list of nucleotides counts produces by the nucleotides_counts 
            function.

        output_file: str
            The name of the file which will be created and used with EST-SFS.

    Returns
    -------
        None
            A message informs the user that output files have been created.
    '''
    # Reformat the nucleotide counts thanks to the format_counts function
    formatted_counts = format_counts(bases_counts)
    # Write the output file
    with open(output_file, "w") as output:
        for line in formatted_counts:
            cols = line.strip().split("\t")
            # Check if the line contains at least 3 columns and each character is either '0', '1', ',', '\t', or ' '
            if len(cols) >= 3 and all(char in ['0', '1', ',', '\t', ' '] for char in line):
                output.write(line + "\n")
    # Print the output message
    print(f"{output_file} created")




def main():
    '''
    This function allows to create the parser for automatically run this script
    in a single command line
    '''
    parser = argparse.ArgumentParser(description = "Script for count the nucleotides occurences \
    for using EST-SFS Software (Keightley et al., 2016)")
    parser.add_argument("-i", "--input_dir", type = str, required = True, \
    help = "The directory containing the BED files (compressed)")
    parser.add_argument("-o", "--output_dir", type = str, required = True, \
    help = "The output file (TXT format) which contains the nucleotides counts (will use for EST-SFS)")
    args = parser.parse_args()

    # Find the BED files into the directory
    bed_files = glob.glob(os.path.join(args.input_dir, "*.bed.gz"))

    for bed_file in bed_files:
        output_file = os.path.join(args.output_dir, os.path.basename(bed_file).replace(".bed.gz", ".data_file.txt"))

        bases_counts = nucleotides_counts(bed_file)
        write_output(bases_counts, output_file) 

# Script execution
if __name__ == "__main__":
    main()