import re
import argparse

def extract_genes_ids_from_gpff(gpff_file: str, tsv_file: str):
    '''
    This function extracts the gene ID or locus tag from the GPFF file 
    (GenBank Protein Feature File) and associates it with the version 
    (header ID of the FASTA protein file). Selected identifiers are extracted 
    from a TSV file containing primary transcripts information.

    Parameters
    ----------
        gpff_file : str
            The file which contains all the proteins information (sequences
            and annotations).
        tsv_file : str
            The file which contains the informations about the primary 
            transcripts.

    Returns
    -------
        dict :
            A dictionary with the "version" values in values and the
            gene ID or locus tag in keys.
    '''
    # Initialize the output dictionary
    versions_names = {}

    # Extract the IDs from the TSV file
    with open(tsv_file, 'r') as tsv:
        # Skip the header
        next(tsv) 
         # Browse the file
        for line in tsv:
            columns = line.strip().split('\t')
            attributes = columns[2].split(';')
            # Access to the attributes
            for attribute in attributes:
                # Get the ID
                if attribute.startswith('Name='):
                    gene_name = attribute.split("=")[1]
                    versions_names[gene_name] = None  
    
    # Access to the "version" of the proteins                
    with open(gpff_file, 'r') as gpff:
        current_version = None
        # Browse the file
        for line in gpff:
            # Get the version
            if line.startswith('VERSION'):
                version_line = line.strip()
                current_version = re.search(r'VERSION\s+(\S+)', version_line).group(1)
            # Access the gene ID or locus tag
            if "locus_tag=" in line or "gene=" in line:
                match = re.search(r'(?:locus_tag|gene)="([^"]+)"', line)
                # Add the IDs into the dictionary
                if match:
                    gene_id = match.group(1)
                    # Add the versions into the dictionary
                    if gene_id in versions_names:
                        versions_names[gene_id] = current_version

    # Remove None values of the dictionary
    versions_names = {k: v for k, v in versions_names.items() if v is not None}
                        
    return versions_names



def extract_sequences_from_fasta(input_fasta : str, output_fasta : str, ids_to_extract : str):
    '''
    This function allows to extract the sequences form the original FATSA file
    of the species to only keep the sequences of the primary transcripts.

    Parameters
    ----------
        input_fasta : str
            The original FASTA file of the species proteins (directly 
            downloaded from the NCBI database).
        output_fasta : str
            The created FASTA file with only the sequences of the primary 
            transcripts.
        ids_to_extract : dict
            The "versions" values of the primary transcripts previously 
            stored into the dictionary.

    Returns
    -------
        The new FASTA file with the sequences of the primary transcripts.
    '''
    # Open the two files
    with open(input_fasta, 'r') as fasta_in, \
     open(output_fasta, 'w') as output:
        write_sequence = False
        current_sequence = []
        # Extract the header and the sequences from the original FASTA file
        for line in fasta_in:
            if line.startswith('>'):
                identifier = line.split()[0][1:]
                # Checks if the identifier is on the values of the dictionary
                if identifier in ids_to_extract:
                    # If so, write the header of the new file
                    if line.startswith(">"):
                        fasta_out.write(line)
                    # Get the proteic sequence
                    write_sequence = True
                    current_sequence.append(line)
                else:
                    write_sequence = False
            elif write_sequence:
                # Remove newline characters and append the sequence to the current sequence list
                current_sequence.append(line.strip())

        # Write the sequence to output file as a single line
        output.write(''.join(current_sequence) + '\n')

    # Message for the user that the output file has been created.
    print("File created : ", output_fasta)




if __name__ == "__main__":
    # Create the parser for execute the script directly in command lines:
    parser = argparse.ArgumentParser\
    (description = 'Get a proteic FASTA file with only the sequences of the primary transcripts')
    parser.add_argument('-gpff', required = True, type = str,
    help = 'The GPFF file of the species (from GeneBank database)')
    parser.add_argument('-tsv', required = True, type = str,
    help = 'The TSV file which contains the information about primary transcipts')
    parser.add_argument('-fasta',  required = True, type = str,
    help = 'The original proteic FASTA file of the sepcies (from GeneBank database)')
    parser.add_argument('-o', required = True, type = str,
    help = 'The proteic FASTA file which only contains the protein sequences of the primary transcripts')

    args = parser.parse_args()

    # Execute the "extract_genes_ids_from_gpff" function :
    versions_names = extract_genes_ids_from_gpff(args.gpff, args.tsv) 
    
    # Execute the "extract_sequences_from_fasta" function :
    extract_sequences_from_fasta(args.fasta, args.o, set(versions_names.values()))