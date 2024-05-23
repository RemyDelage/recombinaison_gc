import re
import argparse

def extract_genes_ids_from_gpff(gpff_file: str, tsv_file: str) -> dict:
    '''
    Extracts the gene ID or locus tag from the GPFF file 
    (GenBank Protein Feature File) and associates it with the version 
    (header ID of the FASTA protein file). Selected identifiers are extracted 
    from a TSV file containing primary transcripts information.

    Parameters
    ----------
        gpff_file : str
            The file which contains all the proteins information (sequences
            and annotations).
        tsv_file : str
            The file which contains the information about the primary 
            transcripts.

    Returns
    -------
        dict :
            A dictionary with the gene ID or locus tag as keys and their 
            corresponding version values.
    '''
    gene_versions = {}

    # Extract the IDs from the TSV file
    with open(tsv_file, 'r') as tsv:
        next(tsv)  # Skip header
        for line in tsv:
            columns = line.strip().split('\t')
            attributes = columns[2].split(';')
            for attribute in attributes:
                if attribute.startswith('Name='):
                    gene_name = attribute.split("=")[1]
                    gene_versions[gene_name] = None  

    # Access the "version" of the proteins                
    with open(gpff_file, 'r') as gpff:
        current_version = None
        for line in gpff:
            if line.startswith('VERSION'):
                version_line = line.strip()
                current_version = re.search(r'VERSION\s+(\S+)', version_line).group(1)
            if "locus_tag=" in line or "gene=" in line:
                match = re.search(r'(?:locus_tag|gene)="([^"]+)"', line)
                if match:
                    gene_id = match.group(1)
                    if gene_id in gene_versions:
                        gene_versions[gene_id] = current_version

    # Remove None values from the dictionary
    gene_versions = {k: v for k, v in gene_versions.items() if v is not None}                   
    return gene_versions

def extract_sequences_from_fasta(input_fasta: str, output_fasta: str, ids_to_extract: dict):
    '''
    This function allows to extract the sequences form the original FASTA file
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
    with open(input_fasta, 'r') as fasta_in, open(output_fasta, 'w') as fasta_out:
        write_sequence = False
        current_sequence = []
        # Extract the header and the sequences from the original FASTA file
        for line in fasta_in:
            if line.startswith('>'):
                identifier = line.split()[0][1:]
                 # Checks if the identifier is on the values of the dictionary
                if identifier in ids_to_extract:
                    # If it's a new sequence, write the previous sequence if exists
                    if current_sequence:
                        fasta_out.write(''.join(current_sequence) + '\n')
                        current_sequence = []
                    # Write the header
                    fasta_out.write(line)
                    write_sequence = True
                else:
                    write_sequence = False
            elif write_sequence:
                current_sequence.append(line.strip())

        if current_sequence:
            fasta_out.write(''.join(current_sequence) + '\n')

    # Inform the user that the output file has been created
    print("File created: ", output_fasta)


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
