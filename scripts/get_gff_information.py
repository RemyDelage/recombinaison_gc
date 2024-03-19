# Import the modules to execute the script in command line
import argparse
import re

def extract_protein_id(line):
    '''
    This function allows to obtain the protein ID from the ggf files.

    Parameters
    ----------
        line : int
            The reading line of the gff file

    Returns
    -------
        - str : The protein ID if it is found
        - None otherwise 
    '''
    # Search in each line if the protein ID is mentionned
    match = re.search(r'protein_id=([^;\n]+)', line)
    if match:
        protein_id = match.group(1)
        #   Remove 'gnl|JCVI|' prefix if present
        protein_id = protein_id.replace('gnl|JCVI|', '')
        return protein_id
    else:
        return None

def extract_gene_id(line):
    '''
    This function allows to obtain the gene ID from the ggf files.

    Parameters
    ----------
        line : int
            The reading line of the gff file

    Returns
    -------
        - str : The gene ID if it is found
        - None otherwise 
    '''
    # Search in each line if the gene ID is present
    match_gene = re.search(r'gene=([^;\n]+)', line)
    if match_gene:
        return match_gene.group(1)
    else:
        # If the gene ID is not founded, checks the locus tag
        match_locus = re.search(r'locus_tag=([^;\n]+)', line)
        if match_locus:
            return match_locus.group(1)
        # If both are missing, printing a dot for missing data
        else:
            return "."


def gff_extraction(gff_file:str, tsv_file:str):
    '''
    This function allows to create a tsv file which contains all the 
    information obtained from the gff file.

    Parameters
    ----------
        gff_file : str
            The gff file which the informations will be extracted.

        tsv_file : str
            The tsv file will be created with the extracted informations

    Returns
    -------
        The tsv file
    '''
    # Open the gff file to extract data and tsv file to write the results
    with open(gff_file, 'r') as gff, open(tsv_file, "w") as tsv:
        # Create the header of the tsv file
        tsv.write("FEATURE\tCHROMOSOME_ID\tPROTEIN_ID\tGENE_ID\tSTRAND\tSTART\tEND\n")
        # Browse gff file
        for line in gff:
            # Ignore the first lines starts with "#"
            if not line.startswith("#"):
                fields = line.strip().split('\t')
                # Getting the feature type
                feature_type = fields[2]
                # Extract the gff information
                if feature_type in ['gene', 'exon', 'mRNA', 'CDS']:
                    protein_id = extract_protein_id(fields[8])
                    gene_id = extract_gene_id(fields[8])
                    if protein_id:
                        feature = fields[2]
                        chromosome = fields[0]
                        strand = fields[6]
                        start = fields[4]
                        end = fields[3]

                        # Write the tsv file
                        tsv.write(f"{feature}\t{chromosome}\t{protein_id}\t{gene_id}\t{strand}\t{start}\t{end}\n")

# Create a parser to run the script directly in command line
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=
                                    'Extract informations from GFF file')

    parser.add_argument('-gff', type = str, 
                    help = 'GFF file to analyse')
    
    parser.add_argument('-tsv', type = str, 
                    help = 'The tsv to store the data')
    
    args = parser.parse_args()

    gff_extraction(args.gff, args.tsv)
