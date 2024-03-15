# Import the modules to execute the script in command line
import argparse
import re

def extract_protein_id(line):
    match = re.search(r'protein_id=([^;\n]+)', line)
    if match:
        protein_id = match.group(1)
        #   Remove 'gnl|JCVI|' prefix if present
        protein_id = protein_id.replace('gnl|JCVI|', '')
        return protein_id
    else:
        return None

def extract_gene_id(line):
    match = re.search(r'gene=([^;\n]+)', line)
    if match:
        return match.group(1)
    else:
        return "."

def gff_extraction(gff_file:str, tsv_file:str):
    # Open the gff file to extract data and tsv file to write the results
    with open(gff_file, 'r') as gff, open(tsv_file, "w") as tsv:
        tsv.write("FEATURE\tCHROMOSOME_ID\tPROTEIN_ID\tGENE_ID\tSTRAND\tSTART\tEND\n")
        for line in gff:
            if not line.startswith("#"):
                fields = line.strip().split('\t')
                feature_type = fields[2]

                if feature_type in ['gene', 'exon', 'mRNA', 'CDS']:
                    protein_id = extract_protein_id(fields[8])
                    gene_id = extract_gene_id(fields[8])
                    if protein_id:
                        feature = fields[2]
                        chromosome = fields[0]
                        strand = fields[6]
                        start = fields[4]
                        end = fields[3]

                        tsv.write(f"{feature}\t{chromosome}\t{protein_id}\t{gene_id}\t{strand}\t{start}\t{end}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=
                                    'Extract informations from GFF file')

    parser.add_argument('-gff', type = str, 
                    help = 'GFF file to analyse')
    
    parser.add_argument('-tsv', type = str, 
                    help = 'The tsv to store the data')
    
    args = parser.parse_args()

    gff_extraction(args.gff, args.tsv)
