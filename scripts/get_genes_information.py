# Import the modules 
import argparse
import re
from get_gff_information import extract_gene_id

def get_gene_information(gff_file:str, tsv_file:str):
    with open(gff_file, 'r') as gff, open(tsv_file, 'w') as tsv:
        tsv.write("FEATURE\tGENE_ID\tSTRAND\tSTART\tEND\tCHROMOSOME_ID\n")
        for line in gff:
            if not line.startswith("#"):
                fields = line.strip().split('\t')
                feature_type = fields[2]

                if feature_type == "gene":
                    gene_id = extract_gene_id(fields[8])
                    if gene_id:
                        feature = fields[2]
                        chromosome = fields[0]
                        strand = fields[6]
                        start = fields[3]
                        end = fields[4]
                        tsv.write(f"{feature}\t{gene_id}\t{strand}\t{start}\t{end}\t{chromosome}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=
                                    'Extract informations from GFF file')

    parser.add_argument('-gff', type = str, 
                    help = 'GFF file to analyse')
    
    parser.add_argument('-tsv', type = str, 
                    help = 'The tsv to store the data')
    
    args = parser.parse_args()

    get_gene_information(args.gff, args.tsv)
