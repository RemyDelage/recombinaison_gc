from Bio import SeqIO
import re
import sys

def extract_genes_names(gff_file: str):
    genes_id = {}
    with open(gff_file, 'r') as gff:
        for line in gff:
            if not line.startswith("#"):  
                info = line.strip().split('\t')
                if len(info) > 8:  
                    attributes = info[8]
                    gene_id_match = re.search(r'Name=([^;]+)', attributes)
                    if gene_id_match:
                        gene_id = gene_id_match.group(1)
                    else:
                        dbxref_match = re.search(r'Dbxref=([^;]+)', attributes)
                        if dbxref_match:
                            gene_id = dbxref_match.group(1)
                        else:
                            continue
                        
                    chromosome_id = info[0]
                    start = int(info[3]) - 1
                    end = info[4]
                    strand = info[6]
                    gene_coordinates = f"{chromosome_id}:{start}-{end}({strand})"
                    genes_id[gene_coordinates] = gene_id
    print(genes_id)
    return genes_id


def clean_fasta_header(header):
    # Supprimer les caractères '>' et les espaces dans les en-têtes FASTA
    return header.lstrip(">").strip()

def update_fasta_headers(gff_file: str, fasta_file: str, output_file: str):
    genes_id = extract_genes_names(gff_file)
    with open(output_file, 'w') as output:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            header = clean_fasta_header(record.id)
            if header in genes_id:
                gene_name = genes_id[header]
                record.id = gene_name
                print("Gene found : ", gene_name)
                output.write(f">{gene_name}\n")
            else:
                print("Gene name not found. Writing original header : ", header)
                output.write(f">{header}\n")
            output.write(f"{record.seq}\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python update_fasta_headers.py <gff_file> <fasta_file> <output_file>")
        sys.exit(1)

    gff_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_file = sys.argv[3]
    update_fasta_headers(gff_file, fasta_file, output_file)
