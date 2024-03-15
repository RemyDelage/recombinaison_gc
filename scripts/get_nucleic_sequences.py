# Import necessary packages
import glob

# Set working directory and needed files
data_dir = './'
gff_info_files = glob.glob(data_dir + '*_gff_informations.tsv')
genes_info_files = glob.glob(data_dir + '*_genes_informations.tsv')
genomes_fasta_files = glob.glob(data_dir + '*_genomic.fna')

def get_protein_id():
    # Lists to store proteins id
    arabidopsis_thaliana_prot_id = []
    arabidopsis_lyrata_prot_id = []
    brassica_rapa_prot_id = []

    # Read Orthogroups.tsv file (choose a specific line)
    with open("Orthogroups.tsv", 'r') as orthogroups:
        for i, line in enumerate(orthogroups):
            if i == 32:
                # Separate the fields of the line
                fields = line.strip().split('\t')
                
                # Column 2 treatment (Arabidopsis thaliana)
                arabidopsis_thaliana_ids = fields[1].split(', ')
                for protein_id in arabidopsis_thaliana_ids:
                    arabidopsis_thaliana_prot_id.append(protein_id)
                    
                # Column 3 treatment (Arabidopsis lyrata)
                arabidopsis_lyrata_ids = fields[2].split(', ')
                for protein_id in arabidopsis_lyrata_ids:
                    arabidopsis_lyrata_prot_id.append(protein_id)
                    
                # Column 4 treatment (Brassica rapa)
                brassica_rapa_ids = fields[3].split(', ')
                for protein_id in brassica_rapa_ids:
                    brassica_rapa_prot_id.append(protein_id)
                    
                break  # Exit the loop after line treatment

    # Testing the function
    print("Tests")
    print("-----------------------------------------------", '\n')
    print("Protein_IDs", '\n')
    print("Arabidopsis thaliana Protein IDs:", arabidopsis_thaliana_prot_id)
    print("Arabidopsis lyrata Protein IDs:", arabidopsis_lyrata_prot_id)
    print("Brassica rapa Protein IDs:", brassica_rapa_prot_id)

    return arabidopsis_thaliana_prot_id, arabidopsis_lyrata_prot_id, brassica_rapa_prot_id

# Execute the function and store the results in global variables
arabidopsis_thaliana_prot_id, arabidopsis_lyrata_prot_id, brassica_rapa_prot_id = get_protein_id()

def get_gene_id():
    # Dictionaries to store the genes IDs associated with the proteins IDs
    arabidopsis_thaliana_genes_id = {} 
    arabidopsis_lyrata_genes_id = {}
    brassica_rapa_genes_id = {}

    # Read the lines of gff info file
    for file in gff_info_files:
        with open(file, 'r') as gff_info:
            for line in gff_info:
                # Separate the fields of the line
                fields = line.strip().split('\t')
                # Get the genes and proteins IDs
                protein_id = fields[2]
                gene_id = fields[3] 
                
                # Checks if the protein_id is present on the gff_informations file
                if protein_id in arabidopsis_thaliana_prot_id:
                    arabidopsis_thaliana_genes_id[protein_id] = gene_id
                elif protein_id in arabidopsis_lyrata_prot_id:
                    arabidopsis_lyrata_genes_id[protein_id] = gene_id
                elif protein_id in brassica_rapa_prot_id:
                    brassica_rapa_genes_id[protein_id] = gene_id

    # Testing the function
    print("-----------------------------------------------")
    print("\nGene IDs")
    print("Arabidopsis thaliana Gene IDs:", arabidopsis_thaliana_genes_id)
    print("Arabidopsis lyrata Gene IDs:", arabidopsis_lyrata_genes_id)
    print("Brassica rapa Gene IDs:", brassica_rapa_genes_id)

get_gene_id()
