# Import necessary packages
import glob
from os.path import join
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import shutil


# Set working directory and needed files
data_dir = './'
gff_info_files = glob.glob(data_dir + '*_gff_informations.tsv')
genes_info_files = glob.glob(data_dir + '*_genes_informations.tsv')
genomes_fasta_files = glob.glob(data_dir + '*_genomic.fna')


#print(accession_number[species])

def get_protein_id():
    '''
    This funtion allows to get the ID of each protein of the orthogroup according 
    to the species.

    Returns
    -------
    lists :
        The proteins ID (one list for each species)
    '''

    # Lists to store proteins id
    arabidopsis_thaliana_prot_id = []
    arabidopsis_lyrata_prot_id = []
    brassica_rapa_prot_id = []

    # Read Orthogroups.tsv file (choose a specific line)
    with open("Orthogroups.tsv", 'r') as orthogroups:
        for i, line in enumerate(orthogroups):
            if i == 26: # One specific orthogroup (--> will be changed for taking account each orthogroup)
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
    #print("Tests")
    #print("-----------------------------------------------", '\n')
    #print("Protein_IDs", '\n')
    #print("Arabidopsis thaliana Protein IDs:", arabidopsis_thaliana_prot_id)
    #print("Arabidopsis lyrata Protein IDs:", arabidopsis_lyrata_prot_id)
    #print("Brassica rapa Protein IDs:", brassica_rapa_prot_id)

    return arabidopsis_thaliana_prot_id, arabidopsis_lyrata_prot_id, brassica_rapa_prot_id

# Execute the function and store the results in global variables
arabidopsis_thaliana_prot_id, arabidopsis_lyrata_prot_id, brassica_rapa_prot_id = get_protein_id()

def get_gene_id():
    '''
    This function allows to get the genes ID associated with the proetins ID 
    previously founded.

    Returns
    -------
    dict
        One dictionary per species with in key the protein ID and in value the
        associated gene ID
    '''
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
    #print("-----------------------------------------------")
    #print("\nGene IDs")
    #print("Arabidopsis thaliana Gene IDs:", arabidopsis_thaliana_genes_id)
    #print("Arabidopsis lyrata Gene IDs:", arabidopsis_lyrata_genes_id)
    #print("Brassica rapa Gene IDs:", brassica_rapa_genes_id)

    return arabidopsis_thaliana_genes_id, arabidopsis_lyrata_genes_id, brassica_rapa_genes_id

# Execute the function
arabidopsis_thaliana_genes_id, arabidopsis_lyrata_genes_id, brassica_rapa_genes_id = get_gene_id()

def get_chromosomes_id():
    '''
    This function allows to obtain the chromosomes ID for obtain the genes
    locarion into the fasta files

    Returns
    -------
    dict
        One dictionary per species, with in key the gene ID and in value the
        chromosome ID where the gene is located
    '''
    # Dictionaries to store the genes IDs associated with the proteins IDs
    arabidopsis_thaliana_chrom_id = {} 
    arabidopsis_lyrata_chrom_id = {}
    brassica_rapa_chrom_id = {}

    # Read all the gff_informations.tsv files
    for file in gff_info_files:
        with open (file, 'r') as gff_info:
            for line in gff_info:
                fields = line.strip().split('\t')
                # Get the chromosomes IDs
                chromosome_id = fields[1]
                protein_id = fields[2]

                if  protein_id in arabidopsis_thaliana_prot_id:
                    arabidopsis_thaliana_chrom_id[protein_id] = chromosome_id
                elif protein_id in arabidopsis_lyrata_prot_id:
                    arabidopsis_lyrata_chrom_id[protein_id] = chromosome_id
                elif protein_id in brassica_rapa_prot_id:
                    brassica_rapa_chrom_id[protein_id] = chromosome_id    
           
    # Testing the function
    #print("-----------------------------------------------")
    #print("\nChromosomes IDs")
    #print("Arabidopsis thaliana Chromosome IDs:", arabidopsis_thaliana_chrom_id)
    #print("Arabidopsis lyrata Chromosome IDs:", arabidopsis_lyrata_chrom_id)
    #print("Brassica rapa Chromosome IDs:", brassica_rapa_chrom_id)

    return arabidopsis_thaliana_chrom_id, arabidopsis_lyrata_chrom_id, brassica_rapa_chrom_id

# Execute the function
get_chromosomes_id()


def get_nucleic_sequences(genes_info_files, genomes_fasta_files):
    '''
    This function allows to obtain the nucleic sequence from the species 
    genomes files (FASTA format) according to the information contains into the
    species gff files.

    Parameters
    ----------
        genes_info_files : str
            The tsv files which contains the information extracted from the gff
            files (one file per species).
        
        genomes_fasta_files : str
            The whole genome FASTA files of the specices which the 
            nucleic sequences will be extracted.

    Returns
    -------
        None.
    '''
    # Obtain the accession number and the species names from the input file names
    for genes_file in genes_info_files:
        genes_file = genes_file.replace('./', '')
        parts = genes_file.split('_')
        accession = join(parts[0], parts[1]).replace('/', '_')
        species = join(parts[2], parts[3]).replace('/', '_')
        
        # Create dictionary to store the results
        gene_id_dict = {}

        # Output directories
        species_output_dir = f'./{species.lower()}_output'
        if not os.path.exists(species_output_dir):
            os.makedirs(species_output_dir)  

        # Open the genes_information.tsv files and extract the informations
        with open(join('./', genes_file), 'r') as genes_info:
            next(genes_info) # Don't read the header
            for line in genes_info:
                info = line.strip().split('\t')
                gene_id = info[1]
                start = int(info[3])
                end = int(info[4])
                strand = info[2]
                chromosome = info[5]

                # Reads each FASTA genomes files
                for genome_file in genomes_fasta_files:
                    with open(genome_file, 'r') as genome_fasta:
                        for record in SeqIO.parse(genome_fasta, "fasta"):
                            if record.id == chromosome:
                                # Taking account of the strand of the gene
                                if strand == "+":
                                    sequence = record.seq[start-1:end]
                                else:
                                    sequence = record.seq[start-1:end].reverse_complement()

                                # Create the output files names
                                output_filename = f"{accession}_{species}_{gene_id}.faa"
                                # Store the results into the dictionary
                                gene_id_dict[gene_id] = output_filename
                                # Determine the output directory
                                output_path = os.path.join(species_output_dir, output_filename)

                                # Write the output file
                                with open(output_path, 'w') as output_file:
                                    SeqIO.write(SeqRecord(sequence, id=gene_id, description=''), output_file, "fasta")
                                    print(f"Gene sequence {gene_id} of {species} written in {output_path}")
                                    

# Function execution
get_nucleic_sequences(genes_info_files, genomes_fasta_files)
