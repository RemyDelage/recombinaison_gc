import re

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

def extract_sequences_from_fasta(input_fasta, output_fasta, ids_to_extract):
    with open(input_fasta, 'r') as fasta_in, open(output_fasta, 'w') as fasta_out:
        write_sequence = False
        current_sequence = []

        for line in fasta_in:
            if line.startswith('>'):
                identifier = line.split()[0][1:]
                if identifier in ids_to_extract:
                    write_sequence = True
                    current_sequence.append(line)
                else:
                    write_sequence = False
                    if current_sequence:
                        fasta_out.write(''.join(current_sequence))
                        fasta_out.write('\n')  # Ajouter un saut de ligne après chaque séquence
                        current_sequence = []
            elif write_sequence:
                current_sequence.append(line.strip())

        if current_sequence:
            fasta_out.write(''.join(current_sequence))
            fasta_out.write('\n')  # Ajouter un saut de ligne après la dernière séquence

    print("File created: ", output_fasta)

# Example usage
versions_names = extract_genes_ids_from_gpff("GCF_000001735.4_TAIR10.1_protein.gpff", 
                                             "GCF_000001735.4_Arabidopsis_thaliana_primary_transcripts_updated.tsv")

extract_sequences_from_fasta("GCF_000001735.4_TAIR10.1_protein.faa",
                             "test3.faa", set(versions_names.values()))