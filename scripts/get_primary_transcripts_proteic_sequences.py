import re

def extract_genes_ids_from_gpff(gpff_file, tsv_file):
    versions_names = {}

    # Extract the genes IDs or proteins IDs from the tsv file
    with open(tsv_file, 'r') as tsv:
        next(tsv)
        for line in tsv:
            columns = line.strip().split('\t')
            attributes = columns[2].split(';')
            for attribute in attributes:
                if attribute.startswith('Name='):
                    gene_name = attribute.split("=")[1]
                    versions_names[gene_name] = None  # Initialiser avec None

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
                    if gene_id in versions_names:
                        versions_names[gene_id] = current_version

    # Supprimer les valeurs None du dictionnaire résultant
    versions_names = {k: v 
    for k, v in versions_names.items() 
    if v is not None}

    #print(versions_names)                        
    return versions_names

# Test de la fonction
versions_names = extract_genes_ids_from_gpff("GCF_000001735.4_TAIR10.1_protein.gpff", 
        "GCF_000001735.4_Arabidopsis_thaliana_primary_transcripts_updated.tsv")
#print(versions_names.values())

def extract_sequences_from_fasta(input_fasta, output_fasta, ids_to_extract):
    with open(input_fasta, 'r') as fasta_in, open(output_fasta, 'w') as fasta_out:
        write_sequence = False
        current_sequence = []

        for line in fasta_in:
            if line.startswith('>'):
                identifier = line.split()[0][1:]
                if identifier in ids_to_extract:
                    if line.startswith(">"):
                        fasta_out.write(line)
                    write_sequence = True
                    current_sequence.append(line)
                else:
                    write_sequence = False
            elif write_sequence:
                # Remove newline characters and append the sequence to the current sequence list
                current_sequence.append(line.strip())

        # Write the sequence to output file as a single line
        fasta_out.write(''.join(current_sequence) + '\n')

    print("File created : ", output_fasta)


extract_sequences_from_fasta("GCF_000001735.4_TAIR10.1_protein.faa",
                            "test3.faa", set(versions_names.values()))
