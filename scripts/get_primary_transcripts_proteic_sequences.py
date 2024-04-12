import re

def get_sequences(gpff_file: str, genes_file: str):
    sequences = dict()
    current_locus_id = None
    current_locus_seq = []
    with open(gpff_file, 'r') as gpff:
        for line in gpff:  # Start reading from FEATURES section
            if line.startswith("VERSION"):
                # Extract locus ID from VERSION line
                locus_id_match = re.search(r'(\w+\.\d+)', line)
                if locus_id_match:
                    current_locus_id = locus_id_match.group(1)
            elif "gene=" in line or "locus_tag=" in line:
                cds_match = re.search(r'/locus_tag="([^"]+)"', line)
                if cds_match:
                    current_locus_tag = cds_match.group(1)
            elif line.startswith("ORIGIN"):
                for seq_line in gpff:
                    if seq_line.startswith('//'):
                        break
                    # Extract DNA sequence and remove spaces
                    sequence_part = re.sub(r'\d+|\s', '', seq_line.strip()).upper()  # Remove line numbers and spaces, convert to uppercase
                    current_locus_seq.append(sequence_part)
                # Add the sequence associated with the current locus ID
                if current_locus_id:
                    sequences[current_locus_id] = ''.join(current_locus_seq)
                    current_locus_seq = []
    return sequences

# Exemple d'utilisation
sequences = get_sequences("GCF_000001735.4_TAIR10.1_protein.gpff", "GCF_000001735.4_Arabidopsis_thaliana_primary_transcripts.tsv")
print(sequences)
