import sys
import re
import argparse

def correspondance_transcript_genes(gff_file:str, output_file:str):
    """
    This function allows to establish the correspondance between the genes IDs 
    and the transcripts IDs

    Parameters
    ----------
        - gff_file : str
            The input GFF file to extract the genes ID and the transcript ID.
        - output_file : str
            The output file which the IDs were strored.
    
    Returns
    -------
        A TSV file which contains the genes IDs and the proteins IDs 
        (or transcripts IDs) associated.
    """
    # Store the unique informations of each gene into a dictionnary
    gene_info = {}  
    # Read the GFF input file
    with open(gff_file, 'r') as gff:
        # Browse the GFF file
        for line in gff:
            # Ignore the files with "#" (commentary lines)
            if not line.startswith("#"):
                # Divide each fields by a tabulation
                fields = line.strip().split('\t')
                if len(fields) >= 9:
                    # Get the attributes
                    attributes = fields[8]
                    # Get the genes IDs
                    gene_match = re.search(r'gene=([^;]+)', attributes)
                    # Store the genes IDs
                    gene_id = gene_match.group(1) if gene_match else None
                    # Get the proteins IDs
                    protein_match = re.search(r'protein_id=([^;]+)', attributes)
                    # Store the proteins IDs
                    protein_id = protein_match.group(1) if protein_match else None
                    # Get the Transcripts IDs
                    transcript_matches = re.findall(r'transcript_id=([^;]+)', attributes)
                    transcript_id = None # Initialisation of the variable
                    # Checks if the transcript_matches variable is not empty
                    if transcript_matches:
                        # Access to the last element "transcript_id = "
                        last_transcript = transcript_matches[-1]
                        # Checks if the the identifier of "transcript_id" not contains ""gnl|JCVI|""
                        if not last_transcript.startswith("gnl|JCVI|"):
                            # Store the transcript ID
                            transcript_id = last_transcript
                    
                    # Checks if a gene ID is found
                    if gene_id:
                        # Checks if the gene ID is in gene_info dictionnay
                        if gene_id not in gene_info:
                            # Add the gene ID into the gene_info dictionnary to store it
                            gene_info[gene_id] = {"transcript_id": transcript_id, "protein_id": protein_id}
                        else:
                            # If the gene ID is already in gene_info, associate the corresponding transcript ID
                            if transcript_id and not gene_info[gene_id]["transcript_id"]:
                                gene_info[gene_id]["transcript_id"] = transcript_id
                            if protein_id and not gene_info[gene_id]["protein_id"]:
                                gene_info[gene_id]["protein_id"] = protein_id

    # Create the output file
    with open(output_file, 'w') as output:
        # Write the header of the output file
        output.write("gene_id\ttranscript_id\n")
        for gene_id, info in gene_info.items():
            # Write the informations obtained into the input file (or a dot if there are no information)
            transcript_id = info["transcript_id"] if info["transcript_id"] else "."
            protein_id = info["protein_id"] if info["protein_id"] else "."
            output.write(f"{gene_id}\t{transcript_id}\n")

if __name__ == "__main__":
    # Arguments for run the script in a command line
    parser = argparse.ArgumentParser(description = "Get the correspondance \
    between genes ID and transcript ID.")
    parser.add_argument('-i', '--input', help = "The input GFF file to \
    extract the genes ID and the transcript ID.")
    parser.add_argument('-o', '--output', help = "The output file which \
    the IDs were strored.")
    args = parser.parse_args()

    # Executue the function
    correspondance_transcript_genes(args.input, args.output)