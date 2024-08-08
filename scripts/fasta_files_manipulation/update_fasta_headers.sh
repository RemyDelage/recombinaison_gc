#!/bin/bash
#SBATCH --job-name=update_fatsa_headers
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --output=update_fasta_headers_arabidopsis_20240325.out

# Initialisation environment
. /local/env/envpython-3.9.5.sh

# Variables
directory="/home/genouest/cnrs_umr6553/rdelage/results/02_genes_fasta/"
gff_files=(${directory}*.gff)
fasta_files=(${directory}*.fa)


# Code
for gff_file in ${gff_files[@]}; do
    filename=$(basename -- "$gff_file")
    prefix="${filename%%_*}"
    output_file="${directory}${prefix}_genes_sequences.fa"
    
    echo "Output file: ${output_file}"

    if [ -f "$output_file" ]; then
        echo "Output file ${output_file} already exists. Skipping..."
        python update_fasta_headers.py "$gff_file" "${fasta_files[@]}" "$output_file"
        echo "${gff_file} treated. ${output_file} created."
    else
        echo "Fasta file not found for ${gff_file}. Skipping..."
    fi
done
