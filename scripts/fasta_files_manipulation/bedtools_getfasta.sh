#!/bin/bash
#BATCH --job-id=bedtools_getfasta
#BATCH --cpus-per-task=8
#BATCH --mem=50G
#BATCH --output=bedtools_getfasta_20240322.out

# Initialisation environment
. /local/env/envbedtools-2.27.1.sh
. /local/env/envpython-3.9.5.sh

# Variables definitions
genomes_fasta="/home/genouest/cnrs_umr6553/rdelage/data/genomes_fasta/"
gff_files_dir="/home/genouest/cnrs_umr6553/rdelage/results/02_genes_fasta/02a_gffutils/"
gff_files="${gff_files_dir}*.gff"
output_dir="/home/genouest/cnrs_umr6553/rdelage/results/02_genes_fasta/"

#Code
date
for file in $gff_files; do
	echo "File read : ${file}"
	# Get the accession number
	accession=$(echo "$file" | awk -F'/' '{print $NF}' | awk -F'_' '{print $1 "_" $2 "_" $3 "_" $4 "_" $5}')	
	# Read the correct fasta file
	fasta_file="${genomes_fasta}${accession}.fa"
        echo "Checking fasta file : ${fasta_file}"
	# Checks the file readed
	if [ -f "$fasta_file" ]; then
		echo "GFF file readed : ${file}"
		echo "FASTA file readed : ${fasta_file}"

		# Run bcftools getfasta
		bedtools getfasta -fi "$fasta_file" -bed "$file" -fo "${output_dir}/${accession}_genes.fa" -s
		echo "File ${output_dir}/${accession}_genes.fa created"
		# Rename the fasta header according to the gene ID
		python update_fasta_headers.py "$file" "${output_dir}/${accession}_genes.fa"
	else
		echo "Any FASTA file found for this species : ${accession}"
	fi
done
date
