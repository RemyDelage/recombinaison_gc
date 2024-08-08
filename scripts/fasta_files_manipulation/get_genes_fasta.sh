#!/bin/bash
#SBATCH --job-name=bedtools_getfasta_genes
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --output=bedtools_getfasta_20240320.out

# Initialisation environment
. /local/env/envbedtools-2.27.1.sh

# Variables
data_dir="/home/genouest/cnrs_umr6553/rdelage/data/"
results_dir="/home/genouest/cnrs_umr6553/rdelage/results/genes_fasta/"
input_gff=$(find "${data_dir}" -name "*.gff")
log_file="${data_dir}get_fasta.log"
test_fasta_file="${data_dir}GCF_000001735.4_Arabidopsis_thaliana_genomic.fa"
test_gff_file="${data_dir}GCF_000001735.4_Arabidopsis_thaliana_genomic.gff"

# Extract sequence identifiers from the FASTA file
fasta_identifiers=$(grep "^>" "$test_fasta_file" | sed 's/>//')


# Read each gff file
for gff_file in $test_gff_file; do
	# Get the files names
	file_name=$(basename "$gff_file" .gff)

	# Get the accessions and the species names
	species_id=$(echo "$file_name" | cut -d '_' -f 1-4)
	
	#Create the BED file variable
	bed_file="${results_dir}${species_id}_genes.bed"

	# read each line of the gff file
	while IFS=$'\t' read -r chromosome source feature start end score strand frame attributes; do
		# Check if the feature is a gene
		if [[ "$feature" == "gene" ]]; then
			# Extract the gene id
			gene_id=$(echo "$attributes" | grep -Po 'Name=[^;]+' | sed 's/Name=//')
			# Create the name of the output files
			output_file="${results_dir}${species_id}_${gene_id}.faa"

			# Create the bed file
			awk -F '[\t;]' '/^#/{next} $3 == "gene" {split($11, gene_id, "="); print $1, $4, $5, $7, gene_id[2]}' "$gff_file" > "$bed_file"
	

			# Run bcftools
			bedtools getfasta -fi "${data_dir}${species_id}_genomic.fa" -bed "$bed_file" -name -s >> "$output_file"

			# Print the name of the created file in log file
			echo "File $output_file created" >> "$log_file"
                fi
        done < "$gff_file"
done
