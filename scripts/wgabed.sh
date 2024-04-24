#!/bin/bash
#SBATCH --job-name=wgabed
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G
#SBATCH --output=wgabed_20240423.log

# Load Python 2 environment for use this script
. /local/env/envpython-2.7.15.sh

# Variables
data_dir="/home/genouest/cnrs_umr6553/rdelage/results/03_cactus/cactus_results/April18/"
results_dir="/home/genouest/cnrs_umr6553/rdelage/results/04_wgabed/"
scripts_dir="/home/genouest/cnrs_umr6553/rdelage/scripts/05_wgabed/WGAbed/"
fasta_file_ref="/home/genouest/cnrs_umr6553/rdelage/results/02_genes_fasta/02a_gffutils/GCF_000001735.4_Arabidopsis_thaliana_genomic.fa"
input_file=$data_dir"arabidopsis_alignment.maf.gz"
uncompressed_input_file=$data_dir"arabidopsis_genomes_alignment.maf"
updated_input_file=$data_dir"arabidopsis_genomes_alignment.maf.gz"
species_names_file=$scripts_dir"data/species_list.txt"

# Chromosomes ID list
chromosomes=($(grep -n ">" "$fasta_file_ref" | awk -F ">" '{print $2}' | awk '{print $1}'))

# Modify the name of the species of the MAF file
gunzip -c "$input_file" > "$uncompressed_input_file"
sed -i 's/GCF_000001735.4_Arabidopsis_thaliana/GCF_000001735_Arabidopsis_thaliana/g' "$uncompressed_input_file"
sed -i 's/GCF_000004255.2_Arabidopsis_lyrata/GCF_000004255_Arabidopsis_lyrata/g' "$uncompressed_input_file"
sed -i 's/GCF_000309985.2_Brassica_rapa/GCF_000309985_Brassica_rapa/g' "$uncompressed_input_file"
gzip -f "$uncompressed_input_file" 

# Cheks if the bgzip package is installed; install it otherwize
if ! command -v bgzip &> /dev/null; then
    pip install bgzip
fi

# Script execution (for each chromosome)
for i in "${chromosomes[@]}"; do
	python "$scripts_dir"maf_to_bed_v2.py -i "$updated_input_file" -r GCF_000001735_Arabidopsis_thaliana -c "$i" -s "$species_names_file" | sort -k1,1 -k2,2n | bgzip -c > "$results_dir"arabidopsis_"$i".wga.bed.gz
done
