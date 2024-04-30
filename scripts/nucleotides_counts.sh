#!/bin/bash
#SBATCH --job-name=nucleotide_counts
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --output=nucleotides_counts_20240429.log

# Initialisation environment
. /local/env/envpython-3.9.5.sh

# Variables
input_dir="/home/genouest/cnrs_umr6553/rdelage/results/05_pybedtools_intersect/"
output_dir="/home/genouest/cnrs_umr6553/rdelage/results/06_est-sfs/"

# Script execution
for file in "${input_dir}*.bed.gz"; do
	filename=$(basename "$file" .bed.gz)
	python ./nucleotides_counts_for_est_sfs.py --input_dir "$input_dir" --output_dir "$output_dir"
done
