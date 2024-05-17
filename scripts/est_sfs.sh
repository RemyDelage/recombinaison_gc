#!/bin/bash
#SBATCH --job-name=est_sfs
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G
#SBATCH --output=est_sfs_20240502.log

# Variables
files_dir="/home/genouest/cnrs_umr6553/rdelage/results/06_est-sfs/"
data_files="${files_dir}*_filtered_snp.txt"
config_file="${files_dir}config_file.txt"
seed_file="${files_dir}seedfile.txt"


for file in ${data_files}; do
        echo "File : ${file}"
        filename=$(basename "$file")
        est-sfs "$config_file" "$data_file" "$seed_file" "${file%_snp.txt}_output_file_sfs.txt" "${file%_snp.txt}_output_file_pvalues.txt"
done
