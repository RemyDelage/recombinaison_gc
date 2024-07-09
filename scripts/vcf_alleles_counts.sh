#!/bin/bash
#SBATCH --job-name=vcftools_counts
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G
#SBATCH --output=vcftools_counts_sorghum_20240709.log

# Environment
. /local/env/envvcftools-0.1.16.sh

# Variables
directory="/home/genouest/cnrs_umr6553/rdelage/data/"
output_dir="/home/genouest/cnrs_umr6553/rdelage/results/vcftools/"
input_vcf="${directory}Sorghum_bicolor_Lozano2021.pop_sfs.no_indels.recode.vcf.gz"

# Script execution
vcftools --gzvcf "$input_vcf" --counts --out "${output_dir}Sorghum_bicolor_snp_vcf"
vcftools --gzvcf "$input_vcf" --freq --out "${output_dir}Sorghum_bicolor_freq"
