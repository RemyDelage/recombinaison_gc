#!/bin/bash
#SBATCH --job-name=vcftools_counts
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G
#SBATCH --output=vcftools_counts_20240529.log

# Variables
directory="/home/genouest/cnrs_umr6553/rdelage/data/"
output_dir="/home/genouest/cnrs_umr6553/rdelage/results/vcftools/"
input_vcf="${directory}Arabidopsis_thaliana_1001genomes.pop_sfs.no_indels.recode.vcf.gz"



# Script execution
vcftools --gzvcf "$input_vcf" --counts --out "${output_dir}Arabidopsis_thaliana_snp_vcf"
vcftools --gzvcf "$input_vcf" --freq --out "${output_dir}Arabidopsis_thaliana_freq"
