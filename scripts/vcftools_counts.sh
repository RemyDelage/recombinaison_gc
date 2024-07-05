#!/bin/bash
#SBATCH --job-name=vcftools_counts
#SBATCH --mem=150G
#SBATCH --cpus-per-task=8
#SBATCH --output=vcftools_oryza_20240607.log

# Environment
. /local/env/envvcftools-0.1.16.sh

# Variables
in_dir="/home/genouest/cnrs_umr6553/rdelage/data/"
out_dir="/home/genouest/cnrs_umr6553/rdelage/results/04_wgabed/04b_Oryza_sativa/"
vcf_file="${in_dir}Oryza_sativa_Wang2018.pop_sfs.no_indels.recode.vcf.gz"

for i in {1..12} ; do
	vcftools --counts --gzvcf "$vcf_file" --chr "$i" --min-alleles 2 --max-alleles 2 --out "Oryza_sativa_Wang2018.${i}"
done
