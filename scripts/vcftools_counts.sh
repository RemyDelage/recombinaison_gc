#!/bin/bash
#SBATCH --job-name=vcftools_counts
#SBATCH --mem=150G
#SBATCH --cpus-per-task=8


#######################################################
# Count alleles for each SNP position from a VCF file #
#          VCFtools (Danecek et al., 2011)            #
#######################################################

# Environment
. /local/env/envvcftools-0.1.16.sh

# Variables
in_dir="/home/genouest/cnrs_umr6553/rdelage/data/" # Set working directory
out_dir="/home/genouest/cnrs_umr6553/rdelage/results/04_wgabed/04b_Oryza_sativa/" # Set output directory
vcf_file="${in_dir}Oryza_sativa_Wang2018.pop_sfs.no_indels.recode.vcf.gz" # Specify the VCF file

# Loop for do the computations on each chromosome (one after the other)
for i in {1..12} ; do
	vcftools --counts --gzvcf "$vcf_file" --chr "$i" --min-alleles 2 --max-alleles 2 --out "Oryza_sativa_Wang2018.${i}"
done
