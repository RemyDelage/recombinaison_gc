#!/bin/bash
#SBATCH --job-name=preprocess_est-sfs
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
#SBATCH --output=preprocessing_est_sfs_populus_20240716.log

###########################################################
# 		Pre-processing for EST-SFS		  # 
# Exection of the complete_pipeline_for_est_sfs.R script  #
###########################################################

# Loading environment
. /local/env/envconda.sh
conda activate dplyr

# Variables

data_dir="/home/genouest/cnrs_umr6553/rdelage/results/04_wgabed/04d_Populus_tremula/"
species="Populus_tremula"
species_genus="populus"
counts_dir="/home/genouest/cnrs_umr6553/rdelage/results/vcftools/Populus_tremula/"

# Get the chromosomes IDs
chromosomes_ids=()
files=($(ls "${data_dir}"*.bed.gz))
for file in "${files[@]}"; do 
	chromosomes_ids+=($(basename "$file" | grep -oP 'chr.*(?=\.wga\.bed\.gz)')) 
done
echo "Chromosomes IDs : ${chromosomes_ids[@]}"

# Script execution
for i in {1..19}; do
	date
	chromosome_id=${chromosomes_ids[$i-1]}
	Rscript complete_pipeline_for_est_sfs.R -s "${species}" -g "${species_genus}" -c ${i} -i "${chromosome_id}" -d "${counts_dir}" -f "${counts_dir}Populus_tremula_Liu2022.${i}.frq.count"
	date
done 
