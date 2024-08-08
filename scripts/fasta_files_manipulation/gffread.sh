#!/bin/bash
#SBATCH --job-name=gffread_filter
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --output=gffread_filter_test_20240325.out

# Initialisation environment
. /local/env/envgffread-0.12.7.sh

# Variables
data_dir="/home/genouest/cnrs_umr6553/rdelage/home/data"
test_gff_input="${data_dir}/GCF_000001735.4_Arabidopsis_thaliana_genomic.gff"
test_gff_output="${data_dir}/GCF_000001735.4_Arabidopsis_thaliana_genomic_filtered.gff"

#Command line
gffread -E -C --keep-genes $test_gff_input -o $test_gff_output
