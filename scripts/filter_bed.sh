#!/bin/bash
#SBATCH --job-name=pybedtools
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G
#SBATCH --output=pybedtools_20240429.log

#Variables
bed_path="/home/genouest/cnrs_umr6553/rdelage/results/04_wgabed/"
vcf_file="./Arabidopsis_thaliana_1001genomes.pop.vcf"
output_path="/home/genouest/cnrs_umr6553/rdelage/results/05_pybedtools_intersect"

python ./filter_bed.py --input_path "$bed_path" --vcf "$vcf_file" --output_path "$output_path"
