#!/bin/bash
#BATCH --job-name=gffutils
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --output=gffutils_20240322.out

date
# Initailisation environment
. /local/env/envconda.sh
conda activate gffutils-0.12-1

# Run the script
python gffutils_filter_gff.py

date
