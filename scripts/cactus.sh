#!/bin/bash
#SBATCH --job-name=cactus
#SBATCH --cpus-per-task=8
#SBATCH --mem=350G
#SBATCH --output=cactus_all_genomes_20240405.log


#######################################################################
# Sequences alignments on reference genome using (Progressive) Cactus #
#######################################################################

# Variables definition
data_dir="/home/genouest/cnrs_umr6553/rdelage/data/"
out_dir="/home/genouest/cnrs_umr6553/rdelage/results/03_cactus/"
species_dir="arabidopsis/"
jobStorePath="${out_dir}tmp4"
inputSeqFile="${data_dir}GenomicSequencesFile_arabidopsis.txt"
sif_file="/scratch/rdelage/cactus.sif"
outputHalFile="arabidopsis_all_genomes_alignment.hal"

# Create the output directory
mkdir -p "${out_dir}/cactus_results"

# Set the maximum memory available on the cluster (recommanded by documentation : https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/progressive.md#running-on-a-cluster)
max_memory="3T"

# Script execution (with Singularity)
singularity exec --bind "$data_dir":/mnt --bind "$out_dir" --pwd "$data_dir" "$sif_file" cactus --consCores 2 --maxMemory "$max_memory" "$jobStorePath" "$inputSeqFile" "$out_dir$outputHalFile"
