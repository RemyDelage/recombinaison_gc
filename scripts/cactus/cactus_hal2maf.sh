#!/bin/bash
#SBATCH --job-name=hal2maf
#SBATCH --cpus-per-task=8
#SBATCH --mem=350G
#SBATCH --output=cactus_hal2maf_sorghum_20240712.log


#######################################################################
# Sequences alignments on reference genome using (Progressive) Cactus #
#######################################################################

# Variables definition
data_dir="/home/genouest/cnrs_umr6553/rdelage/data/"
out_dir="/home/genouest/cnrs_umr6553/rdelage/results/03_cactus/cactus_results/July09/"
species_dir="sorghum/"
jobStorePath="${out_dir}tmp"
inputSeqFile="${data_dir}GenomicSequencesFile_sorghum.txt"
sif_file="/scratch/rdelage/cactus.sif"
outputHalFile="sorghum_all_genomes_alignment.hal"
outputMafFile="sorghum_genomes_alignment.maf.gz"

# Set the maximum memory available on the cluster (recommanded by documentation : https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/progressive.md#running-on-a-cluster)
max_memory="3T"

# cactus_hal2maf execution (convert hal file into maf file)
singularity exec --bind "$data_dir":/mnt --bind "$out_dir" --pwd "$data_dir" "$sif_file" cactus-hal2maf --refGenome GCF_000003195.3_Sorghum_bicolor --noAncestors --chunkSize 500000 --batchCores 2 --filterGapCausingDupes "$jobStorePath" "$out_dir$outputHalFile" "$out_dir$outputMafFile"

# Alignment inspection
singularity exec "$sif_file" halStats "$outputHalFile"
