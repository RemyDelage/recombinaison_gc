#!/bin/bash
#SBATCH --job-name=cactus
#SBATCH --cpus-per-task=8
#SBATCH --mem=350G
#SBATCH --output=cactus_hal2maf_20240417.log


#######################################################################
# Sequences alignments on reference genome using (Progressive) Cactus #
#######################################################################

# Variables definition
data_dir="/home/genouest/cnrs_umr6553/rdelage/data/"
out_dir="/home/genouest/cnrs_umr6553/rdelage/results/03_cactus/"
species_dir="arabidopsis/"
jobStorePath="${out_dir}tmp"
inputSeqFile="${data_dir}GenomicSequencesFile_arabidopsis.txt"
sif_file="/scratch/rdelage/cactus.sif"
outputHalFile="arabidopsis_all_genomes_alignment.hal"
outputMafFile="arabidopsis_genomes_alignment.maf.gz"

# Create the output directory
mkdir -p "${out_dir}/cactus_results"

# Set the maximum memory available on the cluster (recommanded by documentation : https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/progressive.md#running-on-a-cluster)
max_memory="3T"

# cactus execution (with Singularity)
singularity exec --bind "$data_dir":/mnt --bind "$out_dir" --pwd "$data_dir" "$sif_file" cactus --consCores 2 --maxMemory "$max_memory" "$jobStorePath" "$inputSeqFile" "$out_dir$outputHalFile"

# cactus_hal2maf execution (convert hal file into maf file)
singularity exec --bind "$data_dir":/mnt --bind "$out_dir" --pwd "$data_dir" "$sif_file" cactus-hal2maf --refGenome GCF_000001735.4_Arabidopsis_thaliana --noAncestors --chunkSize 500000 --batchCores 2 --filterGapCausingDupes "$jobStorePath" "$out_dir$outputHalFile" "$out_dir$outputMafFile"

# Alignment inspection
singularity exec "$sif_file" halStats "$outputHalFile"
