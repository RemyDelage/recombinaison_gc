#!/bin/bash

#######################################################################
# Sequences alignments on reference genome using (Progressive) Cactus #
#######################################################################

# Variables definition
data_dir="/home/genouest/cnrs_umr6553/rdelage/data/"
jobStorePath="$data_dir""tmp"
inputSeqFile="$data_dir""GenomicSequencesFile_arabidopsis.txt"
out_dir="/home/genouest/cnrs_umr6553/rdelage/results/02_cactus/"
outputHalFile="arabidopsis_alignment.hal"
outputMafFile="arabidopsis_alignment.maf"


# Script execution (with Singularity)
singularity exec --bind "$data_dir":/mnt --bind "$out_dir":/mnt --pwd $PWD:/mnt cactus.sif cactus --consCores 16 "$jobStorePath" "$inputSeqFile" "$out_dir$outputHalFile"

# Convert hal files to maf files for make it readable by many tools (HAL Tool ; use with Singularity)
#singularity exec cactus.sif hal2maf "$out_dir$outputHalFile" "$out_dir$outputMafFile"
