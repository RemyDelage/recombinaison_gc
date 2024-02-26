#!/bin/bash


####################################################
## De-novo orthologs search with OrthoFinder 2    ##
####################################################

# Intitialization environment
. /local/env/envorthofinder-2.5.2.sh

# Variables
out_dir="/home/genouest/cnrs_umr6553/rdelage/results/01_orthofinder"
in_dir="/home/genouest/cnrs_umr6553/rdelage/data"

# Uncompress the fasta files
for file in "$in_dir"/*.fasta.gz; do
	gunzip "$file"
done

# Command line
orthofinder -o "$out_dir" -f "$in_dir"

# Compress the fasta file at the end of the orthofinder execution
for file in "$in_dir"/*.fasta; do
        gzip "$file"
done
