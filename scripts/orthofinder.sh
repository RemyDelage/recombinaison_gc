#!/bin/bash


####################################################
## De-novo orthologs search with OrthoFinder 2    ##
####################################################

# Intitialization environment
. /local/env/envorthofinder-2.5.2.sh
. /local/env/envpython-3.9.5.sh

# Variables
out_dir="/home/genouest/cnrs_umr6553/rdelage/results/01_orthofinder/arabidopsis"
in_dir="/home/genouest/cnrs_umr6553/rdelage/data"

# Uncompress the fasta files
for file in "$in_dir"/*.faa.gz; do
	gunzip "$file"
done

# Use Orthofinder on the primary transcripts only
for file in *faa ; do python /home/genouest/cnrs_umr6553/rdelage/scripts/primary_transcript.py $f ; done

# Command line
orthofinder -t 16 -a 16 -S diamond -M dendroblast -o "$out_dir" -f "$in_dir"

# Compress the files at the end of the script execution
for file in "$in_dir"/*.faa; do
	gzip "$file"
done
