#!/bin/bash

#########################################################################################
# This script modify the fna.gz files into fasta.gz files for use them with orthofinder #
#########################################################################################

# Go to the repository of the fna.gz files
cd /home/redelage/internship/datas

# Treatment of the files
for file in *.fna.gz; do
	# Decompression of the file
	gunzip "$file"
	# Rename the file
	mv "${file%.fna.gz}.fna" "${file%.fna.gz}.fasta"
       	# Compression of the fasta files
	gzip "${file%.fna.gz}.fasta"
#End of the loop
done

#Print end script message
echo "Complete !"
