# recombination_gc

The method that we use to estimate the gBGC (GC-biaised gene conversion) needs to requires a major step, which is the polarization of SNPs. For that, we need to For this we need to go through several stages

## Get NCBI data

This first step allows to retrieve data at genomic and proteic levels for the studied species and for 2 outgroups species. These outgroups are needed to identify ancestral alleles (present in all species genomes) and derived alleles (appearing after speciation events). The two critera of selection of the outgroups are : 
 * Must be fairly close : not too much divergence between sequences ;
 * Must be sufficiently remote : polymorphyism must be observed.

The script allowed to download directly data from the _RefSeq_ or _GenBank_ databases (NCBI) is the following : 

``` python get_ncbi_data.py -i __accession_number__ -n (or -p)``` 

The parameters are :
 * __-i__ : The accession number
 * __-n__ : Download genomic (nucleotidic) files
 * __-p__ : Download proteic files

By default, the script allows to download genomic data. If they are available, two kind of data were downloaded : 
 * The sequences file (FASTA format)
 * The annotation file (GFF format)
