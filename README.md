# recombination_gc

The aim of this project is to compare the strength of gBGC (GC-biaised gene conversion) between two groups of Angiosperm species (flowering plants) to test the hypothseis according to which gBGC is stronger in monocotyledonous species than in dicotyledonous species.

The method that we use to estimate the gBGC needs to requires a major step, which is the polarization of SNPs (Single Nucleotide Polymorphism). For this we need to go through several stages.

## Get NCBI data

This first step allows to retrieve data at genomic and proteic levels for the studied species and for 2 outgroups species. These outgroups are needed to identify ancestral alleles (present in all species genomes) and derived alleles (appearing after speciation events). The two critera of selection of the outgroups are : 
 * Must be fairly close : not too much divergence between sequences
 * Must be sufficiently remote : polymorphyism must be observed

The script allowed to download directly data from the _RefSeq_ or _GenBank_ databases (NCBI) is the following : 

``` python get_ncbi_data.py -i accession_number -n (or -p)``` 

The parameters are :
 * __-i__ : The accession number
 * __-n__ : Download genomic (nucleotidic) files
 * __-p__ : Download proteic files

By default, the script allows to download genomic data. If they are available, two kind of data were downloaded : 
 * The sequences file (FASTA format)
 * The annotation file (GFF format)

Example : ```python get_ncbi_data.py -i GCF_000001735.4 -n```

## Alleles counts (VCFtools)

For the SNPs polarization, we have to determine which are the major and the minor allele in the population. From VCF file, the __counts__ function from the _VCFtools_ software (Danecek _et al._, 2011) allows to count the occurrences of the two alleles at each SNPs sites.
The user have to run the following command :

``` vcftools --counts --gzvcf vcf_file.vcf.gz --chr chrom_num --min-alleles 2 --max-alleles 2 --out out_path/output_filename```

The parameters are :
 * __--counts__ : Specify you wants to use the **counts** function 
 * __--gzvcf__ : The compressed VCF input file containing the SNPs position and the different alleles
 * __--chr__ : The chromosome on which the counts will be done
 * __--min-alleles__ : Filter the minimal alleles number for the counts (must be 2 for the SNPs)
 * __--max-alleles__ : Filter the maximal alleles number for the counts (must be 2 for the SNPs)
 * __--out__ : Define the path where the counts file will be stored and the prefix file name (the suffixe file name will automatically be added : __.frq.counts__)

Example :
``` vcftools --counts --gzvcf Populus_tremula_Liu2022.pop_sfs.no_indels.recode.vcf.gz --chr "chr1" --min-alleles 2 --max-alleles 2 --out "/results/vcftools/Populus_tremula_Liu2022.pop_sfs.1"```

This script can be run on a computing cluster as a SLURM script (on a computing node) thanks to the **vcftools_counts.sh** shell script :

``` sbatch vcftools_counts.sh -o log_file.log ```

Various parameters must be modified according to the working directory of the user (see scripts directory).

## Genomes alignments (Cactus)

To determine which of the two possible alleles at each polymorphic site is ancestral and which is derived, it is necessary to align the genome sequences of the outgroups with the reference sequence of the species under study.

The **whole genome** alignments are allowed by the *Cactus* aligner (Armstrong *et al.*, 2020).

For use *Cactus*, the user must to type this command line : 

``` cactus --consCores 2 --maxMemory "$max_memory" "$jobStorePath" "$inputSeqFile" "$out_dir$outputHalFile" ```

The parameters are the following : 
 * __jobStorePath__ : Path where temporary files generated by Cactus will be stored. These files will be deleted when execution is complete
 * __inputSeqFile__: file containing the paths to the genomic FASTA files for each species to be aligned (the reference species is
specified by an asterisk)
* __outputHalFile__ : output file which contains the alignment of each site.

Example : 

``` cactus --consCores 2 --maxMemory "3T" "home/results/03_cactus/tmp" "/home/data/GenomicSequencesFile_arabidopsis.txt" "home/results/03_cactus/arabidopsis_all_genomes_alignment.hal" ```

The HAL files produced by *Cactus* are too complicated to use in other tools. They can be converted as MAF files thanks to the *hal2maf* option from *HAL* package (Hickey *et al.*,2013), directly included in *Cactus* thanks to this command line :

``` cactus-hal2maf --refGenome "$referenceGenome" --noAncestors --chunkSize 500000 --batchCores 2 --filterGapCausingDupes "$jobStorePath" "$out_dir$outputHalFile" "$out_dir$outputMafFile" ```

All the parameters information can be found in the official *Cactus* documentation (https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/progressive.md).

Exemple :

```cactus-hal2maf --refGenome GCF_000001735.4_Arabidopsis_thaliana --noAncestors --chunkSize 500000 --batchCores 2 --filterGapCausingDupes "home/results/03_cactus/tmp" "home/results/03_cactus/arabidopsis_all_genomes_alignment.hal" "home/results/03_cactus/arabidopsis_all_genomes_alignment.maf"```

If the user works on a cluster, he/she can run the program with SLURM (on a computing node) thanks to the **cactus.sh** shell script. This file must be modified according to the working environment of the user (see scripts directory). After the modification, the simple command to run is : ```sbatch cactus.sh -o log_file.log```.

## SNPs Polarization (est-sfs)
SNP polarization involves determining which of the two alleles found at this site is the ancestral allele (present in all species) and which is the derived allele (appearing in the species genome after speciation). This step is allowed by the *est-sfs* software (Keigthley & Jackson, 2019).

### Preprocessing
For use *est-sfs* software, data must be formatted to match valid input data (see *est-sfs* documentation). These data can be formatted using the script 

```complete_pipeline_for_est_sfs.R```.

Many parameters must be entered : 

 * **-s** : The species name (genus name and species name must be separate by underscore). Example : *Sorghum_bicolor* ;
 * **-g** : The genus name (must be started by lowercase letter). Example : *sorghum* ;
 * **-c** : The chromosome number processed. Example : *1* ;
 * **-i** : The identifyer of chromosome processed. Example : *NC_012870.2* ;
 * **-d** : The directory of the counts files (i.e. VCFtools results directory) ;
 * **-f** : The counts file of the chromosome (i.e. the file produced by VCFtools).

Example of use : 

```Rscript complete_pipeline_for_est_sfs.R -s "Sorghum_bicolor" -g "sorghum" -c 1 -i "NC_012870.2" -d "/home/results/vcftools/Sorghum_bicolor/" -f "Sorghum_bicolor_Lozano2021.1.frq.count"```


This script also needs the *R* packages **optparse** and **dplyr** to be used. If they are not installed on your environment, you can type the command ```install.packages("optparse")``` and ```install.packages("dplyr")``` in *R* for install them on your device. If you using conda environments, you can also install these packages with the commands ```conda install conda-forge::r-optparse``` and ```conda install conda-forge::r-dplyr```.

If you wants to execute this script on many chromosome, you can use the bash script 

```complete_pipeline_for_est_sfs.sh```.

You must to change the parameters inside the script before run it (they are clearely identifiable). The part *get the chromosomes IDs* not needed to be changed. In that script, each loop lap allows to process one chromosome. You need to adapt the number of loop turns to the number of chromosomes you wish to process.

You can also run this script on a computing cluster as a SLURM script (on a computing node) thanks to the command :

``` sbatch complete_pipeline_for_est_sfs.sh -o log_file.log ```.

### *est-sfs* (Keigthley & Jackson, 2019)

Once the data are correctly formated, the SNPs polarization strictly speaking can be done thanks to the *est-sfs* software. This code is written in *C* language and takes in input 3 kind of files : 

 * **data-input-file.txt** : The formatted data (obtained by the previous script)
 * **config_file.txt** : A file where the user provides information on the parameters used (the number of groups, the substitution model and the number of random runs generated by the algorithm (nrandom)).
 * **seedfile.txt** : Not usefull but mandatory.

Many scripts must be into the same directory for using the software (must not be modified) : **est-sfs**, **est-sfs.c**, **Makefile**, **routine-library.c** and **routine-library.h**.

The command use for run *est-sfs* is contain on a bash script. The default command is : 

```est-sfs "config_file.txt" "$data-input-file.txt" "seedfile.txt" "output_file_sfs.txt" "output_file_pvalues.txt"```.

The variables can be defined and modified into the bash script for run it for each chromosome. The execution command line become simply :

```./est-sfs.sh```.

If you work on computing cluster, you can run it into a SLURM script by using : 

```sbatch est-sfs.sh -o log_file.log ```.

In that script, each loop lap allows to process one chromosome.