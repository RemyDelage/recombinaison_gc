####################################################
# SNP count, alignment analysis and data preparation 
# for EST-SFS (Keigthley & Jackson, 2019)
####################################################

setwd("/home/genouest/cnrs_umr6553/rdelage/results/04_wgabed/04d_Populus_tremula/")

# ---------------------------- #
# Environment initialization
# ---------------------------- #

# Package for command line parser
library(optparse)


# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Parse the arguments
option_list <- list(
  make_option(c("-s", "--species"), type="character", default=NULL,
              help="Species name", metavar="character"),
  make_option(c("-g", "--species_genus"), type="character", default=NULL,
              help="Species genus", metavar="character"),
  make_option(c("-c", "--chromosome_num"), type="integer", default=NULL,
              help="Chromosome number", metavar="integer"),
  make_option(c("-i", "--chromosome_id"), type="character", default=NULL,
              help="Chromosome ID", metavar="character"),
  make_option(c("-d", "--counts_dir"), type="character", default=NULL,
              help="Counts files directory", metavar="character"),
  make_option(c("-f", "--counts_file"), type="character", default=NULL,
              help="Counts file", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser, args)

# Assign variables from options
species <- opt$species
species_genus <- opt$species_genus
chromosome_num <- opt$chromosome_num
chromosome_id <- opt$chromosome_id
counts_dir <- opt$counts_dir
counts_file <- opt$counts_file

# Debugging: print the variables to check if they are correctly assigned
cat("species: ", species, "\n")
cat("species_genus: ", species_genus, "\n")
cat("chromosome_num: ", chromosome_num, "\n")
cat("chromosome_id: ", chromosome_id, "\n")
cat("counts_file: ", counts_file, "\n")

# Verify that all required arguments are provided
if (is.null(species) || is.null(species_genus) || is.null(chromosome_num) || is.null(chromosome_id) || is.null(counts_dir) || is.null(counts_file)) {
  stop("All arguments --species, --species_genus, --chromosome_num, --chromosome_id, --counts_dir and --counts_file are required")
}


#'@title Load necessary packages
#' 
#' @description Allows to load all the necessary packages for process data and estimate gBGC values
#' 
#'
#' @return None

load_packages <- function(){
 library(dplyr)
}

load_packages()

# ---------------------------- #
# Population analysis
# ---------------------------- #


#'@title SNP analysis from population
#' 
#'@description Allows to counts the SNPs from the WGAbed file
#'
#'@param species (character). The name of the studied species (separate with underscore) 
#'@param chromosome_num (integer). The number of the studied chromosome 
#'
#'@return Data frame : snp_counts 
#' A data frame with the counts of each SNPs into the population and the information about the mutation
#' 
#'@example counts_snp_population("Arabidopsis_thaliana", 1)

counts_snp_population <- function(species, chromosome_num, counts_file){
  #wgabed_file = "/home/genouest/cnrs_umr6553/rdelage/results/04_wgabed/04a_Arabidopsis_thaliana/Arabidopsis_thaliana_1001genomes.1.frq.count"
  
  # Read the SNP counts file produced by WGAbed.
  wgabed_file = counts_file
  snp_counts = read.table(wgabed_file, header = F, sep = "\t", skip = 1)
  
  # Update the names of the columns.
  colnames(snp_counts) = c("CHROM","POS","N_ALLELES","N_CHR","{ALLELE:COUNT}_REF","{ALLELE:COUNT}_ALT")
  
  # Split SNP count results for easier manipulation.
  snp_counts$REF = gsub(":[0-9]+", "", snp_counts$`{ALLELE:COUNT}_REF`)
  snp_counts$ALT = gsub(":[0-9]+", "", snp_counts$`{ALLELE:COUNT}_ALT`)
  snp_counts$count_REF = gsub("[ATCG]:", "", snp_counts$`{ALLELE:COUNT}_REF`)
  snp_counts$count_ALT = gsub("[ATCG]:", "", snp_counts$`{ALLELE:COUNT}_ALT`)
  
  # Convert the counts of each alleles into numeric format.
  snp_counts$count_A = as.numeric(ifelse(snp_counts$REF == "A", snp_counts$count_REF, ifelse(snp_counts$ALT == "A", snp_counts$count_ALT, 0)))
  snp_counts$count_T = as.numeric(ifelse(snp_counts$REF == "T", snp_counts$count_REF, ifelse(snp_counts$ALT == "T", snp_counts$count_ALT, 0)))
  snp_counts$count_C = as.numeric(ifelse(snp_counts$REF == "C", snp_counts$count_REF, ifelse(snp_counts$ALT == "C", snp_counts$count_ALT, 0)))
  snp_counts$count_G = as.numeric(ifelse(snp_counts$REF == "G", snp_counts$count_REF, ifelse(snp_counts$ALT == "G", snp_counts$count_ALT, 0)))
  
  # Get the number of individuals of the population.
  snp_counts$n_genomes = rowSums(snp_counts[,c("count_A", "count_T", "count_C", "count_G")])
  
  # Get the direction of the mutation for the SNPs apparitions.
  snp_counts$mut = unlist(lapply(1:nrow(snp_counts), function(x) {paste(snp_counts$REF[x], snp_counts$ALT[x], sep = "->")}))
  head(snp_counts$mut)
  table(snp_counts$mut)

  # Determine if it is a WS mutation (from AT to GC).
  sum(snp_counts$mut %in% c("A->G", "T->G", "A->C", "T->C"))
  # Determine if it is a SW mutation (for GC to AT).
  sum(snp_counts$mut %in% c("G->A", "G->T", "C->A", "C->T"))
  
  # Store the SNP counts data frame for reuse it into other functions.
  snp_counts <<- snp_counts
  return(head(snp_counts))
}



# ---------------------------- #
# Alignments analysis
# ---------------------------- #


#'@title Alignments analysis
#' 
#'@description Recovers alignment data (BED format) to complete SNP count table
#'
#'@param species (character). The name of the studied species (according to the bed file name)
#'@param chromosome_id (character). The identifyer of the chromosome (from the FASTA sequence file) 
#'@param chromosome_num (integer). The number of the studied chromosome 
#'
#'@return Data frame : snp_counts_merged 
#' A data frame with all the information about SNPs (position, frequency, alignments...) 
#' 
#'@example cactus_analysis("arabidopsis", "NC_003070.9", 1)


cactus_analysis <- function(species_genus, chromosome_id, chromsome_num){
  # Read the alignment file (produced by Cactus aligner).
  wgabed = read.table(gzfile(paste0(species_genus, "_", chromosome_id, ".wga.bed.gz")))
  # Upload the column names.
  colnames(wgabed) = c("seq", "start", "end", "strand", "species", "aligned_chrom", "aligned_pos", "sequence", "strands", "score")
  wgabed$CHROM = chromosome_num
  
  # Keep only SNPs.
  wgabed$width = wgabed$end - wgabed$start
  cat("wgabed before SNP filtering \n")
  print(nrow(wgabed))
  wgabed = wgabed[which(wgabed$width == 1),]
  cat("wgabed after SNP filtering \n")
  print(nrow(wgabed))
  
  # Convert the datas from 1-based into 0-based.
  wgabed$POS = wgabed$start + 1
  
  # Obtain the allele of each species at this position (in uppercase).
  wgabed$ref = toupper(unlist(lapply(strsplit(wgabed$sequence, ","), function(x) {x[1]})))
  wgabed$outgroup_1 = toupper(unlist(lapply(strsplit(wgabed$sequence, ","), function(x) {x[2]})))
  wgabed$outgroup_2 = toupper(unlist(lapply(strsplit(wgabed$sequence, ","), function(x) {x[3]})))
  
  # Keep the sequences of the reference species and outgroups if they are SNPs. Replace with NA if not.
  wgabed$ref = unlist(lapply(1:nrow(wgabed), function(x) {ifelse(nchar(wgabed$ref[x]) == 1, wgabed$ref[x], NA)}))
  wgabed$outgroup_1 = unlist(lapply(1:nrow(wgabed), function(x) {ifelse(nchar(wgabed$outgroup_1[x]) == 1, wgabed$outgroup_1[x], NA)}))
  wgabed$outgroup_2 = unlist(lapply(1:nrow(wgabed), function(x) {ifelse(nchar(wgabed$outgroup_2[x]) == 1, wgabed$outgroup_2[x], NA)}))
  
  
  # Check for duplicate lines (created by WGAbed).
  length(wgabed$POS)
  length(unique(wgabed$POS))
  head(wgabed[duplicated(wgabed$POS),])

  # Remove duplicated lines.
  cat("remove duplicates \n")
  wgabed = wgabed[!duplicated(wgabed$POS), ]
  length(wgabed$POS)
  length(unique(wgabed$POS))
  
  # Convert CHROM variable to merge the 2 data frames
  wgabed$CHROM = is.character(wgabed$CHROM)
  snp_counts= is.character(snp_counts$CHROM)
 
  # Merge the previous SNP counts data frame with the alignments data.
  snp_counts_merged = left_join(snp_counts, wgabed[,c("CHROM", "POS", "ref", "outgroup_1", "outgroup_2")])

  # TO FIX : 
  # Error in UseMethod("left_join") :
  #    no applicable method for 'left_join' applied to an object of class "logical"
  # Calls: cactus_analysis -> left_join

  
  # Display the first lines of merged data.
  print(head(snp_counts_merged))
  
  # Store the data frame into the global environment (for use it into other functions).
  snp_counts_merged <<- snp_counts_merged
  
  # Save the table for use it after EST-SFS computations
  write.table(snp_counts_merged, paste0(species_genus, "_", chromosome_id, "_snp_counts_merged_before_est-sfs.tsv"), row.names = F, col.names = F, quote = F, sep = "\t")
  
  return(head(snp_counts_merged))
  }



# ---------------------------- #
# Create EST-SFS input file
# ---------------------------- #

create_est_sfs_file <- function(species, chromosome_num){
  
  ## Convert the outgroups into one-hot encoding.
  # Outgroup 1
  snp_counts_merged$A_outgroup_1 = ifelse(snp_counts_merged$outgroup_1 == "A", 1, 0)
  snp_counts_merged$T_outgroup_1 = ifelse(snp_counts_merged$outgroup_1 == "T", 1, 0)
  snp_counts_merged$C_outgroup_1 = ifelse(snp_counts_merged$outgroup_1 == "C", 1, 0)
  snp_counts_merged$G_outgroup_1 = ifelse(snp_counts_merged$outgroup_1 == "G", 1, 0)
  
  # Outgroup 2
  snp_counts_merged$A_outgroup_2 = ifelse(snp_counts_merged$outgroup_2 == "A", 1, 0)
  snp_counts_merged$T_outgroup_2 = ifelse(snp_counts_merged$outgroup_2 == "T", 1, 0)
  snp_counts_merged$C_outgroup_2 = ifelse(snp_counts_merged$outgroup_2 == "C", 1, 0)
  snp_counts_merged$G_outgroup_2 = ifelse(snp_counts_merged$outgroup_2 == "G", 1, 0)
  
  ## Replace NA with 0
  # Outgroup 1
  snp_counts_merged$A_outgroup_1[which(is.na(snp_counts_merged$A_outgroup_1))] = 0
  snp_counts_merged$T_outgroup_1[which(is.na(snp_counts_merged$T_outgroup_1))] = 0
  snp_counts_merged$C_outgroup_1[which(is.na(snp_counts_merged$C_outgroup_1))] = 0
  snp_counts_merged$G_outgroup_1[which(is.na(snp_counts_merged$G_outgroup_1))] = 0

  # Outgroup 2
  snp_counts_merged$A_outgroup_2[which(is.na(snp_counts_merged$A_outgroup_2))] = 0
  snp_counts_merged$T_outgroup_2[which(is.na(snp_counts_merged$T_outgroup_2))] = 0
  snp_counts_merged$C_outgroup_2[which(is.na(snp_counts_merged$C_outgroup_2))] = 0
  snp_counts_merged$G_outgroup_2[which(is.na(snp_counts_merged$G_outgroup_2))] = 0
  
  
  # Format the alleles counts for create the EST-SFS input file (see EST-SFS (version 2.04) manual, Keigthley, 2019)
  # example :  20,0,0,0        0,0,0,1 0,0,0,1
  countRef = unlist(lapply(1:nrow(snp_counts_merged), function(x) {paste(snp_counts_merged[x,c("count_A", "count_C", "count_G", "count_T")], collapse = ",")}))
  countOut1 = unlist(lapply(1:nrow(snp_counts_merged), function(x) {paste(snp_counts_merged[x,c("A_outgroup_2", "C_outgroup_1", "G_outgroup_1", "T_outgroup_1")], collapse = ",")}))
  countOut2 = unlist(lapply(1:nrow(snp_counts_merged), function(x) {paste(snp_counts_merged[x,c("A_outgroup_2", "T_outgroup_2", "G_outgroup_2", "T_outgroup_2")], collapse = ",")}))

  # Create the data fame
  df = data.frame(ref_species = countRef, outgroup_1 = countOut1, outgroup_2 = countOut2)
  
  # Software constraints : 1 million of SNPs maximum allowed
  n_rows <- nrow(df)
  cat("Number of sites :", n_rows, "\n")
  
  if (n_rows < 1000000) {
    write.table(df, paste0(species, "_chrom_", chromosome_num, "_est-sfs-input.txt"), 
                row.names = F, col.names = F, quote = F, sep = " ")
  } else {
    write.table(df[1:(1000000-1),], paste0(species, "_chrom_", chromosome_num, "_est-sfs-input.txt"), 
                row.names = F, col.names = F, quote = F, sep = " ")
  }
  
  return(head(df))
}

counts_snp_population(species, chromosome_num, counts_file)
cactus_analysis(species_genus, chromosome_id, chromsome_num)
create_est_sfs_file(speices, chromsome_num)

