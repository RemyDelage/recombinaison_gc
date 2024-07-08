### Setting working directory (find a way to automatise this) ###

setwd("/home/genouest/cnrs_umr6553/rdelage/test")

#'@title Load necessary packages
#' 
#' @description Allows to load all the necessary packages for process data and estimate gBGC values
#' 
#'
#' @return None

load_packages <- function(){
  library(GenomicRanges)
  library(dplyr)
  library(ggplot2)
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

counts_snp_population <- function(species, chromosome_num){
  #wgabed_file = "/home/genouest/cnrs_umr6553/rdelage/results/04_wgabed/04a_Arabidopsis_thaliana/Arabidopsis_thaliana_1001genomes.1.frq.count"
  wgabed_file = paste0(species,"_1001genomes.", chromosome_num, ".frq.count")
  snp_counts = read.table(wgabed_file, header = F, sep = "\t", skip = 1)
  head(snp_counts)
  
  colnames(snp_counts) = c("CHROM","POS","N_ALLELES","N_CHR","{ALLELE:COUNT}_REF","{ALLELE:COUNT}_ALT")
  head(snp_counts)
  
  snp_counts$REF = gsub(":[0-9]+", "", snp_counts$`{ALLELE:COUNT}_REF`)
  snp_counts$ALT = gsub(":[0-9]+", "", snp_counts$`{ALLELE:COUNT}_ALT`)
  snp_counts$count_REF = gsub("[ATCG]:", "", snp_counts$`{ALLELE:COUNT}_REF`)
  snp_counts$count_ALT = gsub("[ATCG]:", "", snp_counts$`{ALLELE:COUNT}_ALT`)
  
  head(snp_counts)
  
  snp_counts$count_A = as.numeric(ifelse(snp_counts$REF == "A", snp_counts$count_REF, ifelse(snp_counts$ALT == "A", snp_counts$count_ALT, 0)))
  snp_counts$count_T = as.numeric(ifelse(snp_counts$REF == "T", snp_counts$count_REF, ifelse(snp_counts$ALT == "T", snp_counts$count_ALT, 0)))
  snp_counts$count_C = as.numeric(ifelse(snp_counts$REF == "C", snp_counts$count_REF, ifelse(snp_counts$ALT == "C", snp_counts$count_ALT, 0)))
  snp_counts$count_G = as.numeric(ifelse(snp_counts$REF == "G", snp_counts$count_REF, ifelse(snp_counts$ALT == "G", snp_counts$count_ALT, 0)))
  
  snp_counts$n_genomes = rowSums(snp_counts[,c("count_A", "count_T", "count_C", "count_G")])
  
  snp_counts$mut = unlist(lapply(1:nrow(snp_counts), function(x) {paste(snp_counts$REF[x], snp_counts$ALT[x], sep = "->")}))
  head(snp_counts$mut)
  table(snp_counts$mut)

  # WS
  sum(snp_counts$mut %in% c("A->G", "T->G", "A->C", "T->C"))
  # SW
  sum(snp_counts$mut %in% c("G->A", "G->T", "C->A", "C->T"))
  
  snp_counts <<- snp_counts
  return(snp_counts)
}

counts_snp_population("Arabidopsis_thaliana", 1)

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


cactus_analysis <- function(species, chromosome_id, chromsome_num){
  # Read the alignment file (produced by Cactus aligner)
  wgabed = read.table(gzfile(paste0(species, "_", chromosome_id, ".wga.bed.gz")))
  # Upload the column names
  colnames(wgabed) = c("seq", "start", "end", "strand", "species", "aligned_chrom", "aligned_pos", "sequence", "strands", "score")
  wgabed$CHROM = chromosome_num
  
  # Keep only SNPs
  wgabed$width = wgabed$end - wgabed$start
  cat("wgabed before SNP filtering \n")
  nrow(wgabed)
  wgabed = wgabed[which(wgabed$width == 1),]
  cat("wgabed after SNP filtering \n")
  nrow(wgabed)
  
  # Convert the datas from 1-based into 0-based
  wgabed$POS = wgabed$start + 1
  
  # Obtain the allele of each species at this position (in uppercase)
  wgabed$ref = toupper(unlist(lapply(strsplit(wgabed$sequence, ","), function(x) {x[1]})))
  wgabed$outgroup_1 = toupper(unlist(lapply(strsplit(wgabed$sequence, ","), function(x) {x[2]})))
  wgabed$outgroup_2 = toupper(unlist(lapply(strsplit(wgabed$sequence, ","), function(x) {x[3]})))
  
  # Keep the sequences of the reference species and outgroups if they are SNPs. Replace with NA if not
  wgabed$ref = unlist(lapply(1:nrow(wgabed), function(x) {ifelse(nchar(wgabed$ref[x]) == 1, wgabed$ref[x], NA)}))
  wgabed$outgroup_1 = unlist(lapply(1:nrow(wgabed), function(x) {ifelse(nchar(wgabed$outgroup_1[x]) == 1, wgabed$outgroup_1[x], NA)}))
  wgabed$outgroup_2 = unlist(lapply(1:nrow(wgabed), function(x) {ifelse(nchar(wgabed$outgroup_2[x]) == 1, wgabed$outgroup_2[x], NA)}))
  
  
  # Check for duplicate lines (created by WGAbed)
  length(wgabed$POS)
  length(unique(wgabed$POS))
  head(wgabed[duplicated(wgabed$POS),])

  # Remove duplicated lines
  cat("remove duplicates \n")
  wgabed = wgabed[!duplicated(wgabed$POS), ]
  length(wgabed$POS)
  length(unique(wgabed$POS))
  
  # Merge the previous SNP counts data frame with the alignments data
  snp_counts_merged = left_join(snp_counts, wgabed[,c("CHROM", "POS", "ref", "outgroup_1", "outgroup_2")])
  
  # Display the first lines of merged data
  print(head(snp_counts_merged))
  
  # Store the data frame into the global environment (for use it into other functions)
  snp_counts_merged <<- snp_counts_merged
  
  return(snp_counts_merged)
  }

cactus_analysis("arabidopsis", "NC_003070.9", 1)


# ---------------------------- #
# Create EST-SFS input file
# ---------------------------- #

create_est_sfs_file

# 
# 
# # On convertit les outgroup en one-hot encoding
# snp_counts_merged$A_outgroup_1 = ifelse(snp_counts_merged$outgroup_1 == "A", 1, 0)
# snp_counts_merged$T_outgroup_1 = ifelse(snp_counts_merged$outgroup_1 == "T", 1, 0)
# snp_counts_merged$C_outgroup_1 = ifelse(snp_counts_merged$outgroup_1 == "C", 1, 0)
# snp_counts_merged$G_outgroup_1 = ifelse(snp_counts_merged$outgroup_1 == "G", 1, 0)
# 
# snp_counts_merged$A_outgroup_2 = ifelse(snp_counts_merged$outgroup_2 == "A", 1, 0)
# snp_counts_merged$T_outgroup_2 = ifelse(snp_counts_merged$outgroup_2 == "T", 1, 0)
# snp_counts_merged$C_outgroup_2 = ifelse(snp_counts_merged$outgroup_2 == "C", 1, 0)
# snp_counts_merged$G_outgroup_2 = ifelse(snp_counts_merged$outgroup_2 == "G", 1, 0)
# 
# # Replace NA with 0
# snp_counts_merged$A_outgroup_1[which(is.na(snp_counts_merged$A_outgroup_1))] = 0
# snp_counts_merged$T_outgroup_1[which(is.na(snp_counts_merged$T_outgroup_1))] = 0
# snp_counts_merged$C_outgroup_1[which(is.na(snp_counts_merged$C_outgroup_1))] = 0
# snp_counts_merged$G_outgroup_1[which(is.na(snp_counts_merged$G_outgroup_1))] = 0
# 
# snp_counts_merged$A_outgroup_2[which(is.na(snp_counts_merged$A_outgroup_2))] = 0
# snp_counts_merged$T_outgroup_2[which(is.na(snp_counts_merged$T_outgroup_2))] = 0
# snp_counts_merged$C_outgroup_2[which(is.na(snp_counts_merged$C_outgroup_2))] = 0
# snp_counts_merged$G_outgroup_2[which(is.na(snp_counts_merged$G_outgroup_2))] = 0
# 
# 
# head(snp_counts_merged)

# View(snp_counts_merged[1:100,])
