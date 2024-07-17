####################################################
# SNP analysis and gBGC estimations (Glémin, 2022) #
####################################################

setwd("/home/genouest/cnrs_umr6553/rdelage/results/06_est-sfs/06c_Sorghum_bicolor/") # Must be changed

# ---------------------------- #
# Environment initialization
# ---------------------------- #

# Package for command line parser
library(optparse)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Parse the arguments
option_list <- list(
  make_option(c("-s", "--species"), type = "character", default = NULL,
              help = "Species name", metavar = "character"),
  make_option(c("-c", "--chromosome_num"), type = "integer", default = NULL, 
              help = "Chromosome number", metavar = "integer"),
  make_option(c("-f", "--gff_file"), type = "character", default = NULL,
              help = "GFF file of the species", metavar = "character"),
  make_option(c("-g", "--genomic_level"), type = "character", default = NULL,
              help = "Genomic level for the SFS computations (i.e. gene or exon)",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser, args)

# Assign variables from options
species <- opt$species
chromosome_num <- opt$chromosome_num
gff_file <- opt$gff_file
genomic_level <- opt$genomic_level

# Print the variables to check if they are correctly assigned
cat("species: ", species, "\n")
cat("chormosome_num: ", chromosome_num, "\n")
cat("gff_file: ", gff_file, "\n")
cat("genomic_level: ", genomic_level, "\n")

# Check if all required arguments are provided
if (is.null(species) || is.null(chromosome_num) || is.null(gff_file) || is.null(genomic_level)){
  stop("All arguments --species, --chromosome_num, --gff_file, --genomic_level are required.")
}

# ---------------------------- #
# Packages loading
# ---------------------------- #


load_packages <- function(){
  library(GenomicRanges)
  library(dplyr)
  library(ggplot2)
}

load_packages()

# ---------------------------- #
# SNP polarization
# ---------------------------- #

#'@title snp_polarization
#' 
#'@description Enables SNP polarization from EST-SFS results (Keigthley & Jackson, 2019)
#'
#'@param species (character). The name of the studied species (genius name and species name must be separate by underscore)
#'@param chromosome_num (integer). The number of the studied chromosome 
#'
#'@return Data frame : snp_counts_merged 
#' The updated snp_counts_merged data frame with the information about SNP polarization (mutation, ancestal and derivated alleles) 
#' 
#'@example snp_polarization("Arabidopsis_thaliana", 1)


snp_polarization <- function(species, chromosome_num){
  # Import the EST-SFS results table (P-values tables)
  estsfs = read.table(paste0(species, "_chrom_", chromosome_num, "_output_file_pvalues.txt"), header = FALSE, skip = 8, sep = " ")
  estsfs = estsfs[,-8]
  colnames(estsfs) = c("IdxSite", "Code", "P_major_ancestral", "P_trees_A", "P_trees_C", "P_trees_G", "P_trees_T")
  head(estsfs)
  
  
  # Import the snp_counts_merged table (from 'complete_pipeline_for_est_sfs.R' script)
  snp_counts_merged = read.table(paste0("snp_counts_merged_", species, "_chrom_", chromosome_num, ".tsv"), header = TRUE, sep = "\t")
  
  # Add the EST-SFS results into the snp_counts_merged table
  snp_counts_merged$P_major_ancestral = NA
  snp_counts_merged$P_trees_A = NA
  snp_counts_merged$P_trees_C = NA
  snp_counts_merged$P_trees_G = NA
  snp_counts_merged$P_trees_T = NA
  
  snp_counts_merged$P_major_ancestral[estsfs$IdxSite] = estsfs$P_major_ancestral
  snp_counts_merged$P_trees_A[estsfs$IdxSite] = estsfs$P_trees_A
  snp_counts_merged$P_trees_C[estsfs$IdxSite] = estsfs$P_trees_C
  snp_counts_merged$P_trees_G[estsfs$IdxSite] = estsfs$P_trees_G
  snp_counts_merged$P_trees_T[estsfs$IdxSite] = estsfs$P_trees_T
  
  # Display the first rows to check the table
  head(snp_counts_merged)
  
  # Ancestral and derivate allele (test on a specific site)
  cat("Site tested : 10 \n")
  cat(snp_counts_merged[10,])
  
  # Major allele determination
  major = function(x) {
    # Return the major allele for a line x in snp_counts_merged
    alphabet = c("A", "C", "G", "T")
    maj = alphabet[which.max(snp_counts_merged[x, c("count_A", "count_C", "count_G", "count_T")])]
    return(maj)
  }
  cat("Major allele position 10 : \n", major(10))
  
  # Minor allele determination
  minor = function(x) {
    alphabet = c("A", "C", "G", "T")
    s = snp_counts_merged[x, c("count_A", "count_C", "count_G", "count_T")]
    if (all(s == 0)) {
      return(NA) # Add NAs if there is no minor alleles
    } else {
      s[s == 0] = NA # Add NAs if there is no minor alleles
      mino = alphabet[which.min(s)]
      return(mino)
    }
  }
  cat("Minor allele position 10 : \n", minor(10))
  
  # Ancestral allele determination (higher probability among P_trees)
  ancestral = function(x) {
    alphabet = c("A", "C", "G", "T")
    proba_alleles = snp_counts_merged[x, c("P_trees_A", "P_trees_C", "P_trees_G", "P_trees_T")]
    anc = alphabet[which.max(proba_alleles)]
    anc = ifelse(length(anc) == 1, anc, NA)
    return(anc)
  }
  cat("Ancestral allele position 10 : ", ancestral(10), "\n")
  
  # Add the alleles information into the snp_counts_merged table
  snp_counts_merged$Major = unlist(lapply(1:nrow(snp_counts_merged), major))
  snp_counts_merged$Minor = unlist(lapply(1:nrow(snp_counts_merged), minor))
  ## Remove the lines with only one allele (major and minor are the same)
  snp_counts_merged = snp_counts_merged[which(snp_counts_merged$Major != snp_counts_merged$Minor),]
  snp_counts_merged$Ancestral = unlist(lapply(1:nrow(snp_counts_merged), ancestral))
  ## Consider that the other allele is the derivate one
  snp_counts_merged$Derived = ifelse(is.na(snp_counts_merged$Ancestral), NA, ifelse(snp_counts_merged$Ancestral == snp_counts_merged$Major, snp_counts_merged$Minor, snp_counts_merged$Major))
  
  # Put NAs if the allele is not ancestral or derivate
  snp_counts_merged[which(snp_counts_merged$Ancestral != snp_counts_merged$Minor & snp_counts_merged$Ancestral != snp_counts_merged$Major), c("Ancestral", "Derived")] = NA
  
  # Mutations defintitions
  snp_counts_merged$Mutation = NA
  snp_counts_merged$Mutation[which(snp_counts_merged$Ancestral %in% c("A", "T") & snp_counts_merged$Derived %in% c("C", "G"))] = "WS"
  snp_counts_merged$Mutation[which(snp_counts_merged$Ancestral %in% c("C", "G") & snp_counts_merged$Derived %in% c("A", "T"))] = "SW"
  snp_counts_merged$Mutation[which(snp_counts_merged$Ancestral %in% c("A", "T") & snp_counts_merged$Derived %in% c("A", "T"))] = "W"
  snp_counts_merged$Mutation[which(snp_counts_merged$Ancestral %in% c("C", "G") & snp_counts_merged$Derived %in% c("C", "G"))] = "S"
  
  # Save the updated snp_counts_merged table
  write.table(snp_counts_merged, paste0("snp_counts_merged_", species, "_chrom_", chromosome_num, ".tsv"), row.names = FALSE, col.names = TRUE,quote = FALSE, sep = "\t")
  
  # Store the data frame into the global environment (for use it into other functions).
  snp_counts_merged <<- snp_counts_merged
  
  return(head(snp_counts_merged))                    
}

# ------------------------- #
# uSFS computations
# ------------------------- #


#'@title sfs_computations
#' 
#'@description Enables to make unfold SFS (Site frequency spectrum) according to the Bagley et al., 2016 method.
#'The gBGC computation will be done from these uSFS.
#'
#'@param species (character). The name of the studied species (genius name and species name must be separate by underscore)
#'@param chromosome_num (integer). The number of the studied chromosome
#'@param gff_file (character). The annotation file of the studied species. 
#'It is mandatory for determine the genomic level which computations will be made.
#'@param genomic_level (character). Specify the genomic level for the computations (exons or genes)
#'
#'@return none 
#' 
#'@example sfs_computations("Arabidopsis_thaliana", 1, "gff_rho_Arabidopsis_thaliana_1001genomes.rds", "exon")

sfs_computations <- function(species, chromosome_num, gff_file, genomic_level){
  
  # Read the GFF file
  gff = readRDS(gff_file)
  head(gff)
  
  # Create GenomicRanges object for the GFF data
  gff_ranges = makeGRangesFromDataFrame(gff, keep.extra.columns=T)
  head(gff_ranges)
  
  # Create GenomicRanges object for the SNPs data
  snps = snp_counts_merged
  snps$start = snps$POS
  snps$end = snps$POS + 1
  snp_ranges = makeGRangesFromDataFrame(snps)
  head(snp_ranges)
  
  if(genomic_level == "exon"){
    # Subdivide the data for only keep the exon data (with 15 exons maximum)
    exon_ranges = gff_ranges[which(gff_ranges$feature == "exon" & gff_ranges$rank == 1 & gff_ranges$nb_exons < 15)]
    
    # Keep the SNPs founded into the exons
    hits = findOverlaps(snp_ranges, exon_ranges)
    snp_counts_merged = snp_counts_merged[queryHits(hits),]
  
    }else if(genomic_level == "gene"){
    # Subdivide the data for only keep the genes data
    gene_ranges = gff_ranges[which(gff_ranges$feature == "gene")]
    gene_ranges
    
    # Keep the SNPs founded into the genes
    hits = findOverlaps(snp_ranges, gene_ranges)
    snp_counts_merged = snp_counts_merged[queryHits(hits),]
    }
  
  # Filter the alleles to keep only those with a P-trees probability greater than 0.70
  snp_counts_merged$filter = unlist(lapply(1:nrow(snp_counts_merged), function(x) {ifelse(max(snp_counts_merged[x,c("P_trees_A", "P_trees_C", "P_trees_G", "P_trees_T")]) > 0.7, "PASS", "REMOVE")}))
  
  print(table(snp_counts_merged$filter))
  
  snp_counts_merged_filtered_p70 <- snp_counts_merged[which(snp_counts_merged$filter == "PASS"),]
  cat("nrow(snp_counts_merged_filtered_p70) \n", nrow(snp_counts_merged_filtered_p70))
  table(snp_counts_merged_filtered_p70$Derived)
  
  # Counts the number of derived alleles
  snp_counts_merged_filtered_p70$Nb_Derived <- ifelse(snp_counts_merged_filtered_p70$Derived == "A",
                                                      snp_counts_merged_filtered_p70$count_A,
                                                      ifelse(snp_counts_merged_filtered_p70$Derived == "C",
                                                             snp_counts_merged_filtered_p70$count_C,
                                                             ifelse(snp_counts_merged_filtered_p70$Derived == "G",
                                                                    snp_counts_merged_filtered_p70$count_G,
                                                                    snp_counts_merged_filtered_p70$count_T)))
  
  # Compute the total number of alleles in each site
  snp_counts_merged_filtered_p70$Nb_Total_Alleles <- snp_counts_merged_filtered_p70$N_CHR
  
  # Save the snp_counts_merged filtered data
  write.table(snp_counts_merged_filtered_p70, paste0("snp_counts_merged_filtered_p70_", species, "_chrom_", chromosome_num, "_", genomic_level, ".tsv"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  
  # Check the number of each mututaion (WS : AT -> GC, SW : GC -> AT, W : AT <-> AT, S : GC <-> GC)
  table(snp_counts_merged_filtered_p70$Mutation)
  
  # Create 4 data frames according to the mutations types 
  snp_counts_p70_sw <- snp_counts_merged_filtered_p70[which(snp_counts_merged_filtered_p70$Mutation == "SW"),]
  snp_counts_p70_ws <- snp_counts_merged_filtered_p70[which(snp_counts_merged_filtered_p70$Mutation == "WS"),]
  snp_counts_p70_s <- snp_counts_merged_filtered_p70[which(snp_counts_merged_filtered_p70$Mutation == "S"),]
  snp_counts_p70_w <- snp_counts_merged_filtered_p70[which(snp_counts_merged_filtered_p70$Mutation == "W"),]
  
  # Check the number of lines of these data frames (must be the same as the number of mutations)
  cat("nrow(snp_counts_p70_sw) : ", nrow(snp_counts_p70_sw), "\n")
  cat("nrow(snp_counts_p70_ws) : ", nrow(snp_counts_p70_ws), "\n")
  cat("nrow(snp_counts_p70_s) : ", nrow(snp_counts_p70_s), "\n")
  cat("nrow(snp_counts_p70_w) : ", nrow(snp_counts_p70_w), "\n")
  
  ### Making the SFS (Bagley et al., 2016) ####
  # Principle : sub-sampling data to retain an optimal amount of information while removing missing data
  
  subsampling -> function(dataframe){
    # Compute the number of SNPs
    nsnps <- nrow(dataframe)
    # Compute the number of individuals
    nind <- max(dataframe$Nb_Total_Alleles, na.rm = T)/2
    
    cat("Number of SNPs : ", nsnps, "\n")
    cat("Number of individuals : ", nind, "\n")
    
    # We need as input agenotype matrix with missing data and columns are site rows are individuals
    
    # We don't have this matrix (the GT matrix in the vcf) because we already counted the number of alleles
    
    # Pseudo-matrix  with nsnps columns and nind*2 lignes
    md_geno = matrix(NA, nrow = nind*2, ncol = nsnps)
    
    # For each line, create a vector with 
    # n derived alleles (= 1), m ancestral alleles (= 0) et l NAs
    
    for (i in 1:ncol(md_geno)) {
      n = snp_counts_p70_ws$Nb_Derived[i]
      m = snp_counts_p70_ws$Nb_Total[i] - snp_counts_p70_ws$Nb_Derived[i]
      l = nind*2 - snp_counts_p70_ws$Nb_Total[i]
      md_geno[,i] = c(rep(1, n), rep(0, m), rep(NA, l))
    }
    

    # find minimum sample size and store it to use it for each SFS computations
    md_persite <- colSums(is.na(md_geno))
    min_ss <- abs(nind-max(md_persite))
    
    min_ss <<- min_ss
    
    return(min_ss)
  }

  
  making_sfs <- function(mutation, dataframe, chromosome_num, species, genomic_level){
   # Compute the number of SNPs
   nsnps <- nrow(snp_counts)
   # Compute the number of individuals
   nind <- max(snp_counts$Nb_Total_Alleles, na.rm = T)/2
   
   # Create a pseudo-matrix with 'nsnps' columns and 'nind*2' lignes for making the SFS
   md_geno <- matrix(NA, nrow = nind*2, ncol = nsnps)
   
   # Create vector with n derived alleles (= 1), m ancestal alleles (= 0) and l NAs
   for (i in 1:ncol(md_geno)) {
     n <- snp_counts$Nb_Derived[i]
     m <- snp_counts$Nb_Total[i] - snp_counts$Nb_Derived[i]
     l <- nind*2 - snp_counts$Nb_Total[i]
     md_geno[,i] <- c(rep(1, n), rep(0, m), rep(NA, l))
   }
   
   # Sample min_ss genotypes for each site without missing data
   no_ms_list <- lapply(1:ncol(md_geno), function(i) {
     each_site <- md_geno[, i]
     population <- sum(!is.na(each_site))
     sample_size <- min(min_ss, population)
     sample(each_site[!is.na(each_site)], size = sample_size, replace = FALSE)
   })
   
   # Ensure each sampled vector is of the same length
   no_ms_matrix <- do.call(cbind, lapply(no_ms_list, function(x) {
     length_diff <- min_ss - length(x)
     if (length_diff > 0) {
       return(c(x, rep(NA, length_diff)))
     } else {
       return(x)
     }
   }))
   
   # Ensure matrix is numeric
   no_ms_matrix <- apply(no_ms_matrix, 2, as.numeric)
   
   # Remove columns with missing data
   missing_data_indices <- which(is.na(no_ms_matrix), arr.ind = TRUE)
   if (length(missing_data_indices) > 0) {
     no_ms_matrix <- no_ms_matrix[, colSums(is.na(no_ms_matrix)) == 0]
   }
   
   # Create data frame for plotting the SFS
   sfs_df <- data.frame(Nb_sites = seq(0, min_ss),
                        DAF = NA)
   for (i in 1:nrow(sfs_df)) {
     sfs_df$DAF[i] <- sum(colSums(no_ms_matrix) == sfs_df$Nb_sites[i])
   }
   
   # Create the plot for the SFS
   jpeg(paste0("SFS_", mutation_type, "_Chrom_", chromosome_num, "_", species, ".jpg"), width = 1200, height = 700)
   barplot(sfs_df$DAF, xlab="derived allele frequency", ylab="number of sites")
   title(paste0("SFS ", mutation_type, " Chromosome ", chromosome_num, " ", species, " ", genomic_level))
   dev.off()
   cat(paste0("Plot ", "SFS_", mutation_type, "_Chrom_", chromosome_num, "_", species, ".jpg created \n"))
   
   sfs_df <<- paste0("df_", muatation, "_", species, "_", genomic_level)
  }
    
  # Compute the SFS for each mutation type
  
  subsampling(snp_counts_p70_sw) # Allows to get the min_ss (sub-sampling) value
  WS <- compute_sfs("WS", snp_counts_p70_ws, chromosome_num, species, genomic_level)
  SW <- compute_sfs("SW", snp_counts_p70_sw, chromosome_num, species, genomic_level)
  S <- compute_sfs("S", snp_counts_p70_s, chromosome_num, species, genomic_level)
  W <- compute_sfs("W", snp_counts_p70_w, chromosome_num, species, genomic_level)
  
  # Display all the SFS of each mutation
  jpeg(paste0("SFS_ALL_MUT_Chrom_", chromosome_num,"_", species, ".jpg"), width = 1200, height = 700)
  par(mfrow = c(2,2))
  
  barplot(SW$DAF, xlab="derived allele frequency", ylab="number of sites", main = paste0("SW (n = ", nrow(snp_counts_p70_sw), ")"))
  barplot(WS$DAF, xlab="derived allele frequency", ylab="number of sites", main = paste0("WS (n = ", nrow(snp_counts_p70_ws), ")"))
  barplot(W$DAF, xlab="derived allele frequency", ylab="number of sites", main = paste0("W (n = ", nrow(snp_counts_p70_w), ")"))
  barplot(S$DAF, xlab="derived allele frequency", ylab="number of sites", main = paste0("S (n = ", nrow(snp_counts_p70_s), ")"))
  
  par(mfrow = c(1,1))
  dev.off()
  
  # Create data frame for store all the values for establish the SFSs
  df = data.frame(mutation = rep(c("SW", "WS", "W", "S"), each = nrow(SW)),
                  n_sites = c(SW$Nb_sites, WS$Nb_sites, W$Nb_sites, S$Nb_sites),
                  density = c(SW$DAF/sum(SW$DAF), WS$DAF/sum(WS$DAF), W$DAF/sum(W$DAF), S$DAF/sum(S$DAF)))
  
  # Save the data frame
  write.table(df, paste0("df_", species, "_chrom_", chromosome_num, "_", genomic_level, ".tsv"), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
  
  # Plot the SFS for WS and SW mutations on the same plot (to compare them)
  p = ggplot(df[which(df$mutation %in% c("WS", "SW")),], aes(x = n_sites, y = density, colour = mutation, fill = mutation)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme_bw()
  
  # Save these plots
  ggsave(paste0("SFS_WS_SW_Chrom_", chromosome_num, "_", species, ".jpg"), p, width = 14, height =  10)
  
  # Store the necessary data into the global environment for use it into the next function
  
  WS <<- WS
  SW <<- SW
}

# ------------------------- #
# gBGC strength estimation
# ------------------------- #

#'@title gBGC_estimations
#' 
#'@description Enables to make estimate the strength of the gBGC thanks to another R script (from Glémin, 2022).
#'The computations will be done on the data frames produced by the previous function.
#'
#'@param data_mutation1 (character). The first mutation data frame (e.g. : WS)
#'@param data_mutation2 (character). The second mutation data frame (e.g. : SW)
#'
#'@return gBGC.
#'  gBGC strength estimation value.
#' 
#'@example gBGC_estimations("WS", "SW")


gBGC_estimations <- function(data_mutation1, data_mutation2){
  # Use the 'gBGC_estimation_by_least_squares.R' script (Glémin, 2022)
  source("gBGC_estimation_by_least_squares.R")
  
  # Check if the two data have the same length
  cat(nrow(data_mutation1))
  cat(nrow(data_mutation2))
  
  # The computation is the ratio. We needs to remove the 0 values
  rm_zeros = which(data_mutation1$DAF == 0 | data_mutation2$DAF == 0)
  
  # gBGC strength computations
  if(length(rm_zeros) > 0){
    gBGC = least_square(data_mutation1$DAF[- rm_zeros],data_mutation2$DAF[- rm_zeros])
  }else{
    gBGC = least_square(data_mutation1$DAF,data_mutation2$DAF)
  }
  
  cat(paste0("gBGC estimation chromosome ", chromosome_num, " : \n"))
  return(gBGC)
}


# --------------------- #
# Functions execution
# --------------------- #

snp_polarization(species, chromosome_num)
sfs_computations(species, chromosome_num, gff_file, genomic_level)
gBGC_estimations(WS, SW)
