setwd("~/internship/Results_Mar13_arabidopsis/")
orthogroups <- read.table("Orthogroups.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(orthogroups) <- c("Orthogroup",
                           "GCA_000001735_2_Arabidopsis_thaliana", 
                           "GCF_000004255_2_v_1_0_Arabidopsis_lyrata_subsp_lyrata",
                           "GCF_000309985_2_Brassica_rapa")

#install.packages("BiocManager")
#BiocManager::install("Biostrings")
library(Biostings)

orthogroup_file <- "OG0000016.fa"
orthogroup_row <- orthogroups[orthogroups$Orthogroup == "OG0000016", ]
col_names <- names(orthogroups)[-1]

fasta_data <- readDNAStringSet(orthogroup_file)

modify_sequence_id <- function(fasta, col_names, orthogroup_row){
  for (i in seq_along(fasta)){
    print(i)
    sequence_name <- names(fasta)[i]
    gene_id <- sub(">", "", sequence_name)
    
    for (col_name in col_names){
      if (gene_id %in% unlist(strsplit(orthogroup_row[[col_name]], ", "))) {
        new_seq_name <- paste0(col_name, "_", gene_id)
      names(fasta)[i] <- new_seq_name
      
    break
      }
    }
  }
  return(fasta)
}

new_fasta_file <- modify_sequence_id(fasta_data, col_names, orthogroup_row)

writeXStringSet(new_fasta_file, "OG0000016.faa")
  