## Setting working directory
setwd("/home/redelage/internship/data/arabidopsis_thaliana/")

library(dplyr)


getting_species_datas <- function (data_species, ids_file, output_file) {
  at_data <- read.table(data_species, header = TRUE, sep = '\t', stringsAsFactors = TRUE)
  # at_data <- read.table("GCF_000001735.4_Arabidopsis_thaliana_snp_informations.tsv", header = TRUE, sep = '\t', stringsAsFactors = TRUE)
  at_data$TRANSCRIPT_ID <- as.factor(at_data$TRANSCRIPT_ID)



  ids <- read.csv(ids_file, header = TRUE, sep = '\t', stringsAsFactors = TRUE)
  #ids <- read.csv("GCF_000001735.4_Arabidopsis_thaliana_correspondance_genes_transcripts.tsv", header = TRUE, sep = '\t', stringsAsFactors = TRUE)

  matching_indices <- match(at_data$GENE_ID, ids$GENE_ID)
  at_data$TRANSCRIPT_ID <- ids$TRANSCRIPT_ID[matching_indices]

  write.table(at_data, output_file, sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)
  #write.table(at_data, "GCF_000001735.4_Arabidopsis_thaliana_informations.tsv", sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)
}

getting_species_datas("GCF_000001735.4_Arabidopsis_thaliana_snp_informations.tsv",
                      "GCF_000001735.4_Arabidopsis_thaliana_correspondance_genes_transcripts.tsv",
                      "test5.tsv")
