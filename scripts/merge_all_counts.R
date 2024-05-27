library(dplyr)
setwd("~/internship/data/test/")
arabido_counts <- read.table('Arabidopsis_bed_counts.tsv', header = TRUE, sep = '\t', stringsAsFactors = TRUE)
data_info <- read.table('Arabidopsis_thaliana_data_info.tsv', header = TRUE, sep = '\t', stringsAsFactors = TRUE)
data_info <- data_info[-c(434802, 80332, 400241, 150778, 375522),]
data_info <- data_info %>%
  mutate(TRANSCRIPT_ID = ifelse(TRANSCRIPT_ID == "", GENE_ID, TRANSCRIPT_ID))
  

data_info_merged <- data_info %>%
  left_join(arabido_counts, by = c("CHROMOSOME_ID" = "CHROMOSOME_ID", "SNP_POSITION" = "SNP_POSITION"))

data_info_merged[is.na(data_info_merged)] <- 0

write.table(data_info_merged, 'Arabidopsis_thaliana_data_info.tsv', sep = '\t', col.names = TRUE,
            row.names = FALSE, quote = FALSE)
