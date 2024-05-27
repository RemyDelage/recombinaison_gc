library(utils)

directory = "~/internship/data/test/"
pattern = "*_filtered_snp_counts.tsv"

merge_data <- function(directory, pattern = "*_filtered_snp_counts.tsv", species){
  # Setting working directory
  setwd(directory)
  
  # Get the list of the files into the directory
  files <- list.files(pattern = pattern)
  
  # Read the first data frame
  results <- read.table(files[1], header = TRUE, sep = '\t', stringsAsFactors = TRUE)
  
  # Create the progress bar
  progress_bar <- txtProgressBar(min = 0, max = length(files), style = 3)
  
  # Loop for merge datas
  for (i in 2:length(files)){
    tmp <- read.table(files[i], header = TRUE, sep = '\t', stringsAsFactors = TRUE)
    results <- rbind(results, tmp)
    setTxtProgressBar(progress_bar, i)
  }
  
  # Save the table
  write.table(results, paste0(species,"_bed_counts.tsv"), sep = '\t', quote = FALSE, 
              row.names = FALSE, col.names = TRUE)
    
  # Close the progress bar
  close(progress_bar)
  
  return(results)
  
}

merge_data(directory, pattern, species = "Arabidopsis")

