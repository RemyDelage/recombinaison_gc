### Test GenomicRanges

# Set working directory
path <- "/home/redelage/internship/data"
setwd(path)

## Packages Intsallation
packages_installation <- function(){
  if(!require("BiocManager"))
    install.packages("BiocManager")
    BiocManager::install("GenomicRanges")
  
  if(!require("ape"))
    install.packages("ape")
  
  if(!require("vcfR"))
    install.packages("vcfR")
}

# Load the packages
library(GenomicRanges)
library(ape)
library(vcfR)

### Setting variables
# Accession numbers
accession1 <- "GCF_000001735.4"
accession2 <- "GCF_000004255.2"
accession3 <- "GCF_000309985.2"

# Species names
species1 <- "Arabidopsis_thaliana"
species2 <- "Arabidopsis_lyrata"
species3 <- "Brassica_rapa"

### Read the gff files ###
# Complete genome GFF file (from NCBI database)

read_gff_files <- function(species1, species2, species3, 
                          accession1, accession2, accession3, 
                          gff_type, variable_suffix ,path){
  gff_files_list <- list()
  gff_files_list[[species1]] <- read.gff(file.path(path, paste0(accession1, "_", species1, "_", gff_type)))
  gff_files_list[[species2]] <- read.gff(file.path(path, paste0(accession2, "_", species2, "_", gff_type)))
  gff_files_list[[species3]] <- read.gff(file.path(path, paste0(accession3, "_", species3, "_", gff_type)))
  assign(paste0(species1, variable_suffix), gff_files_list[[species1]], envir = .GlobalEnv)
  assign(paste0(species2, variable_suffix), gff_files_list[[species2]], envir = .GlobalEnv)
  assign(paste0(species3, variable_suffix), gff_files_list[[species3]], envir = .GlobalEnv)
}

read_gff_files(species1, species2, species3, accession1, accession2, accession3, "genomic.gff", "_gff", path)


# Filtered GFF files with only genes features

read_gff_files(species1, species2, species3, accession1, accession2, accession3, "genomic_filtered.gff", "_genes_gff", path)

### GRanges objects creations ###


## Arabidopsis thaliana
# All sequences
one_base <- which(Arabidopsis_thaliana_gff$start >= Arabidopsis_thaliana_gff$end)
Arabidopsis_thaliana_gff <- Arabidopsis_thaliana_gff[-one_base,]
## Remove the NA's of the strnads
Arabidopsis_thaliana_gff <- Arabidopsis_thaliana_gff[-c(which(is.na(Arabidopsis_thaliana_gff$strand))),]

gr_at <- GRanges(seqnames = Arabidopsis_thaliana_gff$seqid,
                 ranges = IRanges(
                   start = Arabidopsis_thaliana_gff$start,
                   end = Arabidopsis_thaliana_gff$end),
                 strand = Arabidopsis_thaliana_gff$strand,
                 metadata = data.frame(
                   source = Arabidopsis_thaliana_gff$source,
                   type = Arabidopsis_thaliana_gff$type,
                   attributes = Arabidopsis_thaliana_gff$attributes
                   )
                 )

# genes only
gr_at_filtered <- GRanges(seqnames = Arabidopsis_thaliana_genes_gff$seqid,
                          ranges = IRanges(
                            start = Arabidopsis_thaliana_genes_gff$start,
                            end = Arabidopsis_thaliana_genes_gff$end,
                          ),
                          strand = Arabidopsis_thaliana_genes_gff$strand,
                          metadata = data.frame(
                            source = Arabidopsis_thaliana_genes_gff$source,
                            type = Arabidopsis_thaliana_genes_gff$type,
                            attributes = Arabidopsis_thaliana_genes_gff$attributes
                            )
                          )


## Arabidopsis lyrata
# Remove the NAs for the strand data
Arabidopsis_lyrata_gff <- Arabidopsis_lyrata_gff[-c(which(is.na(Arabidopsis_lyrata_gff$strand))),]

# All sequences
gr_al <- GRanges(seqnames = Arabidopsis_lyrata_gff$seqid,
                 ranges = IRanges(
                   start = Arabidopsis_lyrata_gff$start,
                   end = Arabidopsis_lyrata_gff$end),
                 strand = Arabidopsis_lyrata_gff$strand,
                 metadata = data.frame(
                   source = Arabidopsis_lyrata_gff$source,
                   type = Arabidopsis_lyrata_gff$type,
                   attributes = Arabidopsis_lyrata_gff$attributes)
                 )


# genes only
gr_al_filtered <- GRanges(seqnames = Arabidopsis_lyrata_genes_gff$seqid,
                          ranges = IRanges(
                            start = Arabidopsis_lyrata_genes_gff$start,
                            end = Arabidopsis_lyrata_genes_gff$end,
                          ),
                          strand = Arabidopsis_lyrata_genes_gff$strand,
                          metadata = data.frame(
                            source = Arabidopsis_lyrata_genes_gff$source,
                            type = Arabidopsis_lyrata_genes_gff$type,
                            attributes = Arabidopsis_lyrata_genes_gff$attributes
                          )
                        )


### Brassica rapa
# All sequences
gr_br <- GRanges(seqnames = Brassica_rapa_gff$seqid,
                 ranges = IRanges(start = Brassica_rapa_gff$start, 
                                  end = Brassica_rapa_gff$end),
                 strand = Brassica_rapa_gff$strand,
                 metadata = data.frame(
                   source = Brassica_rapa_gff$source,
                   type = Brassica_rapa_gff$type,
                   attributes = Brassica_rapa_gff$attributes)
                 )

# Genes only
gr_br_filtered <- GRanges(seqnames = Brassica_rapa_genes_gff$seqid,
                            ranges = IRanges(
                              start = Brassica_rapa_genes_gff$start,
                              end = Brassica_rapa_genes_gff$end,
                            ),
                            strand = Brassica_rapa_genes_gff$strand,
                            metadata = data.frame(
                              source = Brassica_rapa_genes_gff$source,
                              type = Brassica_rapa_genes_gff$type,
                              attributes = Brassica_rapa_genes_gff$attributes)
                            )

  
# ## Get the transcripts sequences of the 3 species
# gr_at_transcripts <- gr_at[gr_at$metadata.type != "gene"]
# gr_al_transcripts <- gr_al[gr_al$metadata.type != "gene"]
# gr_br_transcripts <- gr_br[gr_br$metadata.type != "gene"]


# Store the metadatas
# at_metadatas <- mcols(gr_at)
# al_metadatas <- mcols(gr_al)
# br_metadatas <- mcols(gr_br)


## Test function findOverlaps
overlaps_at <- findOverlaps(gr_at, gr_at_filtered, type = "within")
overlaps_al <- findOverlaps(gr_al, gr_al_filtered, type = "within")
overlaps_br <- findOverlaps(gr_br, gr_br_filtered, type = "within")

## Get the metadata associated with overlaps

get_metadata <- function(overlaps, gr_data, gr_filtered_data) {
  query_gr <- queryHits(overlaps)
  subject_gr <- subjectHits(overlaps)
  
  query_metadata <- mcols(gr_data)[query_gr,]
  subject_metadata <- mcols(gr_filtered_data)[subject_gr,]
    primary_transcript <- data.frame(Query_Metadata = query_metadata,
                                     Subject_Metadata = subject_metadata)
  return(primary_transcript)
}

at_overlapping_data <- get_metadata(overlaps_at, gr_at, gr_at_filtered)
colnames(at_overlapping_data) = c("Original_source", "Original_type", "Original_attributes",
                                   "Overlap_source", "Overlap_type", "Overlap_attributes")

al_overlapping_data <- get_metadata(overlaps_al, gr_al, gr_al_filtered)
colnames(al_overlapping_data) = c("Original_source", "Original_type", "Original_attributes",
                                  "Overlap_source", "Overlap_type", "Overlap_attributes")

br_overlapping_data <- get_metadata(overlaps_br, gr_br, gr_br_filtered)
colnames(br_overlapping_data) = c("Original_source", "Original_type", "Original_attributes",
                                  "Overlap_source", "Overlap_type", "Overlap_attributes")


# Keep only the overlapping genes sequences
at_primary_transcripts <- subset(at_overlapping_data, Original_attributes == Overlap_attributes)
al_primary_transcripts <- subset(al_overlapping_data, Original_attributes == Overlap_attributes)
br_primary_transcripts <- subset(br_overlapping_data, Original_attributes == Overlap_attributes)

## Redefine the columns names for megre these columns :

    colnames(at_primary_transcripts) = c("source", "type", "attributes", 
                                       "source", "type", "attributes")
  colnames(al_primary_transcripts) = c("source", "type", "attributes", 
                                       "source", "type", "attributes")
  colnames(br_primary_transcripts) = c("source", "type", "attributes", 
                                       "source", "type", "attributes")

# Merge the same columns of the dataframe
at_primary_transcripts <- at_primary_transcripts[, -c(4, 5, 6)]
al_primary_transcripts <- al_primary_transcripts[, -c(4, 5, 6)]
br_primary_transcripts <- br_primary_transcripts[, -c(4, 5, 6)]

## Save the primary transcripts into dataframes
write.table(at_primary_transcripts, paste0(accession1, "_", species1, "_primary_transcripts.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(al_primary_transcripts, paste0(accession2, "_", species2, "_primary_transcripts.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(br_primary_transcripts, paste0(accession3, "_", species3, "_primary_transcripts.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

