### Test GenomicRanges

# Set working directory
path <- "/home/redelage/internship/data"
setwd(path)

## Packages Intsallation

if(!require("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("GenomicRanges")
install.packages("ape") # For read the GFF files
install.packages("vcfR") # For read the VCF files

# Load the packages
library(GenomicRanges)
library(ape)
library(vcfR)

### Read the gff files ###
# Complete genome GFF file (from NCBI database)
arabidopsis_thaliana_gff <- read.gff(paste0(path,"/GCF_000001735.4_Arabidopsis_thaliana_genomic.gff"))
arabidopsis_lyrata_gff <- read.gff(paste0(path,"/GCF_000004255.2_Arabidopsis_lyrata_genomic.gff"))
brassica_rapa_gff <- read.gff(paste0(path,"/GCF_000309985.2_Brassica_rapa_genomic.gff"))

# Filtered GFF files with only genes features
filtered_at_gff <- read.gff(paste0(path,"/GCF_000001735.4_Arabidopsis_thaliana_genomic_filtered.gff"))
filtered_al_gff <- read.gff(paste0(path, "/GCF_000004255.2_Arabidopsis_lyrata_genomic_filtered.gff"))
filtered_br_gff <- read.gff(paste0(path, "/GCF_000309985.2_Brassica_rapa_genomic_filtered.gff"))

### GRanges objects creations ###


## Arabidopsis thaliana
# All sequences
at_one_base <- which(arabidopsis_thaliana_gff$start >= arabidopsis_thaliana_gff$end)
at_one_base_df <- arabidopsis_thaliana_gff[at_one_base,]
## Remove the NA's of the strnads
arabidopsis_thaliana_gff <- arabidopsis_thaliana_gff[-c(which(is.na(arabidopsis_thaliana_gff$strand))),]
gr_at <- GRanges(seqnames = arabidopsis_thaliana_gff$seqid[-at_one_base],
                 ranges = IRanges(
                   start = arabidopsis_thaliana_gff$start[-at_one_base],
                   end = arabidopsis_thaliana_gff$end[-at_one_base]
                   ),
                 strand = arabidopsis_thaliana_gff$strand[-at_one_base],
                 metadata = data.frame(
                   source = arabidopsis_thaliana_gff$source[-at_one_base],
                   type = arabidopsis_thaliana_gff$type[-at_one_base],
                   attributes = arabidopsis_thaliana_gff$attributes[-at_one_base]
                   )
                 )
summary(gr_at)

# genes only
gr_at_filtered <- GRanges(seqnames = filtered_at_gff$seqid,
                          ranges = IRanges(
                            start = filtered_at_gff$start,
                            end = filtered_at_gff$end,
                          ),
                          strand = filtered_at_gff$strand,
                          metadata = data.frame(
                            source = filtered_at_gff$source,
                            type = filtered_at_gff$type,
                            attributes = filtered_at_gff$attributes
                            )
                          )


## Arabidopsis lyrata
# All sequences
arabidopsis_lyrata_gff <- arabidopsis_lyrata_gff[-c(which(is.na(arabidopsis_lyrata_gff$strand))),]
gr_al <- GRanges(seqnames = arabidopsis_lyrata_gff$seqid,
                 ranges = IRanges(
                   start = arabidopsis_lyrata_gff$start,
                   end = arabidopsis_lyrata_gff$end),
                 strand = arabidopsis_lyrata_gff$strand,
                 metadata = data.frame(
                   source = arabidopsis_lyrata_gff$source,
                   type = arabidopsis_lyrata_gff$type,
                   attributes = arabidopsis_lyrata_gff$attributes)
                 )


# genes only
gr_al_filtered <- GRanges(seqnames = filtered_al_gff$seqid,
                          ranges = IRanges(
                            start = filtered_al_gff$start,
                            end = filtered_al_gff$end,
                          ),
                          strand = filtered_al_gff$strand,
                          metadata = data.frame(
                            source = filtered_al_gff$source,
                            type = filtered_al_gff$type,
                            attributes = filtered_al_gff$attributes
                          )
                        )


### Brassica rapa
# All sequences
gr_br <- GRanges(seqnames = brassica_rapa_gff$seqid,
                 ranges = IRanges(start = brassica_rapa_gff$start, end = brassica_rapa_gff$end),
                 strand = brassica_rapa_gff$strand,
                 metadata = data.frame(
                   source = brassica_rapa_gff$source,
                   type = brassica_rapa_gff$type,
                   attributes = brassica_rapa_gff$attributes)
                 )

# Genes only
  gr_br_filtered <- GRanges(seqnames = filtered_br_gff$seqid,
                            ranges = IRanges(
                              start = filtered_br_gff$start,
                              end = filtered_br_gff$end,
                            ),
                            strand = filtered_br_gff$strand,
                            metadata = data.frame(
                              source = filtered_br_gff$source,
                              type = filtered_br_gff$type,
                              attributes = filtered_br_gff$attributes)
                            )

  
## Get the transcripts sequences of the 3 species
gr_at_transcripts <- gr_at[gr_at$metadata.type != "gene"]
gr_al_transcripts <- gr_al[gr_al$metadata.type != "gene"]
gr_br_transcripts <- gr_br[gr_br$metadata.type != "gene"]


# Store the metadatas
at_metadatas <- mcols(gr_at)
al_metadatas <- mcols(gr_al)
br_metadatas <- mcols(gr_br)

# Group the 3 gr objects of the 3 species
# all_gr <- GRangesList("Arabidopsis_thaliana" = gr_at,
#                       "Arabidopsis_lyrata" = gr_al,
#                       "Brassica_rapa" = gr_br)

# all_gr_filtered <- GRangesList("Arabidopsis_thaliana_filtered" = gr_at_filtered,
#                                "Arabidopsis_lyrata_filtered" = gr_al_filtered,
#                                "Brassica_rapa_filtered" = gr_br_filtered)



## Test function findOverlaps
overlaps_at <- findOverlaps(gr_at_filtered, gr_at_transcripts, type = "within")
overlaps_al <- findOverlaps(gr_al_filtered, gr_al_transcripts, type = "within")
overlaps_br <- findOverlaps(gr_br_filtered, gr_br_transcripts, type = "within")

# Find the primary transcripts (the largest overlapping sequences)
# primary_transcripts_at <- unique(gr_at_transcripts[queryHits(overlaps_at),])
# primary_transcripts_al <- unique(gr_al_transcripts[queryHits(overlaps_al),])
# primary_transcripts_br <- unique(gr_br_transcripts[queryHits(overlaps_br),])


primary_transcripts_at <- findOverlaps(gr_at, gr_at_filtered)
primary_transcripts_al <- findOverlaps(gr_al, gr_al_filtered)
primary_transcripts_br <- findOverlaps(gr_br, gr_br_filtered)

# ## Save the primary transcripts into GFF file
# 
# 
# 
# #Read the VCF file
# arabidopsis_thaliana_vcf <- read.vcfR("Arabidopsis_thaliana_1001genomes.pop.vcf",
#                                       verbose = FALSE)
# 
# # Load the reference genome (fasta format ; ape package)
# reference <- read.FASTA("GCF_000001735.4_arabidopsis_thaliana_genomic.fna", type = "DNA")
# 
# 
