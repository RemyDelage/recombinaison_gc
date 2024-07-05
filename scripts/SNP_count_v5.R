# On va reprendre de zero la polarisation
# Et poser les bases (simple) de ce qu'on fait
# Ca ne veut pas dire qu'il faut jeter tout ce qu'on a fait,
# mais on a decide de simplifier de beaucoup le pipeline en alignant les genomes entiers
# Donc on peut simplifier la partie analyse des SNPs

# Je n'analyse que le chromosome 1
# chromosome 1 = NC_003070.9


# Preparation du VCF
# vcftools --gzvcf Arabidopsis_thaliana_1001genomes.pop_sfs.no_indels.recode.vcf.gz --recode --out Arabidopsis_thaliana_1001genomes.pop_sfs.filtered --max-missing 0.9 --min-alleles 2 --max-alleles 2
# bgzip Arabidopsis_thaliana_1001genomes.pop_sfs.filtered.recode.vcf
# vcftools --counts --gzvcf Arabidopsis_thaliana_1001genomes.pop_sfs.filtered.recode.vcf.gz --chr 1 --out Arabidopsis_thaliana_1001genomes.1

setwd("/home/genouest/cnrs_umr6553/rdelage/test")

#library(vcfR)
library(dplyr)

# Tout d'abord, on a un fichier VCF qui contient l'information des SNPs dans une population d'A. thaliana
# vcf = read.vcfR("/home/genouest/cnrs_umr6553/rdelage/data/Arabidopsis_thaliana_1001genomes.pop_sfs.no_indels.recode.vcf.gz", verbose = FALSE )
# vcf
# 
# queryMETA(vcf)
# 
# head(getFIX(vcf))
# 
# vcf@gt[1:6, 1:4]

## On ne peut pas lire le VCF (fichier trop lourd, non supporté par R)

# Puisque desormais on aligne le genome entier, on va garder tous les SNPs
# Avant on filtrait uniquement les exons car pour OrthoFinder/MACSE on n'alignait que les gènes

# ---------------------------- #
# VCF
# ---------------------------- #

# On compte les allèles de chaque SNP
# Pour cela on a fait 'vcftools --counts --gzvcf Arabidopsis_thaliana_1001genomes.pop.vcf.gz --chr 1 --out Arabidopsis_thaliana_1001genomes.1'

# 'vcftools --counts --gzvcf Arabidopsis_thaliana_1001genomes.pop_sfs.filtered.recode.vcf.gz --chr 1 --out Arabidopsis_thaliana_1001genomes.1'

# On fait aussi 'vcftools --freq --gzvcf Arabidopsis_thaliana_1001genomes.pop.vcf.gz --chr 1 --out Arabidopsis_thaliana_1001genomes.1' pour avoir plus d'infos
# Toute l'info (freq allelique, counts, alleles de ref et alt) sont dans 'Arabidopsis_thaliana_1001genomes.1.frq.count'


# snp_counts = read.table("/home/genouest/cnrs_umr6553/rdelage/results/04_wgabed/04a_Arabidopsis_thaliana/Arabidopsis_thaliana_1001genomes.1.frq.count", header = F, sep = "\t", skip = 1)
# 
# head(snp_counts)
# #head(getFIX(vcf))
# # On voit que le premier allele est le REF et le deuxieme est le ALT
# # On nomme les colonnes
# colnames(snp_counts) = c("CHROM","POS","N_ALLELES","N_CHR","{ALLELE:COUNT}_REF","{ALLELE:COUNT}_ALT")
# head(snp_counts)
# 
# # Mais ce n'est pas encore pratique à utiliser
# snp_counts$REF = gsub(":[0-9]+", "", snp_counts$`{ALLELE:COUNT}_REF`)
# snp_counts$ALT = gsub(":[0-9]+", "", snp_counts$`{ALLELE:COUNT}_ALT`)
# 
# snp_counts$count_REF = gsub("[ATCG]:", "", snp_counts$`{ALLELE:COUNT}_REF`)
# snp_counts$count_ALT = gsub("[ATCG]:", "", snp_counts$`{ALLELE:COUNT}_ALT`)
# 
# head(snp_counts)
# 
# # Maintenant on peux faire la conversion du comptage en quatre colonnes
# snp_counts$count_A = as.numeric(ifelse(snp_counts$REF == "A", snp_counts$count_REF, ifelse(snp_counts$ALT == "A", snp_counts$count_ALT, 0)))
# snp_counts$count_T = as.numeric(ifelse(snp_counts$REF == "T", snp_counts$count_REF, ifelse(snp_counts$ALT == "T", snp_counts$count_ALT, 0)))
# snp_counts$count_C = as.numeric(ifelse(snp_counts$REF == "C", snp_counts$count_REF, ifelse(snp_counts$ALT == "C", snp_counts$count_ALT, 0)))
# snp_counts$count_G = as.numeric(ifelse(snp_counts$REF == "G", snp_counts$count_REF, ifelse(snp_counts$ALT == "G", snp_counts$count_ALT, 0)))
# 
# head(snp_counts)
# 
# # On compte combien de genomes ne sont pas en missing data
# # On a 184 individus donc 368 genomes au max
# snp_counts$n_genomes = rowSums(snp_counts[,c("count_A", "count_T", "count_C", "count_G")])
# head(snp_counts)
# 
# 
# # Comptage des frequences alleliques REF/ALT
# # Si on assume que REF == Ancestral (le travail a deja ete fait)
# snp_counts$mut = unlist(lapply(1:nrow(snp_counts), function(x) {paste(snp_counts$REF[x], snp_counts$ALT[x], sep = "->")}))
# head(snp_counts$mut)
# table(snp_counts$mut)
# 
# # WS
# sum(snp_counts$mut %in% c("A->G", "T->G", "A->C", "T->C"))
# # SW
# sum(snp_counts$mut %in% c("G->A", "G->T", "C->A", "C->T"))
# 
# # On a bien a peu pres autant de WS que SW
# 
# # Recuperation de l'info AA (Ancestral Allele) ?
# # qu'on pourra comparer avec notre polarisation
# # Mais pas d'info AA
# 
# 
# # Benchmark de la methode ?
# 
# 
# 
# # Maintenant on a fini avec la partie SNP (population)
# 
# 
# # ---------------------------- #
# # WGABED - Alignement avec Cactus
# # ---------------------------- #
# # On part de l'alignement des genomes avec Cactus
# 
# # On a aligne des genomes entiers entre eux, genome de ref et deux outgroups,
# # et pour chaque position du genome on a regarde quels sont
# # les nucleotides chez les outgroups
# 
# # On a donc bien UN et UN seul nucleotide pour chaque espece
# # qui est celui du genome FASTA
# # Il peut donc etre different des deux alleles du VCF, meme si ca doit etre rare, si il y a plusieurs mutations au meme site
# 
# # WGABed permet de reformater la sortie de Cactus en data frame au format BED
# # On recupere l'info en sortie de WGAbed
# # Pour info, une ligne = une position dans l'alignement/genome de reference
# 
# # ATTENTION
# # ON prend la sortie de WGAbed non filtree = tous les sites
# # Car desormais on s'interesse a tout le genome, pas seulement les exons
# # 'python maf_to_bed_v2.py -i data/arabidopsis_genomes_alignment.maf.gz -r GCF_000001735_Arabidopsis_thaliana -c NC_003070.9 -s data/species_list | sort -k1,1 -k2,2n | bgzip -c > arabidopsis_NC_003070.9.wga.bed.gz'
# 
# 
# # ATTENTION C'est un .BED, donc 0-based
# # Les VCFs sont 1-based
# 
# # wgabed = read.table(gzfile("data/test.wga.bed.gz"))
# # head(wgabed)
# 
# wgabed = read.table(gzfile("/home/genouest/cnrs_umr6553/rdelage/results/04_wgabed/04a_Arabidopsis_thaliana/arabidopsis_NC_003076.8.wga.bed.gz"))
# colnames(wgabed) = c("seq", "start", "end", "strand", "species", "aligned_chrom", "aligned_pos", "sequence", "strands", "score")
# head(wgabed)
# 
# wgabed$CHROM = "5"
# # Keep only SNPs
# wgabed$width = wgabed$end - wgabed$start
# cat("wgabed before SNP filtering \n")
# nrow(wgabed)
# wgabed = wgabed[which(wgabed$width == 1),]
# cat("wgabed after SNP filtering \n")
# nrow(wgabed)
# 
# 
# # ATTENTION C'est un .BED, donc 0-based
# # Les VCFs sont 1-based
# wgabed$POS = wgabed$start + 1
# head(wgabed)
# 
# # Etape a verifier a la fin : regarder si les alleles du VCF matchent bien avec le genome de ref
# # (on s'attend a avoir beaucoup de REF == ref)
# # Sinon il faudra verifier si le shift 0-based/1-based a bien fonctionné
# 
# 
# 
# # Il faut recuperer l'info suivante :
# # quels sont les nucleotides chez les deux outgroups a la position donnée de chaque SNP dans le vcf
# 
# # On recupere aussi la ref
# wgabed$ref = toupper(unlist(lapply(strsplit(wgabed$sequence, ","), function(x) {x[1]})))
# wgabed$outgroup_1 = toupper(unlist(lapply(strsplit(wgabed$sequence, ","), function(x) {x[2]})))
# wgabed$outgroup_2 = toupper(unlist(lapply(strsplit(wgabed$sequence, ","), function(x) {x[3]})))
# head(wgabed)
# 
# # TODO
# # Enlever les lignes avec ref ou outgroup qui n'est pas une seule base (length > 1)
# # Bugfix. Toutes les colonnes etaient passees a T -> FIXED
# wgabed$ref = unlist(lapply(1:nrow(wgabed), function(x) {ifelse(nchar(wgabed$ref[x]) == 1, wgabed$ref[x], NA)}))
# head(wgabed)
# wgabed$outgroup_1 = unlist(lapply(1:nrow(wgabed), function(x) {ifelse(nchar(wgabed$outgroup_1[x]) == 1, wgabed$outgroup_1[x], NA)}))
# head(wgabed)
# wgabed$outgroup_2 = unlist(lapply(1:nrow(wgabed), function(x) {ifelse(nchar(wgabed$outgroup_2[x]) == 1, wgabed$outgroup_2[x], NA)}))
# head(wgabed)
# 
# # View(wgabed[1:1000,])
# 
# 
# # On merge le wgabed avec le vcf
# # On garde le VCF comme fichier de reference
# head(snp_counts)
# head(wgabed)
# 
# summary(snp_counts)
# summary(wgabed)
# 
# nrow(snp_counts)
# nrow(wgabed)
# 
# snp_counts$CHROM = as.character(snp_counts$CHROM)
# wgabed$CHROM = as.character(wgabed$CHROM)
# 
# # Ici on fait face a une erreur bizarre quand on merge
# # WGAbed a des positions dupliquees (bizarre...)
# table(wgabed$CHROM)
# length(wgabed$POS)
# length(unique(wgabed$POS))
# head(wgabed[duplicated(wgabed$POS),])
# 
# cat("remove duplicates \n")
# wgabed = wgabed[!duplicated(wgabed$POS), ]
# length(wgabed$POS)
# length(unique(wgabed$POS))
# cat("Have to be the same \n")
# 
# # On a enleve toutes les positions dupliquees avant de merger
# head(wgabed)
# snp_counts_merged = left_join(snp_counts, wgabed[,c("CHROM", "POS", "ref", "outgroup_1", "outgroup_2")])
# 
# 
# cat("Needs to be the same number of lines \n")
# nrow(snp_counts)
# nrow(snp_counts_merged)
# 
# 
# head(snp_counts_merged)
# 
# # View(snp_counts_merged[1:100,])
# 
# # A cette etape on peut verifier que REF == ref
# # L'allele de reference du VCF et du FASTA doivent correspondre dans la majorite des cas
# head(snp_counts_merged[,c("REF", "ref")])
# 
# sum(snp_counts_merged$REF == snp_counts_merged$ref, na.rm = T)/nrow(snp_counts_merged)
# #A quoi correspond ce ratio ?
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


# ---------------------------- #
# EST-SFS
# ---------------------------- #
# Maintenant on peut reformater le tableau pour correspondre au format d'input de EST-SFS


# Manuel d'EST-SFS
# "If there are three outgroups, there are 4 space-separated columns. The first column is for the  focal species, and the next three columns are for the outgroups.  Each column is a comma-separated list of the counts of the four bases in the order A, C, G, T.  For example, the first line in the example data file is:
# 20,0,0,0        0,0,0,1 0,0,0,1 0,0,0,1"

# countRef = unlist(lapply(1:nrow(snp_counts_merged), function(x) {paste(snp_counts_merged[x,c("count_A", "count_C", "count_G", "count_T")], collapse = ",")}))
# countOut1 = unlist(lapply(1:nrow(snp_counts_merged), function(x) {paste(snp_counts_merged[x,c("A_outgroup_2", "C_outgroup_1", "G_outgroup_1", "T_outgroup_1")], collapse = ",")}))
# countOut2 = unlist(lapply(1:nrow(snp_counts_merged), function(x) {paste(snp_counts_merged[x,c("A_outgroup_2", "T_outgroup_2", "G_outgroup_2", "T_outgroup_2")], collapse = ",")}))
# 
# 
# df = data.frame(ref_species = countRef,
#                 outgroup_1 = countOut1,
#                 outgroup_2 = countOut2)

df = read.table("/home/genouest/cnrs_umr6553/rdelage/results/06_est-sfs/06b_Oryza_sativa/Oryza_sativa_chrom_7_est-sfs-input.txt")
head(df)
nrow(df)

# PROBLEME. est-sfs n'accepte pas plus de 1,000,000 SNPs en entrée
# Je fais tourner uniquement sur 999,999 SNPs pour tester
# A toi de faire plusieurs batchs de SNPs et de faire tourner en parallele est-sfs,
# puis de merger tous ces batchs
#write.table(df[1:(1000000-1),], "data/est-sfs-input.txt", row.names = F, col.names = F, quote = F, sep = " ")

# On fait tourner EST-SFS
# est-sfs data/config_file.txt data/est-sfs-input.txt data/seedfile.txt data/output-file-sfs.txt data/output-file-pvalues.txt 

# ---------------------------- #
# Recuperer les resultats d'EST-SFS et polariser finalement les alleles
# ---------------------------- #

snp_counts_merged <- read.table("/home/genouest/cnrs_umr6553/rdelage/results/04_wgabed/04b_Oryza_sativa/oryza_snp_counts_merged_chrom7", header=T, sep="\t")
head(snp_counts_merged)

# estsfs = read.table("arabidopsis_NC_003070.9.wga_filtered_output_file_pvalues.txt", header = F, skip = 8, sep = " ")
estsfs = read.table("/home/genouest/cnrs_umr6553/rdelage/results/06_est-sfs/06b_Oryza_sativa/Oryza_sativa_chrom_7_output_file_pvalues.txt", header = F, skip = 8, sep = " ")
estsfs = estsfs[,-8]
colnames(estsfs) = c("IdxSite", "Code", "P_major_ancestral", "P_trees_A", "P_trees_C", "P_trees_G", "P_trees_T")
head(estsfs)
nrow(estsfs)
# nrow(df)
# cat("Have to be the same \n")

# estsfs = estsfs[1:nrow(df),]

# On merge a nouveau les resultats avec snp_counts_merged
snp_counts_merged$P_major_ancestral = NA
snp_counts_merged$P_major_ancestral[estsfs$IdxSite] = estsfs$P_major_ancestral
snp_counts_merged$P_trees_A = NA
snp_counts_merged$P_trees_A[estsfs$IdxSite] = estsfs$P_trees_A
snp_counts_merged$P_trees_C = NA
snp_counts_merged$P_trees_C[estsfs$IdxSite] = estsfs$P_trees_C
snp_counts_merged$P_trees_G = NA
snp_counts_merged$P_trees_G[estsfs$IdxSite] = estsfs$P_trees_G
snp_counts_merged$P_trees_T = NA
snp_counts_merged$P_trees_T[estsfs$IdxSite] = estsfs$P_trees_T

# Identifier l'allele ancestral et le derive
snp_counts_merged[10,]

# Faire le compte
major = function(x) {
    # Return the major allele for a line x in snp_counts_merged
    alphabet = c("A", "C", "G", "T")
    maj = alphabet[which.max(snp_counts_merged[x, c("count_A", "count_C", "count_G", "count_T")])]
    return(maj)
}
cat("Major allele position 10 : \n", major(10))

# minor = function(x) {
#     # Return the minor allele for a line x in snp_counts_merged
#     alphabet = c("A", "C", "G", "T")
#     s = snp_counts_merged[x, c("count_A", "count_C", "count_G", "count_T")]
#     s[s == 0] = NA
#     mino = alphabet[which.min(s)]
#     return(mino)
# }
# minor(10)
minor = function(x) {
  alphabet = c("A", "C", "G", "T")
  s = snp_counts_merged[x, c("count_A", "count_C", "count_G", "count_T")]
  if (all(s == 0)) {
    return(NA)
  } else {
    s[s == 0] = NA
    mino = alphabet[which.min(s)]
    return(mino)
  }
}
cat("Minor allele position 10 : \n", minor(10))

snp_counts_merged$Major = unlist(lapply(1:nrow(snp_counts_merged), major))
snp_counts_merged$Minor = unlist(lapply(1:nrow(snp_counts_merged), minor))

# TODO. BugFix. Il y a certaines lignes avec un seul allele, quie st donc a la fois majeur et mineur
# Retirer ces lignes
snp_counts_merged = snp_counts_merged[which(snp_counts_merged$Major != snp_counts_merged$Minor),]

head(snp_counts_merged)

# On prend pour allele ancestral celui avec la plus grande proba sur l'arbre
ancestral = function(x) {
    # Return the ancestral allele for a line x in snp_counts_merged
    alphabet = c("A", "C", "G", "T")
    proba_alleles = snp_counts_merged[x, c("P_trees_A", "P_trees_C", "P_trees_G", "P_trees_T")]
    anc = alphabet[which.max(proba_alleles)]
    anc = ifelse(length(anc) == 1, anc, NA)
    return(anc)
}
cat("Ancestral allele position 10 : ", ancestral(10), "\n")
snp_counts_merged$Ancestral = unlist(lapply(1:nrow(snp_counts_merged), ancestral))

# Consider then that the other allele in the VCF is derived
snp_counts_merged$Derived = ifelse(is.na(snp_counts_merged$Ancestral), NA, ifelse(snp_counts_merged$Ancestral == snp_counts_merged$Major, snp_counts_merged$Minor, snp_counts_merged$Major))


# Bug to fix, dans le cas ou l'ancestral n'est ni le majeur, ni le mineur -> cas non resolu, mettre valeur NA
# ...
snp_counts_merged[which(snp_counts_merged$Ancestral != snp_counts_merged$Minor & snp_counts_merged$Ancestral != snp_counts_merged$Major), c("Ancestral", "Derived")] = NA

head(snp_counts_merged)

# Ici tu peux placer ton code pour definir les differentes categories de mutation WS, SW, W, S
# ...
snp_counts_merged$Mutation = NA
snp_counts_merged$Mutation[which(snp_counts_merged$Ancestral %in% c("A", "T") & snp_counts_merged$Derived %in% c("C", "G"))] = "WS"
snp_counts_merged$Mutation[which(snp_counts_merged$Ancestral %in% c("C", "G") & snp_counts_merged$Derived %in% c("A", "T"))] = "SW"
snp_counts_merged$Mutation[which(snp_counts_merged$Ancestral %in% c("A", "T") & snp_counts_merged$Derived %in% c("A", "T"))] = "W"
snp_counts_merged$Mutation[which(snp_counts_merged$Ancestral %in% c("C", "G") & snp_counts_merged$Derived %in% c("C", "G"))] = "S"


write.table(snp_counts_merged, "snp_counts_merged_chrom7.txt", row.names = F, col.names = T,quote = F, sep = "\t")


# A la fin on a un tableau avec une ligne par SNP dans le VCF --> Problème : on ne peut pas ouvrir le VCF avec R
# Pour chaque SNP on connait l'allele ancestral et l'allele derive
# vcf_fix = as.data.frame(getFIX(vcf))
# sum(vcf_fix$CHROM == 1)
nrow(snp_counts_merged)

# On va pouvoir construire le SFS a partir de ce tableau
# Pour le construire, tu peux utiliser tous les SNPs
# Mais ce serait bien de filtrer par qualité (regarder les infos qualité dans le VCF)
# queryMETA(vcf)
# ou par qualite de polarisation (fixer des seuils sur les probas P-tree, par exemple > 0.8 pour l'allele ancestral)
# enlever tous les sites qui ont des NA/? dans les outgroups



# Et puis tu peux decouper en fonction des categories biologiques pour lesquelles on veut estimer la gBGC
# Tu peux faire l'intersection entre snp_counts_merged (passé en GenomicRanges)
# et un GenomicRanges qui definit, par exemple:
# - les exons
# - les introns
# - uniquement l'exon 1
# ....
# Tu peux aussi compter les SNPs WS, SW, W et S dans chacune des categories (countOverlap)
# et faire des gradients de WS vs. SW, etc


# Tu peux aussi reflechir a des tests a mettre en place a chaque etape pour verifier que la sortie correspond bien a ce qu'on attend


## Pooler ici les résultats des 5 chromosomes (rbind ou merge)


# ---------------------------- #
# Construire le SFS pour la gBGC
# ---------------------------- #

# A toi de jouer...


# snp_counts_merged <- read.table('Arabidopsis_thaliana_NC_003070.9.tsv', 
#                                 header = TRUE, sep = '\t')

# snp_counts_merged <- read.table('data/snp_counts_merged.txt', 
#                                 header = TRUE, sep = '\t')


## Filter datas according to the probabilty
head(snp_counts_merged)
nrow(snp_counts_merged)
summary(snp_counts_merged)

table(snp_counts_merged$Mutation)


# Je fais ici le sous-jeu de donnees genique
# et je refais tourner le code avec pour ne pas dupliquer
# Mais toi tu dois refactoriser pour faire les deux en meme temps

# On marque les SNPs qui sont geniques
# BiocManager::install("GenomicRanges")

library(GenomicRanges)
#gff = readRDS("data/gff_rho_Arabidopsis_thaliana_1001genomes.rds")
gff = readRDS("/home/genouest/cnrs_umr6553/rdelage/data/gff_rho_Oryza_sativa_Wang2018.rds")
head(gff)

gff_ranges = makeGRangesFromDataFrame(gff, keep.extra.columns=T)
gff_ranges
gene_ranges = gff_ranges[which(gff_ranges$feature == "gene")]
gene_ranges

# # Find overlapping 
# colnames(snp_counts_merged)
# snps = snp_counts_merged
# snps$start = snps$POS
# snps$end = snps$POS + 1
# snp_ranges = makeGRangesFromDataFrame(snps)
# snp_ranges
# 
# 
# hits = findOverlaps(snp_ranges, gene_ranges)
# hits
# 
# nrow(snp_counts_merged)
# snp_counts_merged = snp_counts_merged[queryHits(hits),]
# nrow(snp_counts_merged)
# length(queryHits(hits))
# 
# # Et seulement exon 1 ?
# exon1_ranges = gff_ranges[which(gff_ranges$feature == "exon" & gff_ranges$rank == 1 & gff_ranges$nb_exons < 15)]
# exon1_ranges
# 
# # Find overlapping
# hits = findOverlaps(snp_ranges, exon1_ranges)
# hits
# 
# nrow(snp_counts_merged)
# snp_counts_merged = snp_counts_merged[queryHits(hits),]
# nrow(snp_counts_merged)
# length(queryHits(hits))


# snp_counts_merged_filtered_p80 <- subset(snp_counts_merged, P_MAJOR_ANCESTRAL >= 0.80)

# On ne va pas garder ce filtre, qui enleve trop de SNPs
# De plus ca biaise en enlevant aussi les SNPs dont le major n'est pas l'ancestral
# Si tu veux filtrer il veut mieux filtrer sur les colonnes P-trees
# en gardant uniquement si le max de P-trees{A,C,G,T} > 0.8 par exemple
snp_counts_merged$filter = unlist(lapply(1:nrow(snp_counts_merged), function(x) {ifelse(max(snp_counts_merged[x,c("P_trees_A", "P_trees_C", "P_trees_G", "P_trees_T")]) > 0.7, "PASS", "REMOVE")}))

table(snp_counts_merged$filter)

snp_counts_merged_filtered_p80 <- snp_counts_merged[which(snp_counts_merged$filter == "PASS"),]
nrow(snp_counts_merged_filtered_p80)

table(snp_counts_merged_filtered_p80$Derived)

# snp_counts_merged_filtered_p80 = snp_counts_merged_filtered_p80[-which(is.na(snp_counts_merged_filtered_p80$Derived)),]
# snp_counts_merged_filtered_p80 = snp_counts_merged_filtered_p80[-which(is.na(snp_counts_merged_filtered_p80$Ancestral)),]
nrow(snp_counts_merged_filtered_p80)
table(snp_counts_merged_filtered_p80$Derived)


### DAF Computations --> On filtered data
# Nb derivated alleles for each row
snp_counts_merged_filtered_p80 <- snp_counts_merged[which(snp_counts_merged$filter == "PASS"),]
nrow(snp_counts_merged_filtered_p80)

table(snp_counts_merged_filtered_p80$Derived)
snp_counts_merged_filtered_p80$Nb_Derived <- ifelse(snp_counts_merged_filtered_p80$Derived == "A",
                                                    snp_counts_merged_filtered_p80$count_A,
                                                    ifelse(snp_counts_merged_filtered_p80$Derived == "C",
                                                           snp_counts_merged_filtered_p80$count_C,
                                                           ifelse(snp_counts_merged_filtered_p80$Derived == "G",
                                                                  snp_counts_merged_filtered_p80$count_G,
                                                                  snp_counts_merged_filtered_p80$count_T)))
# Nb total alleles for each row

snp_counts_merged_filtered_p80$Nb_Total_Alleles <- snp_counts_merged_filtered_p80$N_CHR

# DAF for each line

# ATTENTION tu as toujours la meme erreur dans tes boucles
# Il faut mettre le i en indice
# for (i in 1:nrow(snp_counts_merged_filtered_p80)){
#   snp_counts_merged_filtered_p80$DAF <- snp_counts_merged_filtered_p80$Nb_Derived / snp_counts_merged_filtered_p80$Nb_Total_Alleles
# }

# # Boucle corrigee
# for (i in 1:nrow(snp_counts_merged_filtered_p80)){
#   snp_counts_merged_filtered_p80$DAF[i] <- snp_counts_merged_filtered_p80$Nb_Derived[i] / snp_counts_merged_filtered_p80$Nb_Total_Alleles[i]
# }


write.table(snp_counts_merged_filtered_p80, "snp_counts_merged_filtered_p80_chrom7.tsv", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
##### Separate the data according to the type of mutations
library(dplyr)
nrow(snp_counts_merged_filtered_p80)
table(snp_counts_merged_filtered_p80$Mutation)



snp_counts_p80_sw <- snp_counts_merged_filtered_p80[which(snp_counts_merged_filtered_p80$Mutation == "SW"),]
snp_counts_p80_ws <- snp_counts_merged_filtered_p80[which(snp_counts_merged_filtered_p80$Mutation == "WS"),]
snp_counts_p80_s <- snp_counts_merged_filtered_p80[which(snp_counts_merged_filtered_p80$Mutation == "S"),]
snp_counts_p80_w <- snp_counts_merged_filtered_p80[which(snp_counts_merged_filtered_p80$Mutation == "W"),]

nrow(snp_counts_p80_sw)
nrow(snp_counts_p80_ws)
nrow(snp_counts_p80_s)
nrow(snp_counts_p80_w)

# We can get the SFS for all individuals with 0% missing data
# But we need to subsample

# I could use the function sent by Sylvain
# But I will rather use the same method as in the tutorial
# It is easier to understand (go back to the tutorial and the paper)

# Pour rappel de la methode :

# Find the minimum sample size across all sites (call this min_ss ), and for each site resample min_ss individuals, such that we create an SFS with data for all sites. This can be extended to take into account linkage, sampling blocks of data. See the paper Bagley et al (2017) Mol Ecol (https://doi.org/10.1111/mec.13972) and scripts available on dryad to do this block sampling without missing data: http://datadryad.org/resource/doi:10.5061/dryad.vh75r/2

# In the script below I exemplify how we could get the SFS by resampling genotypes such that we have min_ss individuals per site with data. The required info is a genotypic matrix.

# number of SNPs
nsnps <- nrow(snp_counts_p80_ws)
# number of individuals
nind <- max(snp_counts_p80_ws$Nb_Total_Alleles, na.rm = T)/2

cat("Number of SNPs : ", nsnps, "\n")
cat("Number of individuals : ", nind, "\n")
# nsnps
# nind


# We need as input
#  A genotype matrix with missing data and
# columns are site
# rows are individuals

# We don't have this matrix (the GT matrix in the vcf) because we already counted the number of alleles

# I create a pseudo-matrix
# with nsnps columns and nind*2 lignes
md_geno = matrix(NA, nrow = nind*2, ncol = nsnps)

# Pour chaque ligne, je cree un vecteur contenant
# n alleles derives (= 1), m alleles ancestraux (= 0) et l NAs
# ou n = snp_counts_p80_ws$Nb_derived
# m = snp_counts_p80_ws$Nb_Total - snp_counts_p80_ws$Nb_derived
# l = nind*2 - snp_counts_p80_ws$Nb_Total

for (i in 1:ncol(md_geno)) {
    n = snp_counts_p80_ws$Nb_Derived[i]
    m = snp_counts_p80_ws$Nb_Total[i] - snp_counts_p80_ws$Nb_Derived[i]
    l = nind*2 - snp_counts_p80_ws$Nb_Total[i]
    md_geno[,i] = c(rep(1, n), rep(0, m), rep(NA, l))
}


# View(md_geno[,1:10])
# On s'en fiche que les 0, 1 et NAs soient ordonnes
# Parce qu'on va refaire un sampling aleatoire juste apres

# find minimum sample size
md_persite <- colSums(is.na(md_geno))
min_ss <- abs(nind-max(md_persite))

cat("min_ss : ", min_ss, "\n")
# GETMAF
# Given an unfolded SFS get the minor allele frequency spectrum (folded SFS)
# INPUT
#   unfolded.sfs : numeric vector of size (nindividuals*2)+2 with the unfolded SFS
# NOTE! Only works for unfolded.sfs with an odd number of entries!!
getMaf <- function(unfolded.sfs) {
  # set the maf sfs the same as obs.sfs
  maf <- as.numeric(unfolded.sfs)
  # define the length of the SFS and half the SFS
  n <- length(maf)
  nhalf <- ceiling(n/2)
  if(n %% 2 == 0) {
    stop("This getMaf function only works for SFS with odd number of entries (for data from diploid individuals.)")
  }
  # all the entries of maf.sfs > (n/2) set to zero, where n is 2*nindividuals
  maf[(nhalf+1):n] <- 0
  # add the entry 1 to n, 2 to n-1, and so on...
  maf[1:(nhalf-1)] <- unfolded.sfs[1:(nhalf-1)]+
                              unfolded.sfs[n:(nhalf+1)]
  maf
}

# Go through each site, sample min_ss genotypes for each site without missing data
no_ms_geno_list <- lapply(1:ncol(md_geno), function(i) {
  each_site <- md_geno[, i]
  population <- sum(!is.na(each_site))
  sample_size <- min(min_ss, population)
  sample(each_site[!is.na(each_site)], size = sample_size, replace = FALSE)
})

# Ensure each sampled vector is of the same length
no_ms_geno_matrix <- do.call(cbind, lapply(no_ms_geno_list, function(x) {
  length_diff <- min_ss - length(x)
  if (length_diff > 0) {
    return(c(x, rep(NA, length_diff)))
  } else {
    return(x)
  }
}))

# Ensure no_ms_geno_matrix is numeric
no_ms_geno_matrix <- apply(no_ms_geno_matrix, 2, as.numeric)

# Check and display missing data
missing_data_indices <- which(is.na(no_ms_geno_matrix), arr.ind = TRUE)
if (length(missing_data_indices) > 0) {
  # Remove columns with missing data
  no_ms_geno_matrix <- no_ms_geno_matrix[, colSums(is.na(no_ms_geno_matrix)) == 0]
}

# Get the SFS for the dataset without missing data
obs.sfs <- table(colSums(no_ms_geno_matrix, na.rm = TRUE))
head(obs.sfs)

jpeg("SFS_CHROMOSOME7_ORYZA.jpg", width = 1200, height = 800)
barplot(obs.sfs, xlab = "derived allele frequency", ylab = "number of sites")
title("SFS Chromosome 7 - genes")
dev.off()



# # Go through each site, sample min_ss genotypes for each site without missing data
# no_ms_geno <- apply(md_geno, 2, function(each_site) {sample(each_site[!is.na(each_site)], size=min_ss, replace = F)})
# # Check that there is no missing data
# if(sum(is.na(no_ms_geno))>0) {stop("Error: there is still missing data after resampling!")}
# 
# # Get the SFS for the dataset without missing data
# obs.sfs <- table(colSums(no_ms_geno))
# # Plot the unfolded SFS
# # On ne s'occupe pas du folded (qu'on utilise sur Minor Allele Frequency QUAND ON NE CONNAIT PAS L'ETAT ANCESTRAL)
# barplot(obs.sfs, xlab="derived allele frequency", ylab="number of sites")




# head(obs.sfs)
# Voila tu as un SFS,
# Avec un probleme majeur lie aux donnees qui est
# l'absence de SNPs en faible frequence

# Tu dois refaire tourner WGAbed et EST-SFS sur le mnouveau vcf
# Arabidopsis_thaliana_1001genomes.pop_sfs_no_indels.vcf.gz

cat("#### Mutation WS #### \n")

# On refait la procedure pour avoir WS et SW
# number of SNPs
nsnps <- nrow(snp_counts_p80_ws)
# number of individuals
nind <- max(snp_counts_p80_ws$Nb_Total_Alleles, na.rm = T)/2

cat("Number of SNPs : ", nsnps, "\n")
cat("Number of individuals : ", nind, "\n")

WS = matrix(NA, nrow = nind*2, ncol = nsnps)
# Pour chaque ligne, je cree un vecteur contenant
# n alleles derives (= 1), m alleles ancestraux (= 0) et l NAs
# ou n = snp_counts_p80_ws$Nb_derived
# m = snp_counts_p80_ws$Nb_Total - snp_counts_p80_ws$Nb_derived
# l = nind*2 - snp_counts_p80_ws$Nb_Total
for (i in 1:ncol(WS)) {
    n = snp_counts_p80_ws$Nb_Derived[i]
    m = snp_counts_p80_ws$Nb_Total[i] - snp_counts_p80_ws$Nb_Derived[i]
    l = nind*2 - snp_counts_p80_ws$Nb_Total[i]
    WS[,i] = c(rep(1, n), rep(0, m), rep(NA, l))
}
md_persite <- colSums(is.na(WS))
#min_ss <- nind - max(md_persite)
min_ss <- abs(nind-max(md_persite))
cat("min_ss : ", min_ss, "\n")

# # Go through each site, sample min_ss genotypes for each site without missing data
# WS_no_missing <- apply(WS, 2, function(each_site) {sample(each_site[!is.na(each_site)], size=min_ss, replace = F)})
# # Check that there is no missing data
# if(sum(is.na(WS_no_missing))>0) {stop("Error: there is still missing data after resampling!")}
# Get the SFS for the dataset without missing data
# WS <- table(colSums(WS_no_missing))

WS_no_missing_list <- lapply(1:ncol(WS), function(i) {
  each_site <- WS[, i]
  population <- sum(!is.na(each_site))
  sample_size <- min(min_ss, population)
  sample(each_site[!is.na(each_site)], size = sample_size, replace = FALSE)
})

# Ensure each sampled vector is of the same length
WS_no_missing_matrix <- do.call(cbind, lapply(WS_no_missing_list, function(x) {
  length_diff <- min_ss - length(x)
  if (length_diff > 0) {
    return(c(x, rep(NA, length_diff)))
  } else {
    return(x)
  }
}))

WS_no_missing_matrix <- apply(WS_no_missing_matrix, 2, as.numeric)


# Check and display missing data
missing_data_indices <- which(is.na(WS_no_missing_matrix), arr.ind = TRUE)
if (length(missing_data_indices) > 0) {
  # Remove columns with missing data
  WS_no_missing_matrix <- WS_no_missing_matrix[, colSums(is.na(WS_no_missing_matrix)) == 0]
}

WS = data.frame(Nb_sites = seq(0, min_ss),
                DAF = NA)
head(WS)
for (i in 1:nrow(WS)) {
  WS$DAF[i] = sum(colSums(WS_no_missing_matrix) == WS$Nb_sites[i])
}
head(WS)
# Plot the unfolded SFS
# On ne s'occupe pas du folded (qu'on utilise sur Minor Allele Frequency QUAND ONE CONNAIT PAS L'ETAT ANCESTRAL)
# barplot(WS$DAF, xlab="derived allele frequency", ylab="number of sites")

jpeg("SFS_WS_CHROMOSOME7_ORYZA.jpg", width = 1200, height = 800)
barplot(WS$DAF, xlab="derived allele frequency", ylab="number of sites")
title("SFS_WS Chromosome 7 - genes")
dev.off()
cat("Plot 'SFS_WS_CHROMOSOME7_ORYZA.jpg' created\n")

cat("#### Mutation SW #### \n")

nsnps <- nrow(snp_counts_p80_sw)
cat("Number of SNPs : ", nsnps, "\n")
SW = matrix(NA, nrow = nind*2, ncol = nsnps)
# Pour chaque ligne, je cree un vecteur contenant
# n alleles derives (= 1), m alleles ancestraux (= 0) et l NAs
# ou n = snp_counts_p80_ws$Nb_derived
# m = snp_counts_p80_ws$Nb_Total - snp_counts_p80_ws$Nb_derived
# l = nind*2 - snp_counts_p80_ws$Nb_Total
for (i in 1:ncol(SW)) {
    n = snp_counts_p80_sw$Nb_Derived[i]
    m = snp_counts_p80_sw$Nb_Total[i] - snp_counts_p80_sw$Nb_Derived[i]
    l = nind*2 - snp_counts_p80_sw$Nb_Total[i]
    SW[,i] = c(rep(1, n), rep(0, m), rep(NA, l))
}
# # md_persite <- colSums(is.na(WS))
# # min_ss <- nind - max(md_persite)
# # Go through each site, sample min_ss genotypes for each site without missing data
# SW_no_missing <- apply(SW, 2, function(each_site) {sample(each_site[!is.na(each_site)], size=min_ss, replace = F)})
# # Check that there is no missing data
# if(sum(is.na(SW_no_missing))>0) {stop("Error: there is still missing data after resampling!")}
# # Get the SFS for the dataset without missing data
# # SW <- table(colSums(SW_no_missing))

no_ms_SW_list <- lapply(1:ncol(SW), function(i) {
  each_site <- SW[, i]
  population <- sum(!is.na(each_site))
  sample_size <- min(min_ss, population)
  sample(each_site[!is.na(each_site)], size = sample_size, replace = FALSE)
})

# Ensure each sampled vector is of the same length
no_ms_SW_matrix <- do.call(cbind, lapply(no_ms_SW_list, function(x) {
  length_diff <- min_ss - length(x)
  if (length_diff > 0) {
    return(c(x, rep(NA, length_diff)))
  } else {
    return(x)
  }
}))

# Ensure no_ms_geno_matrix is numeric
no_ms_SW_matrix <- apply(no_ms_SW_matrix, 2, as.numeric)

# Check and display missing data
missing_data_indices <- which(is.na(no_ms_SW_matrix), arr.ind = TRUE)
if (length(missing_data_indices) > 0) {
  # Remove columns with missing data
  no_ms_SW_matrix <- no_ms_SW_matrix[, colSums(is.na(no_ms_SW_matrix)) == 0]
}

SW = data.frame(Nb_sites = seq(0, min_ss),
                DAF = NA)
head(SW)
for (i in 1:nrow(SW)) {
  SW$DAF[i] = sum(colSums(no_ms_SW_matrix) == SW$Nb_sites[i])
}
head(SW)

# Plot the unfolded SFS
# On ne s'occupe pas du folded (qu'on utilise sur Minor Allele Frequency QUAND ONE CONNAIT PAS L'ETAT ANCESTRAL)
#barplot(SW$DAF, xlab="derived allele frequency", ylab="number of sites")

jpeg("SFS_SW_CHROMOSOME7_ORYZA.jpg", width = 1200, height = 800)
barplot(SW$DAF, xlab="derived allele frequency", ylab="number of sites")
title("SFS_SW Chromosome 7 - genes")
dev.off()
cat("Plot 'SFS_SW_CHROMOSOME7_ORYZA.jpg' created\n")

cat("#### Mutation S #### \n")

nsnps <- nrow(snp_counts_p80_s)
cat("Number of SNPs : ", nsnps, "\n")
S = matrix(NA, nrow = nind*2, ncol = nsnps)
# Pour chaque ligne, je cree un vecteur contenant
# n alleles derives (= 1), m alleles ancestraux (= 0) et l NAs
# ou n = snp_counts_p80_ws$Nb_derived
# m = snp_counts_p80_ws$Nb_Total - snp_counts_p80_ws$Nb_derived
# l = nind*2 - snp_counts_p80_ws$Nb_Total
for (i in 1:ncol(S)) {
    n = snp_counts_p80_s$Nb_Derived[i]
    m = snp_counts_p80_s$Nb_Total[i] - snp_counts_p80_s$Nb_Derived[i]
    l = nind*2 - snp_counts_p80_s$Nb_Total[i]
    S[,i] = c(rep(1, n), rep(0, m), rep(NA, l))
}
# md_persite <- colSums(is.na(WS))
# min_ss <- nind - max(md_persite)
# Go through each site, sample min_ss genotypes for each site without missing data
# S_no_missing <- apply(S, 2, function(each_site) {sample(each_site[!is.na(each_site)], size=min_ss, replace = F)})
# # Check that there is no missing data
# if(sum(is.na(S_no_missing))>0) {stop("Error: there is still missing data after resampling!")}
# # Get the SFS for the dataset without missing data
# # SW <- table(colSums(SW_no_missing))

no_ms_S_list <- lapply(1:ncol(S), function(i) {
  each_site <- S[, i]
  population <- sum(!is.na(each_site))
  sample_size <- min(min_ss, population)
  sample(each_site[!is.na(each_site)], size = sample_size, replace = FALSE)
})

# Ensure each sampled vector is of the same length
no_ms_S_matrix <- do.call(cbind, lapply(no_ms_S_list, function(x) {
  length_diff <- min_ss - length(x)
  if (length_diff > 0) {
    return(c(x, rep(NA, length_diff)))
  } else {
    return(x)
  }
}))

# Ensure no_ms_geno_matrix is numeric
no_ms_S_matrix <- apply(no_ms_S_matrix, 2, as.numeric)

# Check missing data
missing_data_indices <- which(is.na(no_ms_S_matrix), arr.ind = TRUE)
if (length(missing_data_indices) > 0) {
  # Remove columns with missing data
  no_ms_S_matrix <- no_ms_S_matrix[, colSums(is.na(no_ms_S_matrix)) == 0]
}
S = data.frame(Nb_sites = seq(0, min_ss),
                DAF = NA)
head(S)
for (i in 1:nrow(S)) {
  S$DAF[i] = sum(colSums(no_ms_S_matrix) == S$Nb_sites[i])
}
head(S)

# Plot the unfolded SFS
# On ne s'occupe pas du folded (qu'on utilise sur Minor Allele Frequency QUAND ONE CONNAIT PAS L'ETAT ANCESTRAL)
#barplot(S$DAF, xlab="derived allele frequency", ylab="number of sites")

jpeg("SFS_S_CHROMOSOME7_ORYZA.jpg", width = 1200, height = 800)
barplot(S$DAF, xlab="derived allele frequency", ylab="number of sites")
title("SFS_S Chromosome 7 - genes")
dev.off()
cat("Plot 'SFS_S_CHROMOSOME7_ORYZA.jpg' created\n")


cat("#### Mutation W #### \n")

nsnps <- nrow(snp_counts_p80_w)
cat("Number of SNPs : ", nsnps, "\n")
W = matrix(NA, nrow = nind*2, ncol = nsnps)
# Pour chaque ligne, je cree un vecteur contenant
# n alleles derives (= 1), m alleles ancestraux (= 0) et l NAs
# ou n = snp_counts_p80_ws$Nb_derived
# m = snp_counts_p80_ws$Nb_Total - snp_counts_p80_ws$Nb_derived
# l = nind*2 - snp_counts_p80_ws$Nb_Total
for (i in 1:ncol(W)) {
    n = snp_counts_p80_w$Nb_Derived[i]
    m = snp_counts_p80_w$Nb_Total[i] - snp_counts_p80_w$Nb_Derived[i]
    l = nind*2 - snp_counts_p80_w$Nb_Total[i]
    W[,i] = c(rep(1, n), rep(0, m), rep(NA, l))
}
# md_persite <- colSums(is.na(WS))
# min_ss <- nind - max(md_persite)
# Go through each site, sample min_ss genotypes for each site without missing data
# W_no_missing <- apply(W, 2, function(each_site) {sample(each_site[!is.na(each_site)], size=min_ss, replace = F)})
# # Check that there is no missing data
# if(sum(is.na(W_no_missing))>0) {stop("Error: there is still missing data after resampling!")}
# Get the SFS for the dataset without missing data
# SW <- table(colSums(SW_no_missing))

no_ms_W_list <- lapply(1:ncol(W), function(i) {
  each_site <- W[, i]
  population <- sum(!is.na(each_site))
  sample_size <- min(min_ss, population)
  sample(each_site[!is.na(each_site)], size = sample_size, replace = FALSE)
})

# Ensure each sampled vector is of the same length
no_ms_W_matrix <- do.call(cbind, lapply(no_ms_W_list, function(x) {
  length_diff <- min_ss - length(x)
  if (length_diff > 0) {
    return(c(x, rep(NA, length_diff)))
  } else {
    return(x)
  }
}))

# Ensure no_ms_geno_matrix is numeric
no_ms_W_matrix <- apply(no_ms_W_matrix, 2, as.numeric)

# Check missing data
missing_data_indices <- which(is.na(no_ms_W_matrix), arr.ind = TRUE)
if (length(missing_data_indices) > 0) {
  # Remove columns with missing data
  no_ms_W_matrix <- no_ms_W_matrix[, colSums(is.na(no_ms_W_matrix)) == 0]
}

W = data.frame(Nb_sites = seq(0, min_ss),
                DAF = NA)
head(W)
for (i in 1:nrow(W)) {
  W$DAF[i] = sum(colSums(no_ms_W_matrix) == W$Nb_sites[i])
}
head(W)


jpeg("SFS_W_CHROMOSOME7_ORYZA.jpg", width = 1200, height = 800)
barplot(W$DAF, xlab="derived allele frequency", ylab="number of sites")
title("SFS_W Chromosome 7 - genes")
dev.off()
cat("Plot 'SFS_W_CHROMOSOME7_ORYZA.jpg' created\n")

# Plot the unfolded SFS
# On ne s'occupe pas du folded (qu'on utilise sur Minor Allele Frequency QUAND ONE CONNAIT PAS L'ETAT ANCESTRAL)
#barplot(W$DAF, xlab="derived allele frequency", ylab="number of sites")


# jpeg("SFS.jpg", width = 1200, height = 800)
# jpeg("SFS_GENE.jpg", width = 1200, height = 800)
jpeg("SFS_ALL_MUT_CHROMOSOME7_ORYZA.jpg", width = 1200, height = 800)
par(mfrow = c(2,2))
barplot(SW$DAF, xlab="derived allele frequency", ylab="number of sites", main = paste0("SW (n = ", nrow(snp_counts_p80_sw), ")"))
barplot(WS$DAF, xlab="derived allele frequency", ylab="number of sites", main = paste0("WS (n = ", nrow(snp_counts_p80_ws), ")"))

barplot(W$DAF, xlab="derived allele frequency", ylab="number of sites", main = paste0("W (n = ", nrow(snp_counts_p80_w), ")"))
barplot(S$DAF, xlab="derived allele frequency", ylab="number of sites", main = paste0("S (n = ", nrow(snp_counts_p80_s), ")"))
par(mfrow = c(1,1))
dev.off()



# jpeg("SFS_density.jpg", width = 1200, height = 800)
# jpeg("SFS_density_GENE.jpg", width = 1200, height = 800)
jpeg("SFS_density_CHROMOSOME7_ORYZA.jpg", width = 1200, height = 800)
par(mfrow = c(2,2))
barplot(SW$DAF/sum(SW$DAF), xlab="derived allele frequency", ylab="Density", main = paste0("SW (n = ", nrow(snp_counts_p80_sw), ")"))
barplot(WS$DAF/sum(WS$DAF), xlab="derived allele frequency", ylab="Density", main = paste0("WS (n = ", nrow(snp_counts_p80_ws), ")"))

barplot(W$DAF/sum(W$DAF), xlab="derived allele frequency", ylab="Density", main = paste0("W (n = ", nrow(snp_counts_p80_w), ")"))
barplot(S$DAF/sum(S$DAF), xlab="derived allele frequency", ylab="Density", main = paste0("S (n = ", nrow(snp_counts_p80_s), ")"))
par(mfrow = c(1,1))
dev.off()


df = data.frame(mutation = rep(c("SW", "WS", "W", "S"), each = nrow(SW)),
                n_sites = c(SW$Nb_sites, WS$Nb_sites, W$Nb_sites, S$Nb_sites),
                density = c(SW$DAF/sum(SW$DAF), WS$DAF/sum(WS$DAF), W$DAF/sum(W$DAF), S$DAF/sum(S$DAF)))

write.table(df, "df_oryza_chrom7_genes.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# View(df)

library(ggplot2)
p = ggplot(df[which(df$mutation %in% c("WS", "SW")),], aes(x = n_sites, y = density, colour = mutation, fill = mutation)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme_bw()
p
# ggsave("SFS_WS_SW.jpg", p, width = 14, height =  10)
# ggsave("SFS_WS_SW_GENE.jpg", p, width = 14, height =  10)
ggsave("SFS_WS_SW_CHROMOSOME10_ORYZA.jpg", p, width = 14, height =  10)

df2 = data.frame(n_sites = SW$Nb_sites,
                ratio = SW$DAF/sum(SW$DAF)/WS$DAF/sum(WS$DAF))

ggplot(df2, aes(x = n_sites, y = ratio)) +
    geom_point() +
    theme_bw()


# Sous-jeux de donnees
# Tous les gènes




#-------------------------------------
# gBGC


source("/home/genouest/cnrs_umr6553/rdelage/scripts/gBGC_estimation_by_least_squares.R")
nrow(WS)
nrow(SW)
# The two SFS must have same length
# Comme c'est un ratio il faut enlever les 0
# Ce probleme sera moins grave avec le VCF non filtre car il n'aura probablement pas de 0
rm_zeros = which(WS$DAF == 0 | SW$DAF == 0)

# GENOME ENTIER
if(length(rm_zeros) > 0){
  gBGC = least_square(WS$DAF[- rm_zeros],SW$DAF[- rm_zeros])
}else{
  gBGC = least_square(WS$DAF,SW$DAF)
}
gBGC
# gBGC = least_square(WS$DAF[- rm_zeros],SW$DAF[- rm_zeros])
# gBGC = least_square(WS$DAF,SW$DAF)
# gBGC
# $B
#         x 
# 0.1486005 
# $mutbias
# (Intercept) 
#    1.358993 
# $e1
# 0 
# $e2
# 0.02921832 


# gBGC GENES
# $B
#         x 
# 0.1414339 
# $mutbias
# (Intercept) 
#    1.402046 
# $e1
# 0 
# $e2


# gBGC exon 1
# $B
#        x 
# 0.177427 
# $mutbias
# (Intercept) 
#    1.651414 
# $e1
# 0 
# $e2 
# 0





# gBGC gradients ?
# Ne maarche pas avec seulement 1 chromosome...
# colnames(gff)
# df = aggregate(weighted.mean.rho ~ rank, gff[which(gff$feature == "gene"),], median)
# df
# df$n_WS = NA
# df$n_SW = NA
# df$B = NA
# df$mutbias = NA
# df$e1 = NA
# df$e2 = NA
# head(df)
# 
# for (j in 1:nrow(df)) {
#   cat(j, "\n")
#   # ranges = gff_ranges[which(gff_ranges$feature == "exon" & gff_ranges$rank == df$rank[j] & gff_ranges$rank == df$nb_exons[j] & gff_ranges$nb_exons < 15)]
#   # ranges = gff_ranges[which(gff_ranges$feature == "exon" & gff_ranges$rank == df$rank[j] & gff_ranges$nb_exons < 15)]
#   ranges = gff_ranges[which(gff_ranges$feature == "gene")]
#   ranges
# 
#   # Find overlapping
#   hits = findOverlaps(snp_ranges, ranges)
#   length(hits)
# 
#   snp_counts_merged_subset = snp_counts_merged[queryHits(hits),]
# 
#   snp_counts_merged_subset$filter = unlist(lapply(1:nrow(snp_counts_merged_subset), function(x) {ifelse(max(snp_counts_merged_subset[x,c("P_trees_A", "P_trees_C", "P_trees_G", "P_trees_T")]) > 0.7, "PASS", "REMOVE")}))
#   snp_counts_merged_filtered_p80 <- snp_counts_merged_subset[which(snp_counts_merged_subset$filter == "PASS"),]
# 
#   ### DAF Computations --> On filtered data
#   # Nb derivated alleles for each row
#   snp_counts_merged_filtered_p80$Nb_Derived <- ifelse(snp_counts_merged_filtered_p80$Derived == "A", 
#                                                         snp_counts_merged_filtered_p80$count_A,
#                                                         ifelse(snp_counts_merged_filtered_p80$Derived == "C", 
#                                                               snp_counts_merged_filtered_p80$count_C,
#                                                               ifelse(snp_counts_merged_filtered_p80$Derived == "G", 
#                                                                       snp_counts_merged_filtered_p80$count_G,
#                                                                       snp_counts_merged_filtered_p80$count_T)))
# 
#   snp_counts_merged_filtered_p80$Nb_Total_Alleles <- snp_counts_merged_filtered_p80$N_CHR
# 
#   snp_counts_p80_sw <- snp_counts_merged_filtered_p80[which(snp_counts_merged_filtered_p80$Mutation == "SW"),]
#   snp_counts_p80_ws <- snp_counts_merged_filtered_p80[which(snp_counts_merged_filtered_p80$Mutation == "WS"),]
# 
#   # On refait la procedure pour avoir WS et SW
#   # number of SNPs
#   nsnps <- nrow(snp_counts_p80_ws)
#   # number of individuals
#   nind <- max(c(snp_counts_p80_ws$Nb_Total, snp_counts_p80_sw$Nb_Total), na.rm = T)/2
#   WS = matrix(NA, nrow = nind*2, ncol = nsnps)
#   # Pour chaque ligne, je cree un vecteur contenant
#   # n alleles derives (= 1), m alleles ancestraux (= 0) et l NAs
#   # ou n = snp_counts_p80_ws$Nb_derived
#   # m = snp_counts_p80_ws$Nb_Total - snp_counts_p80_ws$Nb_derived
#   # l = nind*2 - snp_counts_p80_ws$Nb_Total
#   for (i in 1:ncol(WS)) {
#       n = snp_counts_p80_ws$Nb_Derived[i]
#       m = snp_counts_p80_ws$Nb_Total[i] - snp_counts_p80_ws$Nb_Derived[i]
#       l = nind*2 - snp_counts_p80_ws$Nb_Total[i]
#       WS[,i] = c(rep(1, n), rep(0, m), rep(NA, l))
#   }
#   md_persite <- colSums(is.na(WS))
#   #min_ss <- nind - max(md_persite)
#   min_ss <- abs(nind-max(md_persite))
#   
#   # # Go through each site, sample min_ss genotypes for each site without missing data
#   # WS_no_missing <- apply(WS, 2, function(each_site) {sample(each_site[!is.na(each_site)], size=min_ss, replace = F)})
#   # # Check that there is no missing data
#   # if(sum(is.na(WS_no_missing))>0) {stop("Error: there is still missing data after resampling!")}
#   # # Get the SFS for the dataset without missing data
#   # # WS <- table(colSums(WS_no_missing))
#   
#   WS_no_missing_list <- lapply(1:ncol(WS), function(i) {
#     each_site <- WS[, i]
#     population <- sum(!is.na(each_site))
#     sample_size <- min(min_ss, population)
#     sample(each_site[!is.na(each_site)], size = sample_size, replace = FALSE)
#   })
#   
#   # Ensure each sampled vector is of the same length
#   WS_no_missing_matrix <- do.call(cbind, lapply(WS_no_missing_list, function(x) {
#     length_diff <- min_ss - length(x)
#     if (length_diff > 0) {
#       return(c(x, rep(NA, length_diff)))
#     } else {
#       return(x)
#     }
#   }))
#   
#   WS_no_missing_matrix <- apply(WS_no_missing_matrix, 2, as.numeric)
#   
#   
#   # Check and display missing data
#   missing_data_indices <- which(is.na(WS_no_missing_matrix), arr.ind = TRUE)
#   if (length(missing_data_indices) > 0) {
#     # Remove columns with missing data
#     WS_no_missing_matrix <- WS_no_missing_matrix[, colSums(is.na(WS_no_missing_matrix)) == 0]
#   }
#   
#   
#   WS = data.frame(Nb_sites = seq(0, min_ss),
#                   DAF = NA)
#   for (i in 1:nrow(WS)) {
#     WS$DAF[i] = sum(colSums(WS_no_missing) == WS$Nb_sites[i])
#   }
# 
# 
#   nsnps <- nrow(snp_counts_p80_sw)
#   SW = matrix(NA, nrow = nind*2, ncol = nsnps)
#   # Pour chaque ligne, je cree un vecteur contenant
#   # n alleles derives (= 1), m alleles ancestraux (= 0) et l NAs
#   # ou n = snp_counts_p80_ws$Nb_derived
#   # m = snp_counts_p80_ws$Nb_Total - snp_counts_p80_ws$Nb_derived
#   # l = nind*2 - snp_counts_p80_ws$Nb_Total
#   for (i in 1:ncol(SW)) {
#       n = snp_counts_p80_sw$Nb_Derived[i]
#       m = snp_counts_p80_sw$Nb_Total[i] - snp_counts_p80_sw$Nb_Derived[i]
#       l = nind*2 - snp_counts_p80_sw$Nb_Total[i]
#       SW[,i] = c(rep(1, n), rep(0, m), rep(NA, l))
#   }
#   # md_persite <- colSums(is.na(WS))
#   # min_ss <- nind - max(md_persite)
#   # Go through each site, sample min_ss genotypes for each site without missing data
#   # SW_no_missing <- apply(SW, 2, function(each_site) {sample(each_site[!is.na(each_site)], size=min_ss, replace = F)})
#   # # Check that there is no missing data
#   # if(sum(is.na(SW_no_missing))>0) {stop("Error: there is still missing data after resampling!")}
#   # # Get the SFS for the dataset without missing data
#   # # SW <- table(colSums(SW_no_missing))
#   
#   no_ms_SW_list <- lapply(1:ncol(SW), function(i) {
#     each_site <- SW[, i]
#     population <- sum(!is.na(each_site))
#     sample_size <- min(min_ss, population)
#     sample(each_site[!is.na(each_site)], size = sample_size, replace = FALSE)
#   })
#   
#   # Ensure each sampled vector is of the same length
#   no_ms_SW_matrix <- do.call(cbind, lapply(no_ms_SW_list, function(x) {
#     length_diff <- min_ss - length(x)
#     if (length_diff > 0) {
#       return(c(x, rep(NA, length_diff)))
#     } else {
#       return(x)
#     }
#   }))
# 
#   # Ensure no_ms_geno_matrix is numeric
#   no_ms_SW_matrix <- apply(no_ms_SW_matrix, 2, as.numeric)
# 
#   # Check and display missing data
#   missing_data_indices <- which(is.na(no_ms_SW_matrix), arr.ind = TRUE)
#   if (length(missing_data_indices) > 0) {
#     # Remove columns with missing data
#     no_ms_SW_matrix <- no_ms_SW_matrix[, colSums(is.na(no_ms_SW_matrix)) == 0]
#   }
# 
# 
#   SW = data.frame(Nb_sites = seq(0, min_ss),
#                   DAF = NA)
#   for (i in 1:nrow(SW)) {
#     SW$DAF[i] = sum(colSums(SW_no_missing) == SW$Nb_sites[i])
#   }
# 
#   rm_zeros = which(WS$DAF == 0 | SW$DAF == 0)
# 
#   if (length(rm_zeros) > 0) {
#     WS = WS[-rm_zeros,]
#     SW = SW[-rm_zeros,]
#   }
#   gBGC = least_square(WS$DAF,SW$DAF)
#   gBGC
#   df$n_SW[j] = nrow(snp_counts_p80_sw)
#   df$n_WS[j] = nrow(snp_counts_p80_ws)
#   df$B[j] = gBGC$B
#   df$mutbias[j] = gBGC$mutbias
#   df$e1[j] = gBGC$e1
#   df$e2[j] = gBGC$e2
# }
# 
# View(df)
# 