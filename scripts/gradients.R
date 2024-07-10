###################################################################
# Visualize the recombination gradients and the GC rate gradients #
###################################################################

setwd("/home/redelage/internship/data/gff_rho/") # Must to be changed according to the working directory

#'@title Load necessary packages
#' 
#' @description Allows to load all the necessary packages for process data and estimate gBGC values
#' 
#'
#' @return None

load_packages <- function(){
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(viridis) 
}

load_packages()

#'@title Read GFF
#' 
#' @description Allows to read GFF annotation file of the species
#' 
#' @param species (character) : The species which we wants to read the data
#'
#' @return None

read_gff <- function(species){
  if(species == "Arabidopsis_thaliana"){
    gff = readRDS("gff_rho_Arabidopsis_thaliana_1001genomes.rds")
  }
  if(species == "Oryza_sativa"){
    gff = readRDS("gff_rho_Oryza_sativa_Wang2018.rds")
  }
  if(species == "Sorghum_bicolor"){
    gff = readRDS("gff_rho_Sorghum_bicolor_Lozano2021.rds")
  }
  if(species == "Populus_tremula"){
    gff = readRDS("gff_rho_Populus_tremula_Liu2022.rds")
  }
  cds = gff[which(gff$feature == "CDS" & gff$nb_exons < 15 & gff$rank < 15),]
  
  assign(paste0("gff_", species), gff, envir = .GlobalEnv)
  assign(paste0("cds_", species), cds, envir = .GlobalEnv)
}

read_gff("Arabidopsis_thaliana")
read_gff("Oryza_sativa")


### Recombination gradient

recombination_gradient_plot <- function(species, frontsize = 16, dotsize = 0.2, linesize = 1.5){
  species_plot_title <- gsub("_", " ", species)

  cds <- get(paste0("cds_", species))
  typeof(cds)
  
  recomb_rank = aggregate(weighted.mean.rho ~ rank + nb_exons, cds , median)
  recomb_rank$nb_exons = as.factor(recomb_rank$nb_exons)
  
  p1 = ggplot(recomb_rank, aes(x = rank, y = weighted.mean.rho, group = nb_exons, colour = nb_exons)) +
    scale_color_viridis_d() +
    geom_line() +
    geom_point() +
    xlim(1, 14) +
    ylim(0, 1) +
    xlab("CDS part rank") + ylab("Median recombination rate (ρ/kb)") +
    labs(colour = "# exons", title  = species_plot_title) +
    scale_alpha_manual(values = c(0.4, 1)) +
    guides(alpha = "none")+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5))
  
  df_avg = aggregate(weighted.mean.rho ~ rank , recomb_rank , median)
  
  p1 = p1 + geom_line(data = df_avg, aes(x = rank, y = weighted.mean.rho, group = "black", colour = "black"), colour = "black", linewidth = 1.5)
  p1
  
  ggsave(paste0(species,"_recombination_gradient.png"), p1, width = 8, height = 8)
}

recombination_gradient_plot("Arabidopsis_thaliana")





### GC3 gradient

gc_rank <- aggregate(gc3 ~ rank + nb_exons, cds, median)
summary(gc_rank)
gc_rank$nb_exons <- as.factor(gc_rank$nb_exons)

p2 = ggplot(gc_rank, aes(x = rank, y = gc3, group = nb_exons, colour = nb_exons)) +
  scale_color_discrete(breaks = unique(recomb_rank$nb_exons), 
                       labels = unique(recomb_rank$nb_exons), 
                       aesthetics = "colour") +
  scale_colour_viridis_d(option = "inferno") +
  geom_line() +
  geom_point() +
  xlim(1, 14) +
  ylim(0, 1) +
  xlab("CDS part rank") + ylab("GC3 rate") +
  labs(colour = "# exons") +
  scale_alpha_manual(values = c(0.4, 1)) +
  guides(alpha = "none")+
  theme_classic()+
  theme(legend.position = "none")

p2

df2_avg = aggregate(gc3 ~ rank , gc_rank , median)

p2 = p2 + geom_line(data = df2_avg, aes(x = rank, y = gc3, group = "black", colour = "black"), colour = "black", linewidth = 1.5)
p2

# ---- ORYZA SATIVA ---- 

gff_oryza <- readRDS("gff_rho_Oryza_sativa_Wang2018.rds")
head(gff_oryza)
cds_oryza = gff_oryza[which(gff_oryza$feature == "CDS" & gff_oryza$nb_exons < 15 & gff_oryza$rank < 15),]

### Recombination gradient

recomb_rank_oryza = aggregate(weighted.mean.rho ~ rank + nb_exons, cds_oryza, median)
recomb_rank_oryza$nb_exons = as.factor(recomb_rank_oryza$nb_exons)

fontsize = 16
dotsize = 0.2
linesize = 1.5

p1_oryza = ggplot(recomb_rank_oryza, aes(x = rank, y = weighted.mean.rho, group = nb_exons, colour = nb_exons)) +
  scale_color_viridis_d() +
  geom_line() +
  geom_point() +
  xlim(1, 14) +
  ylim(0, 1)+
  xlab("CDS part rank") + ylab("Median recombination rate (ρ/kb)") +
  labs(colour = "# exons", title  = "Oryza sativa") +
  scale_alpha_manual(values = c(0.4, 1)) +
  guides(alpha = "none")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))

df_avg_oryza = aggregate(weighted.mean.rho ~ rank , recomb_rank_oryza , median)

p1_oryza = p1_oryza + geom_line(data = df_avg_oryza, aes(x = rank, y = weighted.mean.rho, group = "black", colour = "black"), colour = "black", linewidth = 1.5)
p1_oryza


### GC3 gradient

gc_rank_oryza <- aggregate(gc3 ~ rank + nb_exons, cds_oryza, median)
gc_rank_oryza$nb_exons <- as.factor(gc_rank_oryza$nb_exons)

p2_oryza = ggplot(gc_rank_oryza, aes(x = rank, y = gc3, group = nb_exons, colour = nb_exons)) +
  scale_color_discrete(breaks = unique(recomb_rank_oryza$nb_exons), 
                       labels = unique(recomb_rank_oryza$nb_exons), 
                       aesthetics = "colour") +
  scale_colour_viridis_d(option = "inferno") +
  geom_line() +
  geom_point() +
  xlim(1, 14) +
  ylim(0, 1) +
  xlab("CDS part rank") + ylab("GC3 rate") +
  labs(colour = "# exons") +
  scale_alpha_manual(values = c(0.4, 1)) +
  guides(alpha = "none")+
  theme_classic()

df2_avg_oryza = aggregate(gc3 ~ rank , gc_rank_oryza , median)

p2_oryza = p2_oryza + geom_line(data = df2_avg_oryza, aes(x = rank, y = gc3, group = "black", colour = "black"), colour = "black", linewidth = 1.5)
p2_oryza






# ---- Display all plots ----

grad_plots = grid.arrange(p1, p1_oryza, p2, p2_oryza, ncol = 2)

ggsave("Recombination_vs_GC3_gradients.jpg", grad_plots, width = 9, height = 10)

######################################


