###################################################################
# Visualize the recombination gradients and the GC rate gradients #
###################################################################

setwd("/home/redelage/internship/data/gff_rho/") # Must to be changed according to the working directory

#' @title Load necessary packages
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

#' @title Read GFF
#' 
#' @description Allows to read GFF annotation file of the species
#' 
#' @param species (character) : The species which we wants to read the data
#'
#' @return None

read_gff <- function(species){
  if(species == "Arabidopsis_thaliana"){
    gff = readRDS("gff_rho_Arabidopsis_thaliana_1001genomes.rds")
  }else if(species == "Oryza_sativa"){
    gff = readRDS("gff_rho_Oryza_sativa_Wang2018.rds")
  }else if(species == "Sorghum_bicolor"){
    gff = readRDS("gff_rho_Sorghum_bicolor_Lozano2021.rds")
  }else if(species == "Populus_tremula"){
    gff = readRDS("gff_rho_Populus_tremula_Liu2022.rds")
  }
  cds = gff[which(gff$feature == "CDS" & gff$nb_exons < 15 & gff$rank < 15),]
  
  assign(paste0("gff_", species), gff, envir = .GlobalEnv)
  assign(paste0("cds_", species), cds, envir = .GlobalEnv)
}

read_gff("Arabidopsis_thaliana")
read_gff("Oryza_sativa")


#' @title Recombination gradient plot
#' 
#' @description Allows to plot the recombination gradients at CDS level 
#' (previously defined into the read GFF function) with the gradients for genes 
#' with sizes from 1 to 14 exons. The median recombination gradient is also plotted.
#' The plot will be stored into a PNG file format.
#' 
#' @param species (character) : The species which we wants to read the data
#' @param xlim (by default : c(1, 14)) : Scale for the x-axis.
#' @param ylim (by default : c(0, 1)) : Scale for the y-axis
#' @param legend (by default : TRUE) : Display the plot legend if the parameter 
#' is fixed on TRUE. Hide it otherwise
#'
#' @return None

recombination_gradient_plot <- function(species, xlim = c(1, 14), ylim = c(0, 1), legend = TRUE){
  species_plot_title <- gsub("_", " ", species)

  cds <- get(paste0("cds_", species))
  
  recomb_rank = aggregate(weighted.mean.rho ~ rank + nb_exons, cds , median)
  recomb_rank$nb_exons = as.factor(recomb_rank$nb_exons)
  
  p1 = ggplot(recomb_rank, aes(x = rank, y = weighted.mean.rho, group = nb_exons, colour = nb_exons)) +
    scale_color_viridis_d() +
    geom_line() +
    geom_point() +
    xlim(xlim) +
    ylim(ylim) +
    xlab("Gene length (exons)") + ylab("Median recombination rate (ρ/kb)") +
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
recombination_gradient_plot("Oryza_sativa", ylim = c(0.045, 0.125))

#' @title GC gradient plot
#' 
#' @description Allows to plot the GC content gradients at CDS level 
#' (previously defined into the read GFF function) with the gradients for genes 
#' with sizes from 1 to 14 exons. GC content corresponds to that at the 3rd codon position.
#' The median GC content gradient is also plotted.
#' The plot will be stored into a PNG file format.
#' 
#' @param species (character) : The species which we wants to read the data.
#' @param xlim (by default : c(1, 14)) : Scale for the x-axis.
#' @param ylim (by default : c(0, 1)) : Scale for the y-axis.
#' @param legend (by default : TRUE) : Display the plot legend if the parameter 
#' is fixed on TRUE. Hide it otherwise
#'
#' @return None

gc_gradient <- function(species, xlim = c(1,14), ylim = c(0,1), legend = TRUE){
  species_plot_title = gsub("_", " ", species)
  cds <- get(paste0("cds_", species))
  
  recomb_rank = aggregate(weighted.mean.rho ~ rank + nb_exons, cds , median)
  recomb_rank$nb_exons = as.factor(recomb_rank$nb_exons)
  
  gc_rank <- aggregate(gc3 ~ rank + nb_exons, cds, median)
  gc_rank$nb_exons <- as.factor(gc_rank$nb_exons)
  summary(gc_rank)
  
  p2 = ggplot(gc_rank, aes(x = rank, y = gc3, group = nb_exons, colour = nb_exons)) +
    scale_color_discrete(breaks = unique(recomb_rank$nb_exons), 
                         labels = unique(recomb_rank$nb_exons), 
                         aesthetics = "colour") +
    scale_colour_viridis_d(option = "inferno") +
    geom_line() +
    geom_point() +
    xlim(1, 14) +
    ylim(0, 1) +
    xlab("Gene length (exons)") + ylab("GC content (3rd codon positions)") +
    scale_alpha_manual(values = c(0.4, 1)) +
    guides(alpha = "none")+
    theme_classic()
  
  df2_avg = aggregate(gc3 ~ rank , gc_rank , median)
  
  p2 = p2 + geom_line(data = df2_avg, aes(x = rank, y = gc3, group = "black", colour = "black"), colour = "black", linewidth = 1.5)
  
  if(!legend){
    p2 = p2 + theme(legend.position = "none")
  } else {
    p2 = p2 + 
      labs(colour = "# exons", title = species_plot_title) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  
  
  ggsave(paste0(species,"_gc3_gradient.png"), p2, width = 8, height = 8)
}

gc_gradient("Arabidopsis_thaliana")
gc_gradient("Oryza_sativa", legend = FALSE)




# ---- Display all plots ----

grad_plots = grid.arrange(p1, p1_oryza, p2, p2_oryza, ncol = 2)

ggsave("Recombination_vs_GC3_gradients.jpg", grad_plots, width = 9, height = 10)

######################################


