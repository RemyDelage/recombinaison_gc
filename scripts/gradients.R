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
#' @param species (character) : The studied species.
#' @param xlim (by default : c(1, 14)) : Scale for the x-axis.
#' @param ylim (by default : c(0, 1)) : Scale for the y-axis.
#' @param legend (by default : TRUE) : Display the plot legend if the parameter 
#' is fixed on TRUE. Hide it otherwise.
#' @param title (by default : TRUE) : Display the plot title (i.e. the species name) 
#' if the parameter is fixed on TRUE. Hide it otherwise.
#'
#' @return None

recombination_gradient_plot <- function(species, xlim = c(1, 14), ylim = c(0, 1), legend = TRUE, title = TRUE){
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
    scale_alpha_manual(values = c(0.4, 1)) +
    guides(alpha = "none")+
    theme_classic()
    
  df_avg = aggregate(weighted.mean.rho ~ rank , recomb_rank , median)
  
  p1 = p1 + geom_line(data = df_avg, aes(x = rank, y = weighted.mean.rho, group = "black", colour = "black"), colour = "black", linewidth = 1.5)

  if(!legend){
    p1 = p1 + theme(legend.position = "none")
  } else {
    p1 = p1 + 
      labs(colour = "# exons")
  }
  
  if(title){
    p1 = p1 + 
      labs(title = species_plot_title) + 
      theme(plot.title = element_text(hjust = 0.5, face = "italic"))
  }
  
  p1
  
  ggsave(paste0(species,"_recombination_gradient.png"), p1, width = 8, height = 8)
}

recombination_gradient_plot("Arabidopsis_thaliana", title = TRUE)
recombination_gradient_plot("Oryza_sativa", ylim = c(0.045, 0.125), legend = TRUE, title = TRUE)


#' @title GC gradient plot
#' 
#' @description Allows to plot the GC content gradients at CDS level 
#' (previously defined into the read GFF function) with the gradients for genes 
#' with sizes from 1 to 14 exons. GC content corresponds to that at the 3rd codon position.
#' The median GC content gradient is also plotted.
#' The plot will be stored into a PNG file format.
#' 
#' @param species (character) : The studied species.
#' @param xlim (by default : c(1, 14)) : Scale for the x-axis.
#' @param ylim (by default : c(0, 1)) : Scale for the y-axis.
#' @param legend (by default : TRUE) : Display the plot legend if the parameter 
#' is fixed on TRUE. Hide it otherwise.
#' @param title (by default : TRUE) : Display the plot title (i.e. the species name) 
#' if the parameter is fixed on TRUE. Hide it otherwise.
#'
#' @return None

gc_gradient_plot <- function(species, xlim = c(1,14), ylim = c(0,1), legend = TRUE, title = TRUE){
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
    xlim(xlim) +
    ylim(ylim) +
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
      labs(colour = "# exons")
  }
  
  if(title){
  p2 = p2 +
    labs(title = species_plot_title) +
    theme(plot.title = element_text(hjust = 0.5, face = "italic"))
  }
  
  p2
  
  ggsave(paste0(species,"_gc3_gradient.png"), p2, width = 8, height = 8)
}

gc_gradient_plot("Arabidopsis_thaliana", title = FALSE, legend = FALSE)
gc_gradient_plot("Oryza_sativa", title = FALSE, legend = TRUE)

Arabidopsis_thaliana_gc_grad_plot
Oryza_sativa_gc_grad_plot

#' @title Plots comparison
#' 
#' @description Allows to display all the plots of the recombination gradients and 
#' GC gradients between two species to compare them.
#' 
#' @param species1 (character) : The first species to compare.
#' @param species2 (character) : The second species to compare.
#' @param ylim_recomb_sp1 (by default : NULL) : Scale for the y-axis for the 
#' recombination gradients of the first species. The default value taken into account 
#' is that of the "Recombination gradient plot" function if this parameter is set to NULL.
#' @param ylim_gc_sp1 (by default : NULL) : Scale for the y-axis for the 
#' GC gradients of the first species. The default value taken into account 
#' is that of the "GC gradient plot" function if this parameter is set to NULL.
#' @param ylim_recomb_sp2 (by default : NULL) : Scale for the y-axis for the 
#' recombination gradients of the first species. The default value taken into account 
#' is that of the "Recombination gradient plot" function if this parameter is set to NULL.
#' @param ylim_gc_sp2 (by default : NULL) : Scale for the y-axis for the 
#' GC gradients of the first species. The default value taken into account 
#' is that of the "GC gradient plot" function if this parameter is set to NULL.
#'
#' @return None


plots_comparison <- function(species1, species2, ylim_recomb_sp1 = NULL, ylim_gc_sp1 = NULL,
                             ylim_recomb_sp2 = NULL, ylim_gc_sp2 = NULL){
  species1_plot_title = gsub("_", " ", species1)
  species2_plot_title = gsub("_", " ", species2)

  recomb_plot_species1 <- recombination_gradient_plot(species1, legend = FALSE, title = TRUE)
  recomb_plot_species2 <- recombination_gradient_plot(species2, legend = TRUE, title = TRUE)
  gc_plot_species1 <- gc_gradient_plot(species1, legend = FALSE, title = FALSE)
  gc_plot_species2 <- gc_gradient_plot(species2, legend = TRUE, title = FALSE)
  
  if(!is.null(ylim_recomb_sp1)){
    recomb_plot_species1 <- recombination_gradient_plot(species1, ylim = ylim_recomb_sp1, legend = TRUE, title = TRUE)
  }

  if(!is.null(ylim_recomb_sp2)){
    recomb_plot_species2 <- recombination_gradient_plot(species2, ylim = ylim_recomb_sp2, legend = TRUE, title = TRUE)
  }

  if(!is.null(ylim_gc_sp1)){
    gc_plot_species1 <- gc_gradient_plot(species1, ylim = ylim_gc_sp1, legend = TRUE, title = FALSE)
  }
  
  if(!is.null(ylim_gc_sp2)){
    gc_plot_species2 <- gc_gradient_plot(species2, ylim = ylim_gc_sp2, legend = TRUE, title = FALSE)
  }


  grad_plots = grid.arrange(recomb_plot_species1, recomb_plot_species2,
                            gc_plot_species1, gc_plot_species2, ncol = 2)

  ggsave(paste0(species1, "_vs_", species2, "_gradients_comparison.png"), grad_plots, width = 9, height = 10)
}

plots_comparison("Arabidopsis_thaliana", "Oryza_sativa", ylim_recomb_sp2 = c(0.045, 0.125), ylim_gc_sp1 = c(0.25, 1), ylim_gc_sp2 = c(0.25, 1))



