#####################################################################
# User script for plotting recombination gradients and GC gradients #
#####################################################################

# Load the script which contains the functions to make the plots
#(must be on the same directory than this script)

source("gradients.R")

# --------------------- #
# Parameters to modify
# --------------------- #

# Directory where the GFF files are stored
directory = "/home/redelage/internship/data/gff_rho/"

# Species to study
species1 = "Arabidopsis_thaliana"
species2 = "Oryza_sativa"

# File where the plots will be stored
output_file = "Arabidopsis_Oryza_gradients.pdf"

# Scale for the y-axis of the plots
ylim_recomb_sp1 = c(0, 1)
ylim_recomb_sp2 = c(0.045, 0.125)
ylim_gc_sp1 = c(0.25, 1)
ylim_gc_sp2 = c(0.25, 1)

# ------------------------------------------------------- #
# Create and display the plots (You don't need to modify)
# ------------------------------------------------------- #

setwd(directory)

read_gff(species1)
read_gff(species2)

generate_and_save_plots(species1 = species1, species2 = species2,
                        ylim_recomb_sp1 = ylim_recomb_sp1, 
                        ylim_recomb_sp2 = ylim_recomb_sp2 , 
                        ylim_gc_sp1 = ylim_gc_sp1, 
                        ylim_gc_sp2 = ylim_gc_sp2, 
                        output_file = output_file)

