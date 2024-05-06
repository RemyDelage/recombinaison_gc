install.packages("tidyverse")
library(dplyr)
library(ggplot2)
at_rds <- readRDS("gff_rho_Arabidopsis_thaliana_1001genomes.rds")
exons <- at_rds %>% filter(feature == "exon")
exons$gc_percentage <- exons$gc_count/(exons$gc_count+exons$at_count)*100

exons$gc_gradient <- c(NA, diff(exons$gc_percent)) / c(NA, diff(exons$start))

exons_clean <- exons %>%
  filter(!is.infinite(gc_gradient) & !is.na(gc_gradient))

mean_gc_gradient <- exons_clean %>%
  group_by(chromosome) %>%
  summarise(mean_gc_gradient = mean(gc_gradient, na.rm = TRUE))

exons_clean_with_mean <- exons_clean %>%
  mutate(strand_numeric = ifelse(strand == "+", 1, -1),
         exon_position = cumsum(start * strand_numeric) - start * strand_numeric + 1) %>%
  left_join(mean_gc_gradient, by = "chromosome")

gc_grad_per_chrom <- ggplot(data = exons_clean_with_mean, aes(x = exon_position, y = gc_gradient, color = chromosome)) +
  geom_line() +
  geom_hline(data = mean_gc_gradient, 
             aes(yintercept = mean_gc_gradient, linetype = "Mean GC gradient"),
             color = "black") +
  facet_wrap(~chromosome, scales = "free") +
  labs(x = "Position", y = "GC gradient", color = "Chromosome",
       title = "GC gradients for each chromosome of Arabidopsis thaliana") +
  theme_classic()

print(gc_grad_per_chrom)

ggsave("test_gc_gradients_per_chromosome.png", plot = gc_grad_per_chrom)
