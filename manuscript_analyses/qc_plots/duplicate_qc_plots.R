#!/usr/bin/env Rscript

library(cowplot)
library(dplyr)
library(ggplot2)
library(optparse)
library(stringr)
library(tidyr)

option_list = list(
  make_option(c("-a", "--mtdnaserver_vl_changes"), type = "character", default = NULL,
              help = "The vl_changes.txt file for mtDNA-Server results", metavar = "character"),
  make_option(c("-b", "--mutect_vl_changes"), type = "character", default = NULL,
              help = "The vl_changes.txt file for Mutect results", metavar = "character"),
  make_option(c("-d", "--plot_directory"), type = "character", default = NULL,
              help = "Directory to which plots should be written", metavar = "character"),
  make_option(c("-t", "--min_het_threshold"), type = "numeric", default = NULL,
              help = "Minimum heteroplasmy level to define a variant as a PASS heteroplasmic variant, filter out call if below", metavar = "float"))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

pdf(NULL)
# NOTE: Plots are for pairs where at least one sample of the pair had a called variant: (0,0) values are not in the dataframe

# Read in data and set the out directory for plots
mtdnaserver <- read.table(opt$mtdnaserver_vl_changes, sep = '\t', quote = "", header = TRUE)
mutect <- read.table(opt$mutect_vl_changes, sep = '\t', quote = "", header = TRUE)
setwd(opt$plot_directory)

# Combine mtDNA-Server and Mutect dataframes
mtdnaserver <- mutate(mtdnaserver, caller = "MtDNA-Server")
mutect <- mutate(mutect, caller = "Mutect")
all <- rbind(mtdnaserver, mutect)

# Filter to samples in both Mutect and mtDNA-Server dataframes (needed because several samples in mtDNA-Server failed to run and want to drop those in the comparison)
all <- filter(all, s1 %in% mutect$s1, s2 %in% mutect$s2, s1 %in% mtdnaserver$s1, s2 %in% mtdnaserver$s2)

# Double check that at least one sample of the pair has a heteroplasmy level above that of the min_het_threshold
# NOTE: will NOT round s1 and s2 afs before comparing to min_het_threshold
all <- filter(all, s1_af > opt$min_het_threshold | s2_af > opt$min_het_threshold)

# Filter to PASS variants (non artifact-prone sites)
all <- filter(all, filter == "PASS")

# Count number of distinct pairs for each comparison: duplicate pairs, mother-child (mc) or Mutect/mtDNA-Server (mm)
# NOTE: at moment only plotting dups
num_dups <- nrow(distinct(filter(all, dataset == "dups"), s1, s2))
num_mc <- nrow(distinct(filter(all,dataset == "mc"), s1, s2))
num_mm <- nrow(distinct(filter(all, dataset == "mm"), s1, s2))

# Format labels for plots
dup_label = paste("Duplicate Pairs (n=", num_dups, " pairs)", sep = "")
mc_label = paste("Mother-Child Pairs (n=", num_mc, " pairs)", sep = "")
mm_label = paste("Mutect vs mtDNA-Server Pairs (n=", num_mm, " pairs)", sep = "")

all$dataset_name <- factor(all$dataset,
                           levels = c("dups", "mc", "mm"),
                           labels = c(dup_label, mc_label, mm_label))


#####################################################
# Generate functions for plotting data
#####################################################
plot_scatter_comparision_for_pairs <- function(df){
  # Plot allele frequency of one sample vs that of its pair as a scatter plot
  # Input 'df' is dataframe to be plotted
  # 'df' should contain 's1_af' for allele frequency of one sample of a given pair
  # and 's2_af' for allele frequency of the corresponding sample in the given pair
  
  snps_plot <- ggplot(df, aes(x = s1_af, y = s2_af)) +
    geom_point(size = 2.0, stroke = 0) +
    geom_abline(slope = 1, linetype = "dashed", colour = "black") +
    theme_bw() +
    theme(strip.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
          strip.text.x = element_text(size = 18, colour = "black"),
          panel.spacing = unit(1.5, "lines"),
          axis.text = element_text(colour = "black", size = 12),
          axis.title = element_blank())
  
  return(snps_plot)
}

plot_density_for_pairs <- function(df){
  # Plot density of allele frequency values 
  # Input 'df' is dataframe to be plotted
  # 'df' should contain 's1_af' for allele frequency of one sample of a given pair to use for displaying density
  
  density_plot <- ggplot(df, aes(x = s1_af)) +
    geom_density() +
    theme_bw() +
    theme(axis.text.y = element_text(colour = "black", size = 12),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(colour = "black", size = 14, face = "bold", hjust = 0.5),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  
  return(density_plot)
}


#####################################################
# Plot duplicate comparison for SNPs
#####################################################
# Pull out duplicate pairs data
dups <- filter(all, dataset == "dups")

# Separate variant column into position, reference, alternate to help identify SNPs vs Indels
dups$pos <- stringr::str_split_fixed(dups$variant, "_", 4)[,2]
dups$ref <- stringr::str_split_fixed(dups$variant, "_", 4)[,3]
dups$alt <- stringr::str_split_fixed(dups$variant, "_", 4)[,4]
dups <- mutate(dups, pos = as.numeric(pos))

snps_only <- filter(dups, nchar(ref) == 1 & nchar(alt) == 1)

# Plot Mutect SNPs
mutect <- filter(snps_only, caller == "Mutect")
mutect_snps_plot <- plot_scatter_comparision_for_pairs(mutect) + theme(plot.margin = unit(c(0.25, 0.25, 1.5, 1.5), "cm"))

# Plot Mutect density
mutect_density_plot <- plot_density_for_pairs(mutect) + labs(y = "Density", title = "Mutect")

# Combine plots for Mutect density and SNPs
mutect_combined_plot <- plot_grid(mutect_density_plot, mutect_snps_plot, align = "v", ncol = 1, rel_heights = c(1, 5))

# Plot mtDNA-Server SNPs
mtdnaserver <- filter(snps_only, caller == "MtDNA-Server")
mtdnaserver_snps_plot <- plot_scatter_comparision_for_pairs(mtdnaserver) + theme(plot.margin = unit(c(0.25, 1.5, 1.5, 0.25), "cm"))

# Plot mtDNA-Server density
mtdnaserver_density_plot <- plot_density_for_pairs(mtdnaserver) + labs(y = "Density", title = "MtDNA-Server")

# Combine plots for mtDNA-Server density and SNPs
mtdnaserver_combined_plot <- plot_grid(mtdnaserver_density_plot, mtdnaserver_snps_plot, align = "v", ncol = 1, rel_heights = c(1,5))

# Combine plots for Mutect SNPs, Mutect density, mtDNA-Server SNPs, mtDNA-Server density
combined_plot <- plot_grid(mutect_combined_plot, mtdnaserver_combined_plot, align = "h")

combined_plot <- combined_plot + draw_label("Heteroplasmic Level\nSample 1", x = .50, y = .04, hjust = .5, vjust = .5,
                                       fontfamily = "", fontface = "bold", colour = "black", size = 14,
                                       lineheight = 0.9, alpha = 1) 

combined_plot <- combined_plot + draw_label("Heteroplasmic Level\nSample 2", x = .035, y = 0.50, hjust = .5, vjust = .5,
                                        fontfamily = "", fontface = "bold", colour = "black", size = 14,
                                        angle = 90, lineheight = 0.9, alpha = 1)

combined_plot <- combined_plot + draw_label("Density", x = .035, y = 0.90, hjust = .5, vjust = .5,
                                            fontfamily = "", fontface = "bold", colour = "black", size = 14,
                                            angle = 90, lineheight = 0.9, alpha = 1)

combined_plot

ggsave(combined_plot, filename = "scatter_dup_snp.png", dpi = 300, width = 10, height = 6, units = "in")
ggsave(combined_plot, filename = "scatter_dup_snp.pdf", dpi = 300, width = 10, height = 6, units = "in")


#####################################################
# Count number of homoplasmic SNPs
#####################################################
# For SNPs in duplicate samples in Mutect, calculate the fraction that are homoplasmic for the paper
s1_hom <- nrow(filter(mutect, s1_af >= 0.95))
s2_hom <- nrow(filter(mutect, s2_af >= 0.95))
s1_called <- nrow(filter(mutect, s1_af > 0))
s2_called <- nrow(filter(mutect, s2_af > 0))
print("Homoplasmic variants (SNP duplicates Mutect):")
s1_hom + s2_hom
print("Variants called (SNP duplicates Mutect):")
s1_called + s2_called
print("Fraction homoplasmic (SNP duplicates Mutect):")
(s1_hom + s2_hom)/(s1_called + s2_called)


#####################################################
# Plot SNPs and indels together for just Mutect
#####################################################
# Do not run this analysis for mtDNA-Server since indel-calling is in beta (analyses for mtDNA-Server should keep SNPs and indels separate)
mutect <- filter(dups, caller == "Mutect")

mutect_snps_and_indels_plot <- plot_scatter_comparision_for_pairs(mutect) +
  theme(axis.title = element_text(size = 14, face = "bold")) +
  labs(x = "Heteroplasmic Level\nSample 1", y = "Heteroplasmic Level\nSample 2")

mutect_snps_and_indels_density_plot <- plot_density_for_pairs(mutect) +
  theme(axis.title.y = element_text(size = 14, face = "bold")) +
  labs(y = "Density", title = "Mutect")

combined_plot <- plot_grid(mutect_snps_and_indels_density_plot, mutect_snps_and_indels_plot, ncol = 1 , align = "v" , rel_heights = c(1,5))
combined_plot
ggsave(combined_plot, filename = "mutect_snps_and_indels.png", dpi = 300, width = 8, height = 8, units = "in")
ggsave(combined_plot, filename = "mutect_snps_and_indels.pdf", dpi = 300, width = 8, height = 8, units = "in")


#####################################################
# Count number of homoplasmic variants for Mutect
#####################################################
s1_hom <- nrow(filter(mutect, s1_af >= 0.95))
s2_hom <- nrow(filter(mutect, s2_af >= 0.95))
s1_called <- nrow(filter(mutect, s1_af > 0))
s2_called <- nrow(filter(mutect, s2_af > 0))
print("Homoplasmic variants (SNP and indel duplicates Mutect):")
s1_hom + s2_hom
print("Variants called (SNP and indel duplicates Mutect):")
s1_called + s2_called
print("Fraction homoplasmic (SNP and indel duplicates Mutect):")
(s1_hom + s2_hom)/(s1_called + s2_called)


##########################################################################
# Plot duplicate comparison for indels (for both Mutect and mtDNA-Server)
##########################################################################
indels_only <- filter(dups, nchar(ref) != 1 | nchar(alt) != 1)

indel_comparison_plot <- plot_scatter_comparision_for_pairs(indels_only) +
  facet_wrap(~caller) +
  geom_point(alpha = .35, size = 2.0, stroke = 0) +
  theme(axis.title = element_text(colour = "black", size = 18, face = "bold")) +
  labs(x = "Heteroplasmic Level\nSample 1", y = "Heteroplasmic Level\nSample 2", title = dup_label)
indel_comparison_plot                               

ggsave(indel_comparison_plot, filename = "scatter_dup_indel.png", dpi = 300, width = 10, height = 6, units = "in")


#####################################################
# Generate jaccard index table
#####################################################
add_jaccard_results <- function(intersect_table, union_table, jaccard_table, row_name) {
  # intersect_table is dataframe containing only the intersect values
  # union_table is dataframe containing only the union values
  # jaccard_table is two-column table to which results should be added
  # row_name is string to use in first column of the row
  
  jaccard_index <- round(nrow(intersect_table)/nrow(union_table), 3) * 100
  f1 = paste(nrow(intersect_table), "/", nrow(union_table)," (", jaccard_index, "%)", sep = "")
  jaccard_table[nrow(jaccard_table) + 1,] <- c(row_name, f1)
  
  return(jaccard_table)
  }

generate_jaccard_table <- function(qc_type, method, name, data){
  # Function to calculate jaccard index (intersect/union) for:
  # total shared variants
  # total shared variants > 0.10
  # total shared variants > 0.50
  # total shared variants > 0.90
  # total shared variants between 0.10 and 0.90
  
  # qc_type is either "dups", "mc", or "mm"
  # method is either "MtDNA-Server" or "Mutect"
  # name is "snps", "indels", or "both"
  # data is a snp-only, indel-only, or complete dataframe
  
  # Filter to the specified QC type 
  df <- filter(data, dataset == qc_type)
  
  # Filter to the specified caller 
  df <- filter(df, caller == method)
  
  # Set up dataframe to store results
  jaccard_table <- as.data.frame(matrix(0, ncol = 2, nrow = 0))
  colnames(jaccard_table) <- c("Metric", "F1")

  # Make the first row contain the column names
  jaccard_table[1,] <- c("Caller", paste(method, qc_type, name, sep="_"))

  # Add a row of jaccard results for each heteroplasmy level comparison
  # Want to display both the fraction and calculate the percentage
  
  # Calculate jaccard index for total shared variants
  intersect_all <- filter(df, s1_af > 0 & s2_af > 0)
  union_all <- filter(df, s1_af > 0 | s2_af > 0)
  jaccard_table <- add_jaccard_results(intersect_all, union_all, jaccard_table, "Total Shared Variants")

  # Calculate jaccard index for total shared variants > 0.90
  intersect_90 <- filter(df, s1_af > 0.90 & s2_af > 0.90)
  union_90 <- filter(df, s1_af > 0.90 | s2_af > 0.90)
  jaccard_table <- add_jaccard_results(intersect_90, union_90, jaccard_table, "Total Shared Variants > 0.90")
  
  # Calculate jaccard index for total shared variants > 0.50
  intersect_50 <- filter(df, s1_af > 0.50 & s2_af > 0.50)
  union_50 <- filter(df, s1_af > 0.50 | s2_af > 0.50)
  jaccard_table <- add_jaccard_results(intersect_50, union_50, jaccard_table, "Total Shared Variants > 0.50")
  
  # Calculate jaccard index for total shared variants > 0.10
  intersect_10 <- filter(df, s1_af > 0.1 & s2_af > 0.1)
  union_10 <- filter(df, s1_af > 0.1 | s2_af > 0.1)
  jaccard_table <- add_jaccard_results(intersect_10, union_10, jaccard_table, "Total Shared Variants > 0.10")
  
  # Calculate jaccard index for total shared variants between 0.1 and 0.90
  intersect_1_90 <- filter(df, s1_af < 0.90 & s2_af < 0.90 & s1_af > 0.1 & s2_af > 0.1)
  union_1_90 <- filter(df, (s1_af < 0.90 & s1_af > 0.1) | (s2_af < 0.90 & s2_af > 0.1))
  jaccard_table <- add_jaccard_results(intersect_1_90, union_1_90, jaccard_table, "Total Shared Variants Between 0.10 and 0.90")

  return(jaccard_table)
  }

# Calculate jaccard index for mtDNA-Server (separately for SNPs and indels, and both for just Mutect)
results_table <- generate_jaccard_table("dups", "MtDNA-Server", "snps", snps_only)
write.table(results_table, file = "dups_snps_mtdnaserver.txt", sep = "\t", quote = FALSE, row.names = FALSE)
dup_jaccard_results = results_table

results_table <- generate_jaccard_table("dups", "MtDNA-Server", "indels", indels_only)
write.table(results_table, file = "dups_indels_mtdnaserver.txt", sep = "\t", quote = FALSE, row.names = FALSE)
dup_jaccard_results = left_join(dup_jaccard_results, results_table, by = "Metric")

# Calculate jaccard index for Mutect (separately for SNPs and indels)
results_table <- generate_jaccard_table("dups", "Mutect", "snps", snps_only)
write.table(results_table, file = "dups_snps_mutect.txt", sep = "\t", quote = FALSE, row.names = FALSE)
dup_jaccard_results = left_join(dup_jaccard_results, results_table, by = "Metric")

results_table <- generate_jaccard_table("dups", "Mutect", "indels", indels_only)
write.table(results_table, file = "dups_indels_mutect.txt", sep = "\t", quote = FALSE, row.names = FALSE)
dup_jaccard_results = left_join(dup_jaccard_results, results_table, by = "Metric")

results_table <- generate_jaccard_table("dups", "Mutect", "both", dups)
write.table(results_table, file = "dups_snps_and_indels_mutect.txt", sep = "\t", quote = FALSE, row.names = FALSE)
dup_jaccard_results = left_join(dup_jaccard_results, results_table, by = "Metric")

# Output final jaccard result table containing the results from all of the comparisons
write.table(dup_jaccard_results, file = "dup_jaccard_results.txt", sep = "\t", quote = FALSE, row.names = FALSE)
