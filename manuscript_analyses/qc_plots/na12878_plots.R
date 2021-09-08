#!/usr/bin/env Rscript

library(dplyr)
library(ggplot2)
library(optparse)
library(stringr)
library(tidyr)

option_list = list(
  make_option(c("-a", "--mtdnaserver_vl_output"), type = "character", default = NULL,
              help = "The vl_output.txt file for mtDNA-Server results", metavar = "character"),
  make_option(c("-b", "--mutect_vl_output"), type = "character", default = NULL,
              help = "The vl_output.txt file for Mutect results", metavar = "character"),
  make_option(c("-d", "--plot_directory"), type = "character", default = NULL,
  help="Directory to which plots should be written", metavar = "character"),
  make_option(c("-t", "--min_het_threshold"), type = "numeric", default = NULL,
              help = "Minimum heteroplasmy level to define a variant as a PASS heteroplasmic variant, error out if below", metavar = "float"))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# Read in and rename columns for mtDNA-Server results
mtdnaserver <- read.table(opt$mtdnaserver_vl_output, sep = '\t', quote = "")
mtdnaserver <- gather(mtdnaserver, "sample", "value", V7:V19)
mtdnaserver <- mtdnaserver %>% rename(chrom = V1, pos = V2, ref = V3, alt = V4, index = V5, filter = V6)

# Read in and rename columns for Mutect results
mutect <- read.table(opt$mutect_vl_output, sep = '\t', quote = "")
mutect <- gather(mutect, "sample", "value", V7:V19)
mutect <- mutect %>% rename(chrom = V1, pos = V2, ref = V3, alt = V4, index = V5, filter = V6)

# Combine mtDNA-Server and Mutect dataframes, filter to PASS variants (non artifact-prone sites)
mtdnaserver <- mutate(mtdnaserver, caller = "MtDNA-Server")
mutect <- mutate(mutect, caller = "Mutect")
data <- rbind(mtdnaserver, mutect)
data <- filter(data, filter == "PASS")

# Split value column into sample and heteroplasmy level columns
data$sample <- stringr::str_split_fixed(data$value, "=", 2)[ , 1]
data$hl <- stringr::str_split_fixed(data$value, "=", 2)[ , 2]
data <- mutate(data, hl = as.numeric(hl))
data <- mutate(data, variant_id = paste(chrom, pos, ref, alt, sep = "_"))

# Double check that there are no non-zero heteroplasmy levels below the expected heteroplasmy level set by min_het_threshold 
data <- mutate(data, lower_than_hl_threshold=ifelse((hl > 0) & (hl < opt$min_het_threshold), TRUE, FALSE))
if (max(data$lower_than_hl_threshold) == 1) {
  stop("Found heteroplasmy level below expectation based on min_het_threshold")
  }

# Make boxplot of heteroplasmy level for each variant found in at least one NA12878 sample split by SNP/indel
data <- mutate(data, variant_type = ifelse(nchar(as.character(ref)) == 1 & nchar(as.character(alt)) == 1, "SNP", "Indel"))
vaf_condordance_variant_type <- ggplot(data, aes(x = reorder(variant_id, pos), y = hl)) + 
  geom_boxplot(outlier.colour = "red") +
  facet_grid(variant_type ~ caller) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 14, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 10, color = "black", angle = 90, hjust = 1),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 18, face = "bold", color = "black")) +
  labs(x = "Variant",  y = "Heteroplasmy Level")
vaf_condordance_variant_type

setwd(opt$plot_directory)
ggsave(vaf_condordance_variant_type, filename = "na12878_vaf.png", dpi = 300, width = 10, height = 6, units = "in")

