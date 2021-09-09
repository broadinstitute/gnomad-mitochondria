#!/usr/bin/env Rscript

library(ggplot2)
library(optparse)
library(plyr)

option_list = list(
  make_option(c("-v", "--vl_file"), type = "character", default = NULL,
              help = "The fine_snp.txt file containing tab-delimited mix, mix_level, and vl", metavar = "character"),
  make_option(c("-d", "--plot_directory"), type = "character", default = NULL,
              help = "Directory to which plots should be written", metavar = "character"))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

pdf(NULL)

data <- read.table(opt$vl_file, header = TRUE, colClasses = c("character", "character", "numeric"))

data$mix_level <- revalue(data$mix_level, c("9" = "90", "5" = "50"))

# Plot distribution of detected heteroplasmy levels for each mixin (expected heteroplasmy level)
heteroplasmy_hist <- ggplot(data, aes(x = vl)) + 
  geom_histogram() +
  facet_grid(. ~ mix_level) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.y = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 14, angle = 90, colour = "black"),
        axis.title = element_text(colour = "black", size = 18, face = "bold"),
        axis.text = element_text(colour = "black", size = 14),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  labs(x = "Heteroplasmy Level (%)", y = "Count")
heteroplasmy_hist

setwd(opt$plot_directory)
ggsave(heteroplasmy_hist, filename = "heteroplasmy_hist_mixin_snps.png", dpi = 300, width = 10, height = 6, units = "in")
