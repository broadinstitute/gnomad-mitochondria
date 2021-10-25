library(dplyr)
library(ggplot2)
library(reshape2)
library(scales)
library(tidyr)

dir.create("figures")

gnomad <- read.delim(file = 'reformated.vcf', header = TRUE, sep = "\t", na.strings = "")

# this file has two rows for variants in two genes (where each gene consequence on separate line), remove duplicates for this analysis
gnomad$var <- paste(gnomad$REF, gnomad$POS, gnomad$ALT, sep = '')
gnomad <- gnomad[!duplicated(gnomad$var), ]

# data frame with heteroplasmy on each line
gnomad_all <- data.frame(gnomad %>% mutate(all_hl = strsplit(as.character(all_hl), ",")) %>% unnest(all_hl))
gnomad_all$all_hl <- gsub("\\[|]|'|'", "", gnomad_all$all_hl)

# cdf plots, for all variants, and max het for unique variants, overlay
gnomad_all$variable <- "gnomad_all.all_hl"
df <- data.frame(gnomad_all[, c("variable", "all_hl")]) 
colnames(df)[2] <- "value"
df2 <- data.frame(gnomad$max_hl)

# need to merge since different lengths, all=TRUE keeps all
data <- merge(df, melt(df2), all = TRUE)

# plot
ggplot(data, aes(x = as.numeric(value), colour = variable)) + stat_ecdf() +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0.0, 1.0, 0.1), labels = percent_format(accuracy = 2)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0.0, 1.0, 0.25)) + 
  labs(x = "Heteroplasmy level", y = "Empirical CDF") +
  theme(axis.title.x = element_text(size = 12, vjust = -0.75, face = "bold"),
        axis.text.x  = element_text(size = 12, angle = 0), 
        axis.title.y = element_text(size = 12, vjust = 0.5, face = "bold"),
        axis.text.y  = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        legend.text = element_text(size = 12)) +
  guides(colour = guide_legend(override.aes = list(size = 2))) + 
  labs(colour = "Data") +
  scale_color_hue(labels = c("\nHeteroplasmy of all \nvariant calls\n", "\nMaximum heteroplasmy \nof all unique variants\n")) +
  geom_vline(xintercept = 0.95, linetype = "dotted", color = "black", size = 0.5)

#options(bitmapType = 'cairo', device = 'png')
ggsave("figures/FigS5d.png", width = 10, height = 5)
ggsave("figures/FigS5d.pdf", width = 10, height = 5)
