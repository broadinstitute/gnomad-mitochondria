library(ggplot2)
library(tidyr)

dir.create("figures")

synthetic_vep <- read.delim(file = 'final_data_files/synthetic_vcf/NC_012920.1_synthetic_vep_splitvarstwogenes.vcf', header = TRUE, sep = "\t", dec = ".")
gnomad <- read.delim(file = 'reformated.vcf', header = TRUE, sep = "\t", na.strings = "")

# count number of expected, first remove m.3107N
synthetic_vep <- synthetic_vep[synthetic_vep$REF != "N", ]
synthetic_vep$count <- 1
synthetic_vep$Consequence <- ifelse(synthetic_vep$BIOTYPE == "Mt_rRNA","rRNA", ifelse(synthetic_vep$BIOTYPE == "Mt_tRNA", "tRNA", as.character(synthetic_vep$Consequence)))
# select variant types for display
synthetic_vep <- synthetic_vep[synthetic_vep$Consequence == "synonymous_variant" | synthetic_vep$Consequence == "missense_variant" | synthetic_vep$Consequence == "stop_gained" | synthetic_vep$Consequence == "tRNA" | synthetic_vep$Consequence == "rRNA" | synthetic_vep$POS < 577 | synthetic_vep$POS > 16023, ]
expected <- aggregate(synthetic_vep$count, by = list(Category = synthetic_vep$Consequence), FUN = sum)
colnames(expected)[c(2)] <- c("Expected")

# count number of observed, heteroplasmic separately to homoplasmic
gnomad$count <- 1
gnomad$Consequence <- ifelse(gnomad$BIOTYPE == "Mt_rRNA" & !is.na(gnomad$BIOTYPE), "rRNA", ifelse(gnomad$BIOTYPE == "Mt_tRNA" & !is.na(gnomad$BIOTYPE), "tRNA", as.character(gnomad$Consequence)))
gnomad <- gnomad[gnomad$Consequence == "synonymous_variant" | gnomad$Consequence == "missense_variant" | gnomad$Consequence == "stop_gained" | gnomad$Consequence == "tRNA" | gnomad$Consequence == "rRNA" | gnomad$POS < 577 | gnomad$POS > 16023, ]

# count number of observed, SNVs only, at HETEROPLASMY only
gnomad_SNVs_het <- gnomad[gnomad$VARIANT_CLASS == "SNV" & gnomad$max_hl < 0.95, ]
obs_het <- aggregate(gnomad_SNVs_het$count, by = list(Category2 = gnomad_SNVs_het$Consequence), FUN = sum)
colnames(obs_het)[c(2)] <- "Observed_het"

# count number of observed, SNVs only, at HOMOPLASMY (>95% heteroplasmy level)
gnomad_SNVs_hom <- gnomad[gnomad$VARIANT_CLASS == "SNV" & gnomad$max_hl >= 0.95, ]
obs_hom <- aggregate(gnomad_SNVs_hom$count, by = list(Category3 = gnomad_SNVs_hom$Consequence), FUN = sum)
colnames(obs_hom)[c(2)] <- "Observed_hom"

# bind and calculate proportion possible ratios
table <- cbind(obs_hom, obs_het, expected)
table$ratio_het <- table$Observed_het / table$Expected
table$ratio_hom <- table$Observed_hom / table$Expected
table$ratio_none <- (1 - table$ratio_hom - table$ratio_het)

# reformat for plotting
plot_table <- table[, c("Category", "ratio_het", "ratio_hom", "ratio_none")] %>% gather(Data, Ratio, ratio_het, ratio_hom, ratio_none)
plot_table$Category <- factor(plot_table$Category, c("intergenic_variant", "synonymous_variant", "missense_variant", "stop_gained", "rRNA", "tRNA"))
plot_table$Data <- factor(plot_table$Data, c("ratio_none", "ratio_het", "ratio_hom"))

plota <- ggplot(data=plot_table, aes(x = Category, y = Ratio, fill = Data)) +
  geom_bar(stat = "identity", width = 0.75, position = "fill", colour = "black") +
  labs(y = "Proportion possible SNVs observed") +
  scale_y_continuous(breaks = seq(0, 1, 0.1), expand = c(0.01, 0.01)) +
  scale_fill_manual(values = c("White", "Grey", "Black"), labels = c("Not observed", "At heteroplasmy only", "At homoplasmy")) +
  scale_x_discrete(labels = c("Control reg. ", "Synonymous", "Missense", "Stop gain", "rRNA", "tRNA")) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 20, vjust = 0.5, face = "bold"),
        axis.text.y  = element_text(size = 20, colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        legend.text = element_text(size = 20),
        legend.title = element_blank(), 
        plot.margin = unit(c(0.5, 0.1, 0.5, 0.5), "cm"),
        legend.position = "top")  
 
save(plota, file = "figures/Fig5a.rdata")