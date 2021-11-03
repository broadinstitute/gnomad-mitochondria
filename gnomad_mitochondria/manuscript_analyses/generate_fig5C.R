library(ggplot2)
library(tidyr)

dir.create("figures")

synthetic_vep <- read.delim(file = 'final_data_files/synthetic_vcf/NC_012920.1_synthetic_vep_splitvarstwogenes.vcf', header = TRUE, sep = "\t", dec = ".") 
gnomad <- read.delim(file = 'reformated.vcf', header = TRUE, sep = "\t", na.strings = "")

# synonymous transitions - format data frame and count number of expected
Ts_list <- c("A>G", "C>T", "G>A", "T>C")
synthetic_vep$count <- 1
synthetic_vep_syn <- synthetic_vep[synthetic_vep$Consequence == "synonymous_variant", ]
synthetic_vep_syn$mut <- ifelse(paste(synthetic_vep_syn$REF, ">", synthetic_vep_syn$ALT, sep = "") %in% Ts_list, paste(synthetic_vep_syn$REF, ">", synthetic_vep_syn$ALT, sep = ""), "Tv")
expected_syn <- aggregate(synthetic_vep_syn$count, by = list(Category = synthetic_vep_syn$mut), FUN = sum)
colnames(expected_syn)[c(1, 2)] <- c("Mutation", "Expected")

# synonymous transitions - format data frame and count observed
gnomad$count <- 1
gnomad$mut <- ifelse(paste(gnomad$REF, ">", gnomad$ALT, sep = "") %in% Ts_list, paste(gnomad$REF, ">", gnomad$ALT, sep = ""), "Tv")
# heteroplasmic
gnomad_SNVs_het_syn <- gnomad[gnomad$Consequence == "synonymous_variant" & gnomad$VARIANT_CLASS == "SNV" & gnomad$max_hl < 0.95, ]
obs_het_syn <- aggregate(gnomad_SNVs_het_syn$count, by = list(Category2 = gnomad_SNVs_het_syn$mut), FUN = sum)
colnames(obs_het_syn)[c(2)] <- "Observed_het"
# homoplasmic
gnomad_SNVs_hom_syn <- gnomad[gnomad$Consequence == "synonymous_variant" & gnomad$VARIANT_CLASS == "SNV" & gnomad$max_hl >= 0.95, ]
obs_hom_syn <- aggregate(gnomad_SNVs_hom_syn$count, by = list(Category3 = gnomad_SNVs_hom_syn$mut), FUN = sum)
colnames(obs_hom_syn)[c(2)] <- "Observed_hom"

# bind and calculate proportion possible observed ratios
table_syn <- cbind(obs_hom_syn, obs_het_syn, expected_syn)
table_syn$ratio_het <- table_syn$Observed_het / table_syn$Expected
table_syn$ratio_hom <- table_syn$Observed_hom / table_syn$Expected
table_syn$ratio_none <- (1 - table_syn$ratio_hom - table_syn$ratio_het)

# reformat for plotting
percent_poss_vs_obs_syn <- table_syn[, c("Mutation", "ratio_het", "ratio_hom", "ratio_none")] %>% gather(Data, Ratio, ratio_het, ratio_hom, ratio_none)
percent_poss_vs_obs_syn$Data <- factor(percent_poss_vs_obs_syn$Data, c("ratio_none", "ratio_het", "ratio_hom"))

plotc <- ggplot(data = percent_poss_vs_obs_syn, aes(x = Mutation, y = Ratio, fill = Data)) +
  geom_bar(stat = "identity", width = 0.75, position = "fill", colour = "black") + 
  labs(x = "Synonymous variants", y = "Proportion possible SNVs observed") +
  scale_y_continuous(breaks = seq(0, 1, 0.1), expand = c(0.01, 0.01)) + 
  scale_fill_manual(values = c("White", "Grey", "Black"), labels = c("Not observed", "Seen at heteroplasmy only", "Seen at homoplasmy")) +
  theme(axis.title.x = element_text(size = 20, face = "bold"), 
        axis.text.x  = element_text(size = 20, colour = "black"),
        axis.title.y = element_text(size = 20, vjust = 0.5, face = "bold"), 
        axis.text.y  = element_text(size = 20, colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
  guides(fill = FALSE)

save(plotc, file = "figures/Fig5c.rdata")
