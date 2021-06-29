library(dplyr)
library(tidyr)
library(ggplot2)

dir.create("figures")

gnomad_vep_plasmy <- read.delim(file='reformated.vcf', header=TRUE, sep = "\t", na.strings = "")

#select variants and rename for plotting 
gnomad_vep_plasmy$Consequence_clean <- ifelse(gnomad_vep_plasmy$Consequence=="stop_gained",print("Stop gain\nSNVs"),
                                              ifelse(gnomad_vep_plasmy$VARIANT_CLASS=="SNV" & gnomad_vep_plasmy$BIOTYPE=="Mt_tRNA" & !is.na(gnomad_vep_plasmy$BIOTYPE),print("tRNA\nSNVs"),
                                                     ifelse(gnomad_vep_plasmy$VARIANT_CLASS=="SNV" & gnomad_vep_plasmy$BIOTYPE=="Mt_rRNA" & !is.na(gnomad_vep_plasmy$BIOTYPE),print("rRNA\nSNVs"),
                                                            ifelse(gnomad_vep_plasmy$Consequence=="missense_variant",print("Missense\nSNVs"),
                                                                   ifelse(gnomad_vep_plasmy$Consequence=="synonymous_variant",print("Synonymous\nSNVs"),
                                                                          ifelse(gnomad_vep_plasmy$VARIANT_CLASS!="SNV" & gnomad_vep_plasmy$BIOTYPE=="protein_coding" & !is.na(gnomad_vep_plasmy$BIOTYPE),print("Protein\nindels"),
                                                                                 ifelse(gnomad_vep_plasmy$VARIANT_CLASS!="SNV" & gnomad_vep_plasmy$BIOTYPE=="Mt_tRNA" & !is.na(gnomad_vep_plasmy$BIOTYPE),print("tRNA\nindels"),
                                                                                        ifelse(gnomad_vep_plasmy$VARIANT_CLASS!="SNV" & gnomad_vep_plasmy$BIOTYPE=="Mt_rRNA" & !is.na(gnomad_vep_plasmy$BIOTYPE),print("rRNA\nindels"),
                                                                                               ifelse(gnomad_vep_plasmy$Consequence=="intergenic_variant" & (gnomad_vep_plasmy$POS<577 | gnomad_vep_plasmy$POS>16023),print("Control reg.\nSNVs & indels"),
                                                                                                      print("other"))))))))))
#Only heteroplasmic
gnomad_vep_plasmy_het <- gnomad_vep_plasmy[gnomad_vep_plasmy$max_hl<0.95 & gnomad_vep_plasmy$Consequence_clean!="other",]
gnomad_vep_plasmy_het$Consequence_clean <- factor(gnomad_vep_plasmy_het$Consequence_clean, levels=c("Control reg.\nSNVs & indels","Synonymous\nSNVs","Missense\nSNVs","Stop gain\nSNVs","rRNA\nSNVs","tRNA\nSNVs", "Protein\nindels","rRNA\nindels","tRNA\nindels"))

plotb <- ggplot(data=gnomad_vep_plasmy_het, aes(x=Consequence_clean, y=as.numeric(max_hl))) +
  geom_boxplot(fill="Grey",outlier.size=1) + 
  labs(y="Maximum heteroplasmy level") + 
  scale_y_continuous(breaks=seq(0,1.0,0.1), expand =c(0.01,0.01), labels = scales::percent_format(accuracy = 1)) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_text(size=18),
        axis.title.y = element_text(size=20, vjust=0.5, face="bold"),
        axis.text.y  = element_text(size=20),
        axis.ticks.length = unit(0.15, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        plot.margin = unit(c(0.5,0.5,0.5,0.1), "cm"))

save(plotb,file = "figures/Fig5b.rdata")
