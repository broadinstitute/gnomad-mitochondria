library(dplyr)
library(forcats)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(scales)
library(tidyr)

dir.create("figures")

gnomad <- read.delim(file = 'reformated.vcf', header = TRUE, sep = "\t", na.strings = "")

# adding in the one cfrm indel variant manually inspected to be in gnomAD, and not recognized due to left vs right alignment
# listed as m.7471C>CC in mitomap, equivalent to m.7465A>AC in gnomad - note in literature also called 7472insC, m.7471dupC
gnomad_cfrm <- data.frame(gnomad[(gnomad$Mitomap_dz_status == "Cfrm" & !is.na(gnomad$Mitomap_dz_status)) | (gnomad$POS == 7465 & gnomad$ALT == "AC"), ] %>% 
                            mutate(all_hl = strsplit(as.character(all_hl), ",")) %>% unnest(all_hl))
gnomad_cfrm$all_hl <- gsub("\\[|]|'|'", "", gnomad_cfrm$all_hl)
gnomad_cfrm$variant <- paste("m.", gnomad_cfrm$POS, gnomad_cfrm$REF, ">", gnomad_cfrm$ALT, sep = "")
gnomad_cfrm$count <- gnomad_cfrm$AC_hom + gnomad_cfrm$AC_het
gnomad_cfrm$variant <- factor(gnomad_cfrm$variant, levels = c(unique(gnomad_cfrm[order(gnomad_cfrm$count, -gnomad_cfrm$POS), c("variant")]))) #order by count for plot

# total carrier frequency - write to file
write.table(sub("^", "Total carrier frequency of pathogenic variants in gnomAD (>10% heteroplasmy): 1 in ", comma(ceiling(1 / ((sum(gnomad_cfrm[!duplicated(gnomad_cfrm$variant), "count"])) / 56434)), accuracy = 1)),
            file = "figures/total_pathogenic_carrier_freq.txt", append = FALSE, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# dotplot of heteroplasmy level
plot1 <- ggplot(data = gnomad_cfrm, aes(y = variant, x = as.numeric(all_hl))) + 
  geom_jitter(cex = 2.5, position = position_jitter(h = 0.1)) +
  labs(x = "Heteroplasmy level") +
  scale_x_continuous(breaks = seq(0, 1.05, 0.1), expand = c(0.01, 0.01), labels = scales::percent_format(accuracy = 1)) + 
  theme(axis.title.x = element_text(size = 13, face = "bold"), 
        axis.text.x  = element_text(size = 14), 
        axis.title.y = element_blank(),
        axis.text.y  = element_text(size = 13),
        axis.ticks.length = unit(0.15, "cm"),
        panel.grid.major = element_line(colour = "light grey", size = 0.2), 
        panel.grid.minor = element_line(colour = "light grey", size = 0.2),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.margin = unit(c(0.1, 0.3, 0.1, 0.15), "cm")) + 
  geom_vline(xintercept = c(0.95), linetype = "solid", color = "dark grey")

# barplot colored by haplogroups #regenerate the cfrm table, this time unnest by haplogroups
gnomad_cfrm <- data.frame(gnomad[(gnomad$Mitomap_dz_status == "Cfrm" & !is.na(gnomad$Mitomap_dz_status)) | (gnomad$POS == 7465 & gnomad$ALT == "AC"), ] %>% 
                            mutate(all_haplogroups = strsplit(as.character(all_haplogroups), ",")) %>% unnest(all_haplogroups))
gnomad_cfrm$all_haplogroups <- gsub("\\[|]|'|'| |", "", gnomad_cfrm$all_haplogroups)
gnomad_cfrm$variant <- paste("m.", gnomad_cfrm$POS, gnomad_cfrm$REF, ">", gnomad_cfrm$ALT, sep = "")
gnomad_cfrm$count <- gnomad_cfrm$AC_hom + gnomad_cfrm$AC_het
gnomad_cfrm$variant <- factor(gnomad_cfrm$variant, levels = c(unique(gnomad_cfrm[order(gnomad_cfrm$count, -gnomad_cfrm$POS), c("variant")]))) #order by count for plot
gnomad_cfrm$all_haplogroups <- factor(gnomad_cfrm$all_haplogroups, levels = c("L0", "L1", "L5", "L2", "L6", "L4", "L3", "M", "C", "Z", "E", "G", "Q", "D", "N", "Y", "A", "O", "S", "F", "B", "P", "I", "W", "X", "R", "HV", "V", "H", "J", "T", "U", "K"))
# to match subset haplogroups displayed on plot
haplo_colors = c("#e20000", "#f05032", "#fc7e61", #reds - african haplogroups
                 "#00910e", "#3f9f35", "#63af55", "#82be74", "#9fcd93", "#badbb0", "#d2e8cc", "#e7f2e3", #greens - asian haplogroups
                 "#00508e", "#275a95", "#3c659d", "#4e70a4", "#5c79ab", "#6983b1", "#788fb9", "#879ac0", "#96a6c8", "#a7b4d1") #blues - european haplogroups
# other displays on plot
gnomad_cfrm$cf <- sub("^", "1 in ",comma(ceiling(1 / ((gnomad_cfrm$count) / gnomad_cfrm$AN)), accuracy = 1))
gnomad_cfrm$Mitomap_dz_homoplasmy <- ifelse(gnomad_cfrm$variant == "m.7465A>AC", "+", as.character(gnomad_cfrm$Mitomap_dz_homoplasmy)) #manually adding this datapoint for misaligned indel in 
# clean up associated disease phenotypes listed in MITOMAP for plot
gnomad_cfrm$Mitomap_disease <- ifelse(gnomad_cfrm$variant == "m.7465A>AC", "AMDF, other", as.character(gnomad_cfrm$Mitomap_disease)) #manually adding this in for misaligned indel in #original is "PEM / AMDF / Motor neuron disease-like"
gnomad_cfrm$Mitomap_disease <- gsub("Leigh Syndrome|Leigh Disease|Leigh syndrome|- LD", "LS", gnomad_cfrm$Mitomap_disease)
gnomad_cfrm$Mitomap_disease <- gsub("NARP-like disease", "NARP", gnomad_cfrm$Mitomap_disease)
gnomad_cfrm$Mitomap_disease <- gsub("Mitochondrial myopathy, lactic acidosis and sideroblastic anemia \\(MLASA)", "MLASA", gnomad_cfrm$Mitomap_disease)
gnomad_cfrm$Mitomap_disease <- gsub("\\/Deafness",", DEAF", gnomad_cfrm$Mitomap_disease)
gnomad_cfrm$Mitomap_disease <- ifelse(grepl("Multiorgan failure / |\\ / HCM|LDYT / |MICM\\+DEAF / |\\ / MILS|Encephalopathy / | Autism /| Other|\\ / Depressive mood disorder / leukoencephalopathy / HiCM|; autism spectrum intellectual disability; possibly antiatherosclerotic|\\ / Ataxia\\+Lipomas|\\ / carotid atherosclerosis risk|\\ / Progressive Dystonia|FBSN / |BSN / |\\ / other|\\ / FSGS / ASD / Cardiac\\+multi-organ dysfunction|\\ / dystonia|\\/ DMDF / MIDD |\\ / CPEO|\\ / IgG nephropathy", gnomad_cfrm$Mitomap_disease),
                                      paste(gsub("Multiorgan failure / |\\ / HCM|LDYT / |MICM\\+DEAF / |\\ / MILS|Encephalopathy / | Autism /| Other|\\ / Depressive mood disorder / leukoencephalopathy / HiCM|; autism spectrum intellectual disability; possibly antiatherosclerotic|\\ / Ataxia\\+Lipomas|\\ / carotid atherosclerosis risk|\\ / Progressive Dystonia|FBSN / |BSN / |\\ / other|\\ / FSGS / ASD / Cardiac\\+multi-organ dysfunction|\\ / dystonia|\\/ DMDF / MIDD |\\ / CPEO|\\ / IgG nephropathy", "", gnomad_cfrm$Mitomap_disease), ", other", sep = ""),as.character(gnomad_cfrm$Mitomap_disease))
gnomad_cfrm$Mitomap_disease <- gsub("\\ /|;", ",", gnomad_cfrm$Mitomap_disease)
gnomad_cfrm$Mitomap_disease <- gsub("myopathy,", "MM,", gnomad_cfrm$Mitomap_disease)

plot2 <- ggplot(data = gnomad_cfrm, aes(variant, fill = all_haplogroups)) + 
  geom_bar(stat = "Count") + coord_flip(clip = "off") + 
  labs(y = "Number of carriers") +
  scale_y_continuous(breaks = c(1, 10, 20, 30, 40, 50, 60, 70), expand = c(0.01, 0.01)) +
  scale_fill_manual(values = haplo_colors) + 
  theme(axis.title.x = element_text(size = 13, face = "bold"),
        axis.text.x  = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.length = unit(0.15, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position=c(.525, .175),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        plot.margin = unit(c(0.1, 13, 0.1, 0.15), "cm"),
        plot.tag.position = c(1.67, 0.015)) + 
  geom_text(aes(label = cf, hjust = 0), y = 76, colour = "#595757", size = 4.75) + 
  geom_text(aes(label = Mitomap_disease, hjust = 0), y = 100, colour = "#595757", size = 4.75) + 
  labs(tag = expression(bold("Carrier freq.    Associated diseases             Hom. reported"))) +
  geom_text(aes(label = Mitomap_dz_homoplasmy, hjust = 0), y = 168, colour = "#595757", size = 4.75) +
  guides(fill = guide_legend(title = "Haplogroup", ncol = 3))

p1.grob <- ggplotGrob(plot1)
p2.grob <- ggplotGrob(plot2)
maxheights <- grid::unit.pmax(p1.grob$heights[2:5], p2.grob$heights[2:5])
p1.grob$heights[2:5] <- as.list(maxheights)
p2.grob$heights[2:5] <- as.list(maxheights)

#options(bitmapType = 'cairo', device = 'png')
ggsave("figures/Figure6.png", grid.arrange(p1.grob, p2.grob, widths=c(2, 2.5)), width = 16, height = 8) 

dev.off()

