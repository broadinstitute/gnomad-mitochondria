library(ggbeeswarm)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(scales)
library(tidyr)

dir.create("figures")

# plot A
mitomap_disease <- read.delim(file = 'final_data_files/other_databases/MITOMAP_disease_02012021.cgi', header = TRUE, sep = "\t", na.strings = "")

cfrm <- mitomap_disease[mitomap_disease$status == "Cfrm", ]
cfrm$variant <- paste("m.", cfrm$pos, cfrm$ref, ">", cfrm$alt, sep = "")

# manually inspected and updated cfrm indels to expected variant call in gnomad (left vs right alignment calls)
# m.3271T>: #equivalent to m.3270CT>C
# m.3902ACCTTGC>GCAAGGT #equivalent to m.3901GACCTTGC>GGCAAGGT ? #inversion
# m.5537A>AT #same
# m.7471C>CC #equivalent to m.7465A>AC #note in literature also called 7472insC, m.7471dupC
# m.9205TA>: #equivalent to m.9204ATA>A

cfrm$variant <- ifelse(cfrm$variant == "m.3271T>:", print("m.3270CT>C"),
                       ifelse(cfrm$variant == "m.7471C>CC", print("m.7465A>AC"),
                              ifelse(cfrm$variant == "m.9205TA>:", print("m.9204ATA>A"), as.character(cfrm$variant))))

gnomad <- read.delim(file = 'reformated.vcf', header = TRUE, sep = "\t", na.strings = "")
gnomad$variant <- paste("m.", gnomad$POS, gnomad$REF, ">", gnomad$ALT, sep = "")

# annotate confirmed disease variants with the maximum heteroplasmy observed in gnomad
merged_table <- merge(cfrm, gnomad[, c("variant", "max_hl")], by = 'variant', all = TRUE)
cfrm <- merged_table[merged_table$status == "Cfrm" & !is.na(merged_table$status) & !duplicated(merged_table$variant), ]
cfrm$dz_plasmy <- ifelse(cfrm$heteroplasmy == "+" & cfrm$homoplasmy == "-", "het_only", "at_hom")

# count number of observed, heteroplasmic separately to homoplasmic
cfrm$count <- 1 #hack to allow quick count

# num not observed
cfrm_none <- cfrm[is.na(cfrm$max_hl), ]
not_obs <- aggregate(cfrm_none$count, by = list(Category = cfrm_none$dz_plasmy), FUN = sum)
colnames(not_obs)[c(2)] <- "Not_observed"

# count number of observed, SNVs only, at HETEROPLASMY only
cfrm_het <- cfrm[cfrm$max_hl < 0.95 & !is.na(cfrm$max_hl), ]
obs_het <- aggregate(cfrm_het$count, by = list(Category1 = cfrm_het$dz_plasmy), FUN = sum)
colnames(obs_het)[c(2)] <- "Observed_het"

# count number of observed, SNVs only, at HOMOPLASMY (>95% heteroplasmy)
cfrm_hom <- cfrm[cfrm$max_hl >= 0.95 & !is.na(cfrm$max_hl), ]
obs_hom <- aggregate(cfrm_hom$count, by = list(Category2 = cfrm_hom$dz_plasmy), FUN = sum)
colnames(obs_hom)[c(2)] <- "Observed_hom"

# bind and calculate ratios
table <- cbind(not_obs, obs_hom, obs_het) 
plot_table  <- table[, c("Category", "Not_observed", "Observed_het", "Observed_hom")] %>% gather(Data, Number, Not_observed, Observed_het, Observed_hom)
plot_table$Category <- factor(plot_table$Category, levels = c("het_only", "at_hom"))

plota <- ggplot(data = plot_table, aes(x = Category, y = as.numeric(Number), fill = Data)) +
  geom_bar(stat = "identity", position = "fill", width = 0.75, colour = "black") +
  labs(y = "Proportion observed", x = "MITOMAP disease association") +
  scale_y_continuous(breaks = seq(0, 1, 0.1), expand = c(0.01, 0.01)) +
  scale_fill_manual(values = c("White", "Grey", "Black"), labels = c("Not observed", "At heteroplasmy only ", "At homoplasmy")) +
  scale_x_discrete(labels = c("Heteroplasmy only\nn=56", "Homoplasmy\nn=38")) +
  theme(axis.title.x = element_text(size = 14, face = "bold"),
        axis.text.x  = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 14, vjust = 0.5, face = "bold"),
        axis.text.y  = element_text(size = 14, colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14), 
        plot.margin = unit(c(0.5, 10, 0.5, 10), "cm")) +
  labs(fill = "gnomAD:")

# plot B
# adding in the one cfrm indel variant manually inspected to be in gnomAD, and not recognized due to left vs right alignment
# listed as m.7471C>CC in mitomap #equivalent to m.7465A>AC #note in literature also called 7472insC, m.7471dupC
gnomad_cfrm_long <- gnomad[gnomad$Mitomap_dz_status == "Cfrm" & !is.na(gnomad$Mitomap_dz_status) | (gnomad$POS == 7465 & gnomad$ALT == "AC"), c("AF_hom", "AF_het", "variant", "max_hap_AF_hom", "max_pop_AF_hom", "Mitomap_dz_homoplasmy")] %>% 
  gather(category, AF, AF_hom, AF_het, max_hap_AF_hom, max_pop_AF_hom)
gnomad_cfrm_long$Mitomap_dz_homoplasmy <- ifelse(gnomad_cfrm_long$variant == "m.7465A>AC", "+", as.character(gnomad_cfrm_long$Mitomap_dz_homoplasmy)) #manually adding this for the misaligned indel
gnomad_cfrm_long$category <- factor(gnomad_cfrm_long$category, levels = c("AF_het", "AF_hom", "max_hap_AF_hom", "max_pop_AF_hom"))

options(scipen = 999) 

mysqrt_trans <- function() { #this is a hack to enable y=0 tick mark on plot per https://gist.github.com/DarwinAwardWinner/21652acf017880c271e95cc2e35574f4
  trans_new("mysqrt", 
            transform = base::sqrt,
            inverse = function(x) ifelse(x < 0, 0, x ^ 2),
            domain = c(0, Inf))
}

plotb <- ggplot(data = gnomad_cfrm_long, aes(y = AF, x = category, color = Mitomap_dz_homoplasmy)) + 
  geom_beeswarm(size = 1) +
  labs(y = "Allele frequency") +
  scale_y_continuous(trans = "mysqrt", breaks = c(0, 1/50000, 5/50000, 0.001, 0.002, 0.005), expand = c(0.002, 0.002)) +
  scale_x_discrete(labels = c("AF het", "AF hom", "Max. haplogroup-specific\nAF hom", "Max. population-specific\nAF hom")) +
  theme(axis.title.y = element_text(size = 14, face = "bold"), 
        axis.text.y  = element_text(size = 14), 
        axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 14),
        axis.ticks.length = unit(0.15, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
        axis.line = element_line(colour = "black"),
        plot.tag.position = c(0.875, 0.56),
        plot.tag = element_text(color = "dark grey", size = 14),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 13)) + 
  labs(color = "MITOMAP disease\nassociation") +
  scale_color_manual(labels = c("Heteroplasmy only", "Homoplasmy"), values = c("Orange", "Navy")) + 
  geom_hline(yintercept = c(2e-05), linetype = "dashed", color = "dark grey") + 
  geom_hline(yintercept = c(0.005), linetype = "dashed", color = "dark grey") +
  labs(tag = "BS1\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n                  PM2_supporting") +
  geom_text_repel(aes(label = ifelse(AF > 0.00002 & category != "AF_het", as.character(variant), '')), force = 0.5, nudge_x = 0.15, direction = "y", hjust = 0, segment.size = 0.2, show.legend = FALSE) +
  geom_text_repel(aes(label = ifelse(AF > 0.00002 & category == "AF_het", as.character(variant), '')), force = 0.5, nudge_x= -0.25, direction = "y" , hjust = 1, segment.size = 0.2, show.legend = FALSE)

ggarrange(plota, plotb, ncol = 1, nrow = 2, heights = c(5.5, 10), labels = c("A", "B"))

#options(bitmapType = 'cairo', device = 'png')
ggsave("figures/FigureS7.png", width = 15, height = 13)

dev.off()

#% cfrm variants with AF_hom <0.00002
((94 - nrow(gnomad_cfrm_long[gnomad_cfrm_long$category == "AF_hom" & gnomad_cfrm_long$AF > 0.00002, ])) / 94)
#% cfrm variants with AF_het <0.00002
((94 - nrow(gnomad_cfrm_long[gnomad_cfrm_long$category == "AF_het" & gnomad_cfrm_long$AF > 0.00002, ])) / 94)
