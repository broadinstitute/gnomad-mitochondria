library(cowplot)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(plyr)
library(tidyr)

dir.create("figures")

synthetic_vep <- read.delim(file = 'final_data_files/synthetic_vcf/NC_012920.1_synthetic_vep_splitvarstwogenes.vcf', header = TRUE, sep = "\t", dec = ".") 
gnomad <- read.delim(file = 'reformated.vcf', header = TRUE, sep = "\t")

# possible - number of bases in RNA genes, or number of codons in each protein-coding gene, first remove m.3107N
synthetic_vep <- synthetic_vep[synthetic_vep$REF != "N", ]
genes_exp <- rbind(ddply(synthetic_vep[synthetic_vep$BIOTYPE == "Mt_tRNA" | synthetic_vep$BIOTYPE == "Mt_rRNA", ], ~SYMBOL, summarise, number_of_distinct_orders = length(unique(POS))),
                   ddply(synthetic_vep[synthetic_vep$BIOTYPE == "protein_coding", ], ~SYMBOL, summarise, number_of_distinct_orders = length(unique(Protein_position))))
colnames(genes_exp)[c(1, 2)] <- c("GENE_SYMBOL", "Expected")

# observed - how many bases in the RNA have a SNV in gnomad - or how many codons in the protein have a non-synonymous change
# number of bases/codons with only HETEROPLASMIC changes, note max_hl_nonsyn_SNV_codon only calculated from non-synonymous variants
genes_obs_het <- rbind(ddply(gnomad[gnomad$VARIANT_CLASS == "SNV" & (gnomad$BIOTYPE == "Mt_tRNA" | gnomad$BIOTYPE == "Mt_rRNA") & gnomad$max_hl_SNV_base < 0.95 & gnomad$max_hl_SNV_base > 0, ], ~SYMBOL, summarise, number_of_distinct_orders = length(unique(POS))),
                       ddply(gnomad[gnomad$VARIANT_CLASS == "SNV" & gnomad$BIOTYPE == "protein_coding" & gnomad$max_hl_nonsyn_SNV_codon < 0.95 & gnomad$max_hl_nonsyn_SNV_codon > 0, ], ~SYMBOL, summarise, number_of_distinct_orders = length(unique(Protein_position))))
colnames(genes_obs_het)[c(2)] <- "Observed_het"
# number of bases/codons with only HOMOPLASMIC changes
genes_obs_hom <- rbind(ddply(gnomad[gnomad$VARIANT_CLASS == "SNV" & (gnomad$BIOTYPE == "Mt_tRNA" | gnomad$BIOTYPE == "Mt_rRNA") & gnomad$max_hl_SNV_base >= 0.95, ], ~SYMBOL, summarise, number_of_distinct_orders = length(unique(POS))),
                       ddply(gnomad[gnomad$VARIANT_CLASS == "SNV" & gnomad$BIOTYPE == "protein_coding" & gnomad$max_hl_nonsyn_SNV_codon >= 0.95, ], ~SYMBOL, summarise, number_of_distinct_orders = length(unique(Protein_position))))
colnames(genes_obs_hom)[c(2)] <- "Observed_hom"

# bind and calculate proportion possible observed
constraint_bygene <- cbind(genes_exp, genes_obs_het, genes_obs_hom)
constraint_bygene$ratio_het <- constraint_bygene$Observed_het / constraint_bygene$Expected
constraint_bygene$ratio_hom <- constraint_bygene$Observed_hom / constraint_bygene$Expected
constraint_bygene$ratio_none <- (1 - constraint_bygene$ratio_het - constraint_bygene$ratio_hom)
constraint_bygene$BIOTYPE <- ifelse(grepl("MT-R", constraint_bygene$GENE_SYMBOL), print("Mt_rRNA"),
                                    ifelse(grepl("MT-T", constraint_bygene$GENE_SYMBOL), print("Mt_tRNA"),
                                           as.character("protein_coding")))
constraint_bygene <- constraint_bygene[, c("GENE_SYMBOL", "BIOTYPE", "ratio_het", "ratio_hom", "ratio_none")] %>% gather(Data, Ratio, ratio_het, ratio_hom, ratio_none)

# plot RNA and protein separately, ordered by ratio for no variants proportion
constraint_bygene_RNA <- constraint_bygene[constraint_bygene$BIOTYPE != "protein_coding", ]
list <- constraint_bygene_RNA[constraint_bygene_RNA$Data == "ratio_none", ]
list <- list[order(-list$Ratio), ]
constraint_bygene_RNA$GENE_SYMBOL <- factor(constraint_bygene_RNA$GENE_SYMBOL, c(as.character(list$GENE_SYMBOL)))
constraint_bygene_RNA$Data <- factor(constraint_bygene_RNA$Data, level = c("ratio_none", "ratio_het", "ratio_hom"))

# tRNA
plote <- ggplot(constraint_bygene_RNA[constraint_bygene_RNA$BIOTYPE == "Mt_tRNA", ],aes(x = GENE_SYMBOL, y = Ratio, fill = Data))+
  geom_bar(stat = "identity", width = 0.75, position = "fill", colour = "black") + 
  labs(y = "Proportion bases with SNVs") +
  scale_fill_manual(values = c("White", "Grey", "Black"), labels = c("Not observed", "Seen at heteroplasmy only", "Seen at homoplasmy")) +
  theme(axis.title.x = element_blank(), 
        axis.text.x  = element_text(size = 20, angle = 35, hjust = 1),
        axis.title.y = element_text(size = 20, vjust = 0.5, face = "bold"),
        axis.text.y  = element_text(size = 20),
        plot.margin = unit(c(0.5, 0, 0.9, 0.5), "cm"),
        axis.ticks.length = unit(0.15, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  guides(fill = FALSE) 

# rRNA
plotf <- ggplot(constraint_bygene_RNA[constraint_bygene_RNA$BIOTYPE == "Mt_rRNA", ],aes(x = GENE_SYMBOL, y = Ratio, fill = Data))+
  geom_bar(stat = "identity", width = 0.75, position = "fill", colour = "black") +
  scale_fill_manual(values = c("White", "Grey", "Black"), labels = c("Not observed", "Seen at heteroplasmy only", "Seen at homoplasmy")) +
  theme(axis.title.x = element_blank(), 
        axis.text.x  = element_text(size = 20, angle = 35, hjust = 1),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(size = 20),
        plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm"),
        axis.ticks.length = unit(0.15, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  guides(fill = FALSE)

# plot protein
constraint_bygene_prot <- constraint_bygene[constraint_bygene$BIOTYPE == "protein_coding", ]
list <- constraint_bygene_prot[constraint_bygene_prot$Data == "ratio_none", ]
list <- list[order(-list$Ratio), ]
constraint_bygene_prot$GENE_SYMBOL <- factor(constraint_bygene_prot$GENE_SYMBOL, c(as.character(list$GENE_SYMBOL)))
constraint_bygene_prot$Data <- factor(constraint_bygene_prot$Data, level = c("ratio_none", "ratio_het", "ratio_hom"))

plotd <- ggplot(constraint_bygene_prot, aes(x = GENE_SYMBOL, y = Ratio, fill = Data))+
  geom_bar(stat = "identity", width = 0.75, position = "fill", colour = "black") + 
  labs(y = "Proportion codons with non-syn. SNVs") +
  scale_fill_manual(values = c("White", "Grey", "Black"), labels = c("Not observed", "At heteroplasmy only", "At homoplasmy")) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 20, angle = 35, hjust = 1),
        axis.title.y = element_text(size = 20, vjust = 0.5, face = "bold"),
        axis.text.y  = element_text(size = 20),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.75), "cm"),
        axis.ticks.length = unit(0.15, "cm"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  guides(fill = FALSE)

save(plotd, plote, plotf, file = "figures/Fig5def.rdata")