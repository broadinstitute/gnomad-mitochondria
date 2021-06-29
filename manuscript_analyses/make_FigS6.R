library(ggplot2)
library(tidyr)
library(ggpubr)

dir.create("figures")

synthetic_vep <- read.delim(file='final_data_files/synthetic_vcf/NC_012920.1_synthetic_vep_splitvarstwogenes.vcf', header=TRUE, sep = "\t", dec = ".")
gnomad <- read.delim(file='reformated.vcf', header=TRUE, sep = "\t", na.strings = "")

#plot A - proportion possible observed by transitions vs transversions
#count number of expected, first remove m.3107N
synthetic_vep <- synthetic_vep[synthetic_vep$REF!="N",]
synthetic_vep$count <- 1
synthetic_vep$type <- ifelse(synthetic_vep$REF=="G"&synthetic_vep$ALT=="A","Transitions",ifelse(synthetic_vep$REF=="C"&synthetic_vep$ALT=="T","Transitions",ifelse(synthetic_vep$REF=="A"&synthetic_vep$ALT=="G","Transitions",ifelse(synthetic_vep$REF=="T"&synthetic_vep$ALT=="C","Transitions","Transversions"))))
synthetic_vep$Consequence <- ifelse(synthetic_vep$BIOTYPE=="Mt_rRNA","rRNA",ifelse(synthetic_vep$BIOTYPE=="Mt_tRNA","tRNA",as.character(synthetic_vep$Consequence)))
#select variants for plot
synthetic_vep <- synthetic_vep[synthetic_vep$Consequence=="synonymous_variant" | synthetic_vep$Consequence=="missense_variant" | synthetic_vep$Consequence=="stop_gained" | synthetic_vep$Consequence=="tRNA" | synthetic_vep$Consequence=="rRNA" | synthetic_vep$POS<577 | synthetic_vep$POS>16023,]
expected <- aggregate(synthetic_vep$count, by=list(Category=synthetic_vep$Consequence,synthetic_vep$type), FUN=sum)
colnames(expected)[c(2,3)] <- c("Type","Expected")

#count number of observed, heteroplasmic separately to homoplasmic
gnomad$count <- 1
gnomad$type <- ifelse(gnomad$REF=="G"&gnomad$ALT=="A","Transitions",ifelse(gnomad$REF=="C"&gnomad$ALT=="T","Transitions",ifelse(gnomad$REF=="A"&gnomad$ALT=="G","Transitions",ifelse(gnomad$REF=="T"&gnomad$ALT=="C","Transitions","Transversions"))))
gnomad$Consequence <- ifelse(gnomad$BIOTYPE=="Mt_rRNA" & !is.na(gnomad$BIOTYPE),"rRNA",ifelse(gnomad$BIOTYPE=="Mt_tRNA" & !is.na(gnomad$BIOTYPE),"tRNA",as.character(gnomad$Consequence)))
gnomad <- gnomad[gnomad$Consequence=="synonymous_variant" | gnomad$Consequence=="missense_variant" | gnomad$Consequence=="stop_gained" | gnomad$Consequence=="tRNA" | gnomad$Consequence=="rRNA" | gnomad$POS<577 | gnomad$POS>16023,]

#count number of observed, SNVs only, at HETEROPLASMY only
gnomad_SNVs_het <- gnomad[gnomad$VARIANT_CLASS=="SNV" & gnomad$max_hl<0.95,]
obs_het <- aggregate(gnomad_SNVs_het$count, by=list(Category2=gnomad_SNVs_het$Consequence,gnomad_SNVs_het$type), FUN=sum)
colnames(obs_het)[c(3)] <- "Observed_het"

#count number of observed, SNVs only, at HOMOPLASMY (>95% heteroplasmy level)
gnomad_SNVs_hom <- gnomad[gnomad$VARIANT_CLASS=="SNV" & gnomad$max_hl>=0.95,]
obs_hom <- aggregate(gnomad_SNVs_hom$count, by=list(Category3=gnomad_SNVs_hom$Consequence,gnomad_SNVs_hom$type), FUN=sum)
colnames(obs_hom)[c(3)] <- "Observed_hom"

#added in, so have the hypervariable sequences (HVS) as separate column, coordinates per MITOMAP
#HVS expected
synthetic_vep <- synthetic_vep[(synthetic_vep$POS>16023 & synthetic_vep$POS<16384) | (synthetic_vep$POS>56 & synthetic_vep$POS<373) | (synthetic_vep$POS>437 & synthetic_vep$POS<575),]
expected2 <- aggregate(synthetic_vep$count, by=list(Category=synthetic_vep$Consequence,synthetic_vep$type), FUN=sum)
colnames(expected2)[c(2,3)] <- c("Type","Expected")
expected2$Category[expected2$Category=="intergenic_variant"] <- "HVS"

#HVS observed
gnomad <- gnomad[(gnomad$POS>16023 & gnomad$POS<16384) | (gnomad$POS>56 & gnomad$POS<373) | (gnomad$POS>437 & gnomad$POS<575),] #note this dataframe includes all HVS already, hence can restrict
gnomad_SNVs_het <- gnomad[gnomad$VARIANT_CLASS=="SNV" & gnomad$max_hl<0.95,]
obs_het2 <- aggregate(gnomad_SNVs_het$count, by=list(Category2=gnomad_SNVs_het$Consequence,gnomad_SNVs_het$type), FUN=sum)
colnames(obs_het2)[c(3)] <- "Observed_het"
obs_het2$Category2[obs_het2$Category2=="intergenic_variant"] <- "HVS"

gnomad_SNVs_hom <- gnomad[gnomad$VARIANT_CLASS=="SNV" & gnomad$max_hl>=0.95,]
obs_hom2 <- aggregate(gnomad_SNVs_hom$count, by=list(Category3=gnomad_SNVs_hom$Consequence,gnomad_SNVs_hom$type), FUN=sum)
colnames(obs_hom2)[c(3)] <- "Observed_hom"
obs_hom2$Category3[obs_hom2$Category3=="intergenic_variant"] <- "HVS"

#rbind to add in HVS, also have to add stop gain in as 0 for Tv
expected <- rbind(expected,expected2)
obs_hom <- rbind(obs_hom,c("stop_gained","Transversions",0),obs_hom2)
obs_het <- rbind(obs_het,c("stop_gained","Transversions",0),obs_het2)

#bind and calculate proportion possible observed ratios
table <- cbind(obs_hom[order(obs_hom$Category3),],obs_het[order(obs_het$Category2),],expected[order(expected$Category),]) #need to sort 
table$ratio_het <- as.numeric(table$Observed_het)/table$Expected #need as numeric after the rbind
table$ratio_hom <- as.numeric(table$Observed_hom)/table$Expected
table$ratio_none <- (1 - table$ratio_hom - table$ratio_het)

#reformat for plotting
plot_table <- table[,c("Category","Type","ratio_het","ratio_hom","ratio_none")] %>% gather(Data,Ratio,ratio_het,ratio_hom,ratio_none)
plot_table$Category <- factor(plot_table$Category, c("intergenic_variant","HVS","synonymous_variant","missense_variant","stop_gained","rRNA","tRNA"))
plot_table$Data <- factor(plot_table$Data,c("ratio_none","ratio_het","ratio_hom"))

plota <- ggplot(data=plot_table, aes(x=Category, y=Ratio, fill=Data)) +
  geom_bar(stat="identity", width=0.75, position="fill",colour="black") +
  labs(y="Proportion possible SNVs observed") +
  scale_y_continuous(breaks=seq(0,1,0.1), expand =c(0.01,0.01)) +
  scale_fill_manual(values=c("White", "Grey","Black"),labels = c("Not observed","At heteroplasmy only","At homoplasmy")) +
  scale_x_discrete(labels=c("Control reg.","HVS","Synonymous ","Missense","Stop gain","rRNA","tRNA")) +
  facet_wrap(~ Type) + 
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_text(size=15, colour="black"),
        axis.title.y = element_text(size=18, vjust=0.5,face="bold"),
        axis.text.y  = element_text(size=15, colour="black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        legend.text = element_text(size=18),
        legend.title = element_blank(), 
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        legend.position="top",
        strip.text.x = element_text(size = 14))

#plot B - proportion observed at heteroplasmy vs homoplasmy - reread in the gnomad file
gnomad <- read.delim(file='reformated.vcf', header=TRUE, sep = "\t", na.strings = "")
#select and rename variants for plot
gnomad$Consequence_clean <- ifelse(gnomad$Consequence=="stop_gained",print("Stop gain\nSNVs"),
                                   ifelse(gnomad$VARIANT_CLASS=="SNV" & gnomad$BIOTYPE=="Mt_tRNA" & !is.na(gnomad$BIOTYPE),print("tRNA\nSNVs"),
                                          ifelse(gnomad$VARIANT_CLASS=="SNV" & gnomad$BIOTYPE=="Mt_rRNA" & !is.na(gnomad$BIOTYPE),print("rRNA\nSNVs"),
                                                 ifelse(gnomad$Consequence=="missense_variant",print("Missense\nSNVs"),
                                                        ifelse(gnomad$Consequence=="synonymous_variant",print("Synonymous\nSNVs"),
                                                               ifelse(gnomad$VARIANT_CLASS!="SNV" & gnomad$BIOTYPE=="protein_coding" & !is.na(gnomad$BIOTYPE),print("Protein\nindels"),
                                                                      ifelse(gnomad$VARIANT_CLASS!="SNV" & gnomad$BIOTYPE=="Mt_tRNA" & !is.na(gnomad$BIOTYPE),print("tRNA\nindels"),
                                                                             ifelse(gnomad$VARIANT_CLASS!="SNV" & gnomad$BIOTYPE=="Mt_rRNA" & !is.na(gnomad$BIOTYPE),print("rRNA\nindels"),
                                                                                    ifelse(gnomad$VARIANT_CLASS=="SNV" & gnomad$Consequence=="intergenic_variant" & (gnomad$POS<577 | gnomad$POS>16023),print("Control reg.\nSNVs"),
                                                                                           ifelse(gnomad$VARIANT_CLASS!="SNV" & gnomad$Consequence=="intergenic_variant" & (gnomad$POS<577 | gnomad$POS>16023),print("Control reg.\nindels"),
                                                                                                  print("other")))))))))))
#count number of observed, heteroplasmic separately to homoplasmic
gnomad$count <- 1
#count number of observed, SNVs only, at HETEROPLASMY only
gnomad_het <- gnomad[gnomad$max_hl<0.95 & gnomad$Consequence_clean!="other",]
obs_het <- aggregate(gnomad_het$count, by=list(Category=gnomad_het$Consequence_clean), FUN=sum)
colnames(obs_het)[c(2)] <- "Observed_het"
#count number of observed, SNVs only, at HOMOPLASMY (>95% heteroplasmy)
gnomad_hom <- gnomad[gnomad$max_hl>=0.95 & gnomad$Consequence_clean!="other",]
obs_hom <- aggregate(gnomad_hom$count, by=list(Category2=gnomad_hom$Consequence_clean), FUN=sum)
colnames(obs_hom)[c(2)] <- "Observed_hom"

#bind and calculate proportion of observed ratios
table <- cbind(obs_hom,obs_het)
plot_table  <- table[,c("Category","Observed_het","Observed_hom")] %>% gather(Data,Number,Observed_het,Observed_hom)
plot_table$Category <- factor(plot_table$Category, levels=c("Control reg.\nSNVs","Synonymous\nSNVs","Missense\nSNVs","Stop gain\nSNVs","rRNA\nSNVs","tRNA\nSNVs", "Control reg.\nindels","Protein\nindels","rRNA\nindels","tRNA\nindels"))

plotb <- ggplot(data=plot_table, aes(x=Category, y=Number, fill=Data)) +
  geom_bar(stat="identity", position="fill", width=0.75,colour="black") +
  labs(y="Proportion observed") +
  scale_y_continuous(breaks=seq(0,1,0.1), expand =c(0.01,0.01)) +
  scale_fill_manual(values=c("Grey", "Black"),labels = c("At heteroplasmy only ","At homoplasmy")) +
  scale_x_discrete(labels=c(sprintf("Control reg.   \nSNVs\nn=%s",sum(plot_table[plot_table$Category=="Control reg.\nSNVs",c("Number")])),
                            sprintf("Synonymous\nSNVs\nn=%s",sum(plot_table[plot_table$Category=="Synonymous\nSNVs",c("Number")])),
                            sprintf("Missense\nSNVs\nn=%s",sum(plot_table[plot_table$Category=="Missense\nSNVs",c("Number")])),
                            sprintf("Stop gain\nSNVs\nn=%s",sum(plot_table[plot_table$Category=="Stop gain\nSNVs",c("Number")])),
                            sprintf("rRNA\nSNVs\nn=%s",sum(plot_table[plot_table$Category=="rRNA\nSNVs",c("Number")])),
                            sprintf("tRNA\nSNVs\nn=%s",sum(plot_table[plot_table$Category=="tRNA\nSNVs",c("Number")])), 
                            sprintf("Control reg.\nindels\nn=%s",sum(plot_table[plot_table$Category=="Control reg.\nindels",c("Number")])),
                            sprintf("Protein\nindels\nn=%s",sum(plot_table[plot_table$Category=="Protein\nindels",c("Number")])),
                            sprintf("rRNA\nindels\nn=%s",sum(plot_table[plot_table$Category=="rRNA\nindels",c("Number")])),
                            sprintf("tRNA\nindels\nn=%s",sum(plot_table[plot_table$Category=="tRNA\nindels",c("Number")])))) +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_text(size=17, colour="black"),
        axis.title.y = element_text(size=18, vjust=0.5,face="bold"),
        axis.text.y  = element_text(size=18, colour="black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        legend.text = element_text(size=22),
        legend.title = element_text(size=22), 
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        legend.position="top") +
  labs(fill="Legend:")

#plot C - in silicos for missense (APOGEE)
#count number of observed, SNVs only, at HETEROPLASMY only
gnomad_mis_het <- gnomad[gnomad$max_hl<0.95 & gnomad$Consequence=="missense_variant" & gnomad$VARIANT_CLASS=="SNV",]
mis_obs_het <- aggregate(gnomad_mis_het$count, by=list(Category=gnomad_mis_het$APOGEE), FUN=sum)
colnames(mis_obs_het)[c(2)] <- "Observed_het"

#count number of observed, SNVs only, at HOMOPLASMY (>95% heteroplasmy)
gnomad_mis_hom <- gnomad[gnomad$max_hl>=0.95 & gnomad$Consequence=="missense_variant" & gnomad$VARIANT_CLASS=="SNV",]
mis_obs_hom <- aggregate(gnomad_mis_hom$count, by=list(Category2=gnomad_mis_hom$APOGEE), FUN=sum) #NOTE there is 1 N/A
colnames(mis_obs_hom)[c(2)] <- "Observed_hom"

#bind and calculate proportion of observed ratios
mis_table <- cbind(mis_obs_hom,mis_obs_het)
mis_plot_table  <- mis_table[,c("Category","Observed_het","Observed_hom")] %>% gather(Data,Number,Observed_het,Observed_hom)

plotc <- ggplot(data=mis_plot_table, aes(x=Category, y=Number, fill=Data)) +
  geom_bar(stat="identity", position="fill", width=0.75,colour="black") +
  labs(y="Proportion observed",x="\nAPOGEE missense in silico prediction") +
  scale_y_continuous(breaks=seq(0,1,0.1), expand =c(0.01,0.01)) +
  scale_fill_manual(values=c("Grey", "Black"),labels = c("At heteroplasmy only ","At homoplasmy")) +
  scale_x_discrete(labels=c(sprintf("Neutral\nn=%s",sum(mis_plot_table[mis_plot_table$Category=="Neutral",c("Number")])),
                            sprintf("Pathogenic\nn=%s",sum(mis_plot_table[mis_plot_table$Category=="Pathogenic",c("Number")])))) + 
  theme(axis.title.x = element_text(size=18, colour="black",face="bold"),
        axis.text.x  = element_text(size=18, colour="black"),
        axis.title.y = element_text(size=18, vjust=0.5,face="bold"),
        axis.text.y  = element_text(size=18, colour="black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        plot.margin = unit(c(0.5,1,0.5,0.5), "cm")) +
  guides(fill=FALSE)

#plot D - in silicos for tRNA (MitoTIP and HmtVAR)
gnomad$tRNA_interp <- ifelse(grepl("pathogenic",gnomad$Mitotip) & grepl("pathogenic",gnomad$Hmtvar),print("Both pathogenic"),
                             ifelse(grepl("benign",gnomad$Mitotip) & grepl("polymorphic",gnomad$Hmtvar),print("Both benign"),
                                    print("Conflicting"))) #this includes the few variants with Hmtvar=="null"

#count number of observed, SNVs only, at HETEROPLASMY only
gnomad_tRNA_het <- gnomad[gnomad$max_hl<0.95 & gnomad$BIOTYPE=="Mt_tRNA" & !is.na(gnomad$BIOTYPE) & gnomad$VARIANT_CLASS=="SNV",]
tRNA_obs_het <- aggregate(gnomad_tRNA_het$count, by=list(Category=gnomad_tRNA_het$tRNA_interp), FUN=sum)
colnames(tRNA_obs_het)[c(2)] <- "Observed_het"

#count number of observed, SNVs only, at HOMOPLASMY (>95% heteroplasmy)
gnomad_tRNA_hom <- gnomad[gnomad$max_hl>=0.95 & gnomad$BIOTYPE=="Mt_tRNA" & !is.na(gnomad$BIOTYPE) & gnomad$VARIANT_CLASS=="SNV",]
tRNA_obs_hom <- aggregate(gnomad_tRNA_hom$count, by=list(Category2=gnomad_tRNA_hom$tRNA_interp), FUN=sum)
colnames(tRNA_obs_hom)[c(2)] <- "Observed_hom"

#bind and calculate proportion of observed ratios
tRNA_table <- cbind(tRNA_obs_hom,tRNA_obs_het)
tRNA_plot_table  <- tRNA_table[,c("Category","Observed_het","Observed_hom")] %>% gather(Data,Number,Observed_het,Observed_hom)
tRNA_plot_table$Category <- factor(tRNA_plot_table$Category, levels=c("Both benign","Conflicting","Both pathogenic"))

plotd <- ggplot(data=tRNA_plot_table, aes(x=Category, y=Number, fill=Data)) +
  geom_bar(stat="identity", position="fill", width=0.75,colour="black") +
  labs(y="Proportion observed",x="\nMitoTIP & HmtVAR tRNA SNV in silico predictions") +
  scale_y_continuous(breaks=seq(0,1,0.1), expand =c(0.01,0.01)) +
  scale_fill_manual(values=c("Grey", "Black"),labels = c("At heteroplasmy only ","At homoplasmy")) +
  scale_x_discrete(labels=c(sprintf("Both benign\nn=%s",sum(tRNA_plot_table[tRNA_plot_table$Category=="Both benign",c("Number")])),
                            sprintf("Conflicting\nn=%s",sum(tRNA_plot_table[tRNA_plot_table$Category=="Conflicting",c("Number")])),
                            sprintf("Both pathogenic\nn=%s",sum(tRNA_plot_table[tRNA_plot_table$Category=="Both pathogenic",c("Number")])))) + 
  theme(axis.title.x = element_text(size=18, colour="black",face="bold"),
        axis.text.x  = element_text(size=18, colour="black"),
        axis.title.y = element_text(size=18, vjust=0.5,face="bold"),
        axis.text.y  = element_text(size=18, colour="black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        plot.margin = unit(c(0.5,0.5,0.5,1), "cm")) +
  guides(fill=FALSE)

ggarrange(plota,plotb,ggarrange(plotc,plotd,labels = c("C","D"),ncol = 2, nrow = 1,font.label = list(size = 22),widths=c(1,1.3)),labels=c("A","B",""),ncol=1,nrow=3,font.label = list(size = 22))

#options(bitmapType = 'cairo', device = 'png')
ggsave("figures/FigureS6.png", width=17, height=15)