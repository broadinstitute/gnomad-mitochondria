library(ggplot2)

############################
# input files: compapison of AC and AF in gnomAD v3, MITOMAP (downloaded 2021-12-01, and HelixMTdb downloaded 2021-12-01
############################
df <- read.delim(file = 'final_data_files/other_databases/AC.HelixMTDb.MITOMAP.12012021.txt', stringsAsFactors=FALSE)

############################
# comparison gnomad vs mitomap
############################
mat=df[(df$gnomAD_AC_hom + df$MITOMAP_AC_hom) > 0,]

# plot AF gnomAD vs mitomap just for SNVs, since MITMAP right-aligns indels
matSNV=mat[mat$Type == "SNV",]
set=(matSNV$gnomAD_AC_hom > 0) & (matSNV$MITOMAP_AC_hom >0)
matSNV[set,"category"]= paste("shared N=",sum(set),sep="")
set=(matSNV$gnomAD_AC_hom > 0) & (matSNV$MITOMAP_AC_hom ==0)
matSNV[set,"category"]= paste("gnomAD-only N=",sum(set),sep="")
set=(matSNV$gnomAD_AC_hom == 0) & (matSNV$MITOMAP_AC_hom >0)
matSNV[set,"category"]= paste("MITOMAP-only N=",sum(set),sep="")
matSNV$category=factor(matSNV$category,levels=rev(unique(matSNV$category)))
pdf("plots/FigS5E.pdf", width=5, height=4)
ggplot(matSNV,aes(y=MITOMAP_AF_hom,x=gnomAD_AF_hom,color=category)) +
  geom_point(size=1) +
  ggtitle("SNV allele frequency gnomAD vs MITOMAP") +
  labs(color = "Homoplasmic SNV")+
  theme_classic() +
  xlab("gnomAD v3 homoplasmic allele frequency") +
  ylab("MITOMAP homoplasmic allele frequency")+
  scale_color_manual(values=c("cyan","black", "magenta"))
dev.off()

# statistics
# how many more homoplasmic SNVs does gnomAD have over mitomap? 1186
sum((mat$gnomAD_AC_hom > 0) & (mat$Type == "SNV") & (mat$MITOMAP_AC_hom ==0))
# what % increase of new SNVs does gnomAD have over MITOMAP? 0.099
signif(sum((mat$gnomAD_AC_hom > 0) & (mat$Type == "SNV") & (mat$MITOMAP_AC_hom ==0))/sum(mat$MITOMAP_AC_hom >0),2)

# calculate Pearson (0.98) and Spearman (0.80) correlation for the shared homoplasmic SNVs
matSNV=matSNV[(matSNV$gnomAD_AC_hom > 0) & (matSNV$MITOMAP_AC_hom >0),]
cor(x=matSNV$gnomAD_AF_hom,y=matSNV$MITOMAP_AF_hom,method="pearson")
cor(x=matSNV$gnomAD_AF_hom,y=matSNV$MITOMAP_AF_hom,method="spearman")

############################
# comparison gnomad vs helix homoplasmic
############################
mat=df[(df$gnomAD_AC_hom + df$HelixMTdb_AF_hom ) > 0,]

# plot AF gnomAD vs helix hom
set=(mat$gnomAD_AC_hom > 0) & (mat$HelixMTdb_AC_hom >0)
mat[set,"category"]= paste("shared N=",sum(set),sep="")
set=(mat$gnomAD_AC_hom > 0) & (mat$HelixMTdb_AC_hom ==0)
mat[set,"category"]= paste("gnomAD-only N=",sum(set),sep="")
set=(mat$gnomAD_AC_hom == 0) & (mat$HelixMTdb_AC_hom >0)
mat[set,"category"]= paste("HelixMTdb-only N=",sum(set),sep="")
mat$category=factor(mat$category,levels=rev(unique(mat$category)))
pdf("plots/FigS5F.pdf", width=5, height=4)
ggplot(mat,aes(y=HelixMTdb_AF_hom,x=gnomAD_AF_hom,color=category)) +
  geom_point(size=1) +
  ggtitle("Variant allele frequency gnomAD vs HelixMTdb") +
  labs(color = "Homoplasmic variant")+
  theme_classic() +
  xlab("gnomAD v3 homoplasmic allele frequency") +
  ylab("HelixMTdb homoplasmic allele frequency")+
  scale_color_manual(values=c("cyan","black", "magenta"))
dev.off()

# statistics
# how many more homoplasmic variants does gnomAD have over helixMTdb? 712
sum((mat$gnomAD_AC_hom > 0)  & (mat$HelixMTdb_AC_hom ==0))
# what % increase of new SNVs does gnomAD have over HelixMTdb? 0.06
signif(sum((mat$gnomAD_AC_hom > 0) & (mat$HelixMTdb_AC_hom ==0))/sum(mat$HelixMTdb_AC_hom >0),2)

# calculate Pearson (0.97) and Spearman (0.88) correlation for the shared homoplasmic SNVs
mat=mat[(mat$gnomAD_AC_hom > 0) & (mat$HelixMTdb_AC_hom >0),]
cor(x=mat$gnomAD_AF_hom,y=mat$HelixMTdb_AF_hom,method="pearson")
cor(x=mat$gnomAD_AF_hom,y=mat$HelixMTdb_AF_hom,method="spearman")

############################
# comparison gnomad vs helix het
############################
mat=df[(df$gnomAD_AC_het + df$HelixMTdb_AF_het ) > 0,]

# plot AF gnomAD vs other database
set=(mat$gnomAD_AC_het > 0) & (mat$HelixMTdb_AC_het >0)
mat[set,"category"]= paste("shared N=",sum(set),sep="")
set=(mat$gnomAD_AC_het > 0) & (mat$HelixMTdb_AC_het ==0)
mat[set,"category"]= paste("gnomAD-only N=",sum(set),sep="")
set=(mat$gnomAD_AC_het == 0) & (mat$HelixMTdb_AC_het >0)
mat[set,"category"]= paste("HelixMTdb-only N=",sum(set),sep="")
mat$category=factor(mat$category,levels=rev(unique(mat$category)))
pdf("plots/FigS5G.pdf", width=5, height=4)
ggplot(mat,aes(y=HelixMTdb_AF_het,x=gnomAD_AF_het,color=category)) +
  geom_point(size=1) +
  ggtitle("Variant allele frequency gnomAD vs HelixMTdb") +
  labs(color = "Heteroplasmic variant")+
  theme_classic() +
  xlab("gnomAD v3 heteroplasmic allele frequency") +
  ylab("HelixMTdb heteroplasmic allele frequency")+
  scale_color_manual(values=c("cyan","black", "magenta")) +
  scale_x_continuous(limits=c(5e-6,0.03),trans='log10') + 
  scale_y_continuous(limits=c(5e-6,0.03),trans='log10') 
dev.off()

# statistics
# how many more heteroplasmic variants does gnomAD have over helixMTdb? 1198
sum((mat$gnomAD_AC_het > 0) & (mat$HelixMTdb_AC_het ==0))
# what % increase of new SNVs does gnomAD have over HelixMTdb? 0.12
signif(sum((mat$gnomAD_AC_het > 0) & (mat$HelixMTdb_AC_het ==0))/sum(mat$HelixMTdb_AC_het >0),2)

# calculate Pearson (0.45) and Spearman (0.61) correlation for the shared homoplasmic SNVs
mat=mat[(mat$gnomAD_AC_het > 0) & (mat$HelixMTdb_AC_het >0),]
cor(x=mat$gnomAD_AF_het,y=mat$HelixMTdb_AF_het,method="pearson")
cor(x=mat$gnomAD_AF_het,y=mat$HelixMTdb_AF_het,method="spearman")


