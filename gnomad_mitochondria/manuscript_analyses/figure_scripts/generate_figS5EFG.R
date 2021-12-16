library(ggplot2)
library(reshape2)

############################
# read in 3 databases and merge;
# For each database, create universal ID as M-POS-REF-ALT
############################
gnomad=read.delim("final_data_files/gnomad/gnomad.genomes.v3.1.sites.chrM.reduced_annotations.tsv",stringsAsFactors=FALSE)
mitomap=read.delim("final_data_files/other_databases/MITOMAP_polymorphisms_02022021.cgi",stringsAsFactors=FALSE)
helix=read.delim("final_data_files/other_databases/HelixMTdb_20200327.tsv",stringsAsFactors=FALSE)

# gnomad: get pass sites only
gnomad=gnomad[gnomad$filters == "PASS",]
gnomad$REF.ALT=paste(gnomad$ref,gnomad$alt,sep="-")
gnomad$ID=paste("M",gnomad$position,gnomad$REF.ALT,sep="-")
gnomad$DB="gnomAD"
gnomad=gnomad[,c("ID","REF.ALT","AC_hom","AF_hom","AC_het","AF_het","DB")]

# mitomap
#note mitomap ":" indicates deletion, replace with "del" so it does not look like a SNV
mitomap[mitomap$alt == ":","alt"]="del"
mitomap$REF.ALT=paste(mitomap$ref,mitomap$alt,sep="-")
mitomap$ID=paste("M",mitomap$pos,mitomap$REF.ALT,sep="-")
mitomap$AC_hom=mitomap$gbcnt
mitomap$AF_hom=mitomap$gbfreq/100
mitomap$AC_het=0
mitomap$AF_het=0
mitomap$DB="MITOMAP"
mitomap=mitomap[,c("ID","REF.ALT","AC_hom","AF_hom","AC_het","AF_het","DB")]

# helix
mylen=nchar(helix$alleles)
tmp1=gsub("chrM:","M-",helix$locus)
helix$REF.ALT=gsub(",","-",substr(helix$alleles,2,mylen-1),fixed=TRUE)
helix$ID=paste(tmp1,helix$REF.ALT,sep="-")
helix$AC_hom=helix$counts_hom
helix$AF_hom=helix$AF_hom
helix$AC_het=helix$counts_het
helix$AF_het=helix$AF_het
helix$DB="HelixMTdb"
helix=helix[,c("ID","REF.ALT","AC_hom","AF_hom","AC_het","AF_het","DB")]

# merge these 3 databases together in a "long form" table
dflong=rbind(gnomad,mitomap,helix)
# annotate Type as SNV (if both REF and ALT are single characters, otherwise indel)
dflong$Type="indel"
dflong[nchar(dflong$REF.ALT)==3,"Type"]="SNV"
# remove REF.ALT column
dflong=dflong[,c("ID","Type","AC_hom","AF_hom","AC_het","AF_het","DB")]

# now convert long table "dflong" to wide table "df" which has one row per variant
dcastRename<-function(mat,mystr) {
  tmp=dcast(mat,ID + Type ~ DB,value.var=mystr)
  colnames(tmp)[3:5]=paste(colnames(tmp)[3:5],mystr,sep=".")
  tmp[is.na(tmp)]=0
  tmp
}
tmp.AC_hom=dcastRename(dflong,"AC_hom")
tmp.AF_hom=dcastRename(dflong,"AF_hom")
tmp.AC_het=dcastRename(dflong,"AC_het")
tmp.AF_het=dcastRename(dflong,"AF_het")

df=cbind(tmp.AC_hom,tmp.AF_hom[,3:5],tmp.AC_het[,3:5],tmp.AF_het[,3:5])

# remove entries with 0 counts
df=df[df$HelixMTdb.AC_hom + df$MITOMAP.AC_hom + df$gnomAD.AC_hom + df$HelixMTdb.AC_het + df$MITOMAP.AC_het + df$gnomAD.AC_het > 0,]

mycolors=c("black","cyan","magenta")

############################
# how many variants unique to gnomad not in mitomap/helix
############################
dim(df[(df$gnomAD.AC_hom + df$gnomAD.AC_het > 0) & (df$MITOMAP.AC_hom + df$HelixMTdb.AC_hom + df$HelixMTdb.AC_het == 0),])

############################
# comparison gnomad vs mitomap
############################
mat=df[(df$gnomAD.AC_hom + df$MITOMAP.AC_hom) > 0,]

# plot AF gnomAD vs mitomap just for SNVs, since MITMAP right-aligns indels
matSNV=mat[mat$Type == "SNV",]
set=(matSNV$gnomAD.AC_hom > 0) & (matSNV$MITOMAP.AC_hom >0)
matSNV[set,"category"]= paste("shared N=",sum(set),sep="")
set=(matSNV$gnomAD.AC_hom > 0) & (matSNV$MITOMAP.AC_hom ==0)
matSNV[set,"category"]= paste("gnomAD-only N=",sum(set),sep="")
set=(matSNV$gnomAD.AC_hom == 0) & (matSNV$MITOMAP.AC_hom >0)
matSNV[set,"category"]= paste("MITOMAP-only N=",sum(set),sep="")
matSNV$category=factor(matSNV$category,levels=rev(sort(unique(matSNV$category))))
pdf("plots/FigS5E.pdf", width=5, height=4)
ggplot(matSNV,aes(y=MITOMAP.AF_hom,x=gnomAD.AF_hom,color=category)) +
  geom_point(size=1) +
  ggtitle("SNV allele frequency gnomAD vs MITOMAP") +
  labs(color = "Homoplasmic SNV")+
  theme_classic() +
  xlab("gnomAD v3 homoplasmic allele frequency") +
  ylab("MITOMAP homoplasmic allele frequency")+
  scale_color_manual(values=mycolors)
dev.off()

# statistics
# how many shared SNVs = 7682
sum((mat$gnomAD.AC_hom > 0) & (mat$Type == "SNV") & (mat$MITOMAP.AC_hom >0))
# how many more homoplasmic SNVs does gnomAD have over mitomap? 1222
sum((mat$gnomAD.AC_hom > 0) & (mat$Type == "SNV") & (mat$MITOMAP.AC_hom ==0))
# what % increase of new SNVs does gnomAD have over MITOMAP? 0.10
signif(sum((mat$gnomAD.AC_hom > 0) & (mat$Type == "SNV") & (mat$MITOMAP.AC_hom ==0))/sum(mat$MITOMAP.AC_hom >0),2)

# calculate Pearson (0.98) and Spearman (0.80) correlation for the shared homoplasmic SNVs
matSNV=matSNV[(matSNV$gnomAD.AC_hom > 0) & (matSNV$MITOMAP.AC_hom >0),]
cor(x=matSNV$gnomAD.AF_hom,y=matSNV$MITOMAP.AF_hom,method="pearson")
cor(x=matSNV$gnomAD.AF_hom,y=matSNV$MITOMAP.AF_hom,method="spearman")

############################
# comparison gnomad vs helix homoplasmic
############################
mat=df[(df$gnomAD.AC_hom + df$HelixMTdb.AF_hom ) > 0,]

# plot AF gnomAD vs helix hom
set=(mat$gnomAD.AC_hom > 0) & (mat$HelixMTdb.AC_hom >0)
mat[set,"category"]= paste("shared N=",sum(set),sep="")
set=(mat$gnomAD.AC_hom > 0) & (mat$HelixMTdb.AC_hom ==0)
mat[set,"category"]= paste("gnomAD-only N=",sum(set),sep="")
set=(mat$gnomAD.AC_hom == 0) & (mat$HelixMTdb.AC_hom >0)
mat[set,"category"]= paste("HelixMTdb-only N=",sum(set),sep="")
mat$category=factor(mat$category,levels=rev(sort(unique(mat$category))))
pdf("plots/FigS5F.pdf", width=5, height=4)
ggplot(mat,aes(y=HelixMTdb.AF_hom,x=gnomAD.AF_hom,color=category)) +
  geom_point(size=1) +
  ggtitle("Variant allele frequency gnomAD vs HelixMTdb") +
  labs(color = "Homoplasmic variant")+
  theme_classic() +
  xlab("gnomAD v3 homoplasmic allele frequency") +
  ylab("HelixMTdb homoplasmic allele frequency")+
  scale_color_manual(values=mycolors)
dev.off()

# statistics
# how many shared SNVs = 8497
sum((mat$gnomAD.AC_hom > 0) & (mat$HelixMTdb.AC_hom >0))
# how many more homoplasmic variants does gnomAD have over helixMTdb? 712
sum((mat$gnomAD.AC_hom > 0)  & (mat$HelixMTdb.AC_hom ==0))
# what % increase of new SNVs does gnomAD have over HelixMTdb? 0.06
signif(sum((mat$gnomAD.AC_hom > 0) & (mat$HelixMTdb.AC_hom ==0))/sum(mat$HelixMTdb.AC_hom >0),2)

# calculate Pearson (0.97) and Spearman (0.88) correlation for the shared homoplasmic variants
mat=mat[(mat$gnomAD.AC_hom > 0) & (mat$HelixMTdb.AC_hom >0),]
cor(x=mat$gnomAD.AF_hom,y=mat$HelixMTdb.AF_hom,method="pearson")
cor(x=mat$gnomAD.AF_hom,y=mat$HelixMTdb.AF_hom,method="spearman")

############################
# comparison gnomad vs helix het
############################
mat=df[(df$gnomAD.AC_het + df$HelixMTdb.AF_het ) > 0,]

# plot AF gnomAD vs other database
set=(mat$gnomAD.AC_het > 0) & (mat$HelixMTdb.AC_het >0)
mat[set,"category"]= paste("shared N=",sum(set),sep="")
set=(mat$gnomAD.AC_het > 0) & (mat$HelixMTdb.AC_het ==0)
mat[set,"category"]= paste("gnomAD-only N=",sum(set),sep="")
set=(mat$gnomAD.AC_het == 0) & (mat$HelixMTdb.AC_het >0)
mat[set,"category"]= paste("HelixMTdb-only N=",sum(set),sep="")
mat$category=factor(mat$category,levels=rev(sort(unique(mat$category))))
pdf("plots/FigS5G.pdf", width=5, height=4)
ggplot(mat,aes(y=HelixMTdb.AF_het,x=gnomAD.AF_het,color=category)) +
  geom_point(size=1) +
  ggtitle("Variant allele frequency gnomAD vs HelixMTdb") +
  labs(color = "Heteroplasmic variant")+
  theme_classic() +
  xlab("gnomAD v3 heteroplasmic allele frequency") +
  ylab("HelixMTdb heteroplasmic allele frequency")+
  scale_color_manual(values=mycolors) + 
  scale_x_continuous(limits=c(5e-6,0.03),trans='log10') + 
  scale_y_continuous(limits=c(5e-6,0.03),trans='log10') 
dev.off()

# statistics
# how many shared variants = 4334
sum((mat$gnomAD.AC_het > 0) & (mat$HelixMTdb.AC_het >0))
# how many more heteroplasmic variants does gnomAD have over helixMTdb? 1311
sum((mat$gnomAD.AC_het > 0) & (mat$HelixMTdb.AC_het ==0))
# what % increase of new variants does gnomAD have over HelixMTdb? 0.13
signif(sum((mat$gnomAD.AC_het > 0) & (mat$HelixMTdb.AC_het ==0))/sum(mat$HelixMTdb.AC_het >0),2)

# calculate Pearson (0.45) and Spearman (0.67) correlation for the shared homoplasmic SNVs
mat=mat[(mat$gnomAD.AC_het > 0) & (mat$HelixMTdb.AC_het >0),]
cor(x=mat$gnomAD.AF_het,y=mat$HelixMTdb.AF_het,method="pearson")
cor(x=mat$gnomAD.AF_het,y=mat$HelixMTdb.AF_het,method="spearman")


