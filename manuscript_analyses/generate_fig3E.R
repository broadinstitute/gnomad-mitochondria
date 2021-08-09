#!/broad/software/free/Linux/redhat_7_x86_64/pkgs/r_3.5.0-bioconductor/bin/Rscript
############################
## This script generates paper statistics and figure panels
#
# Input files:
#  all.sample2var.passplus.txt: variants passing filters 
#  summary.unfiltered.var.vep.txt: VEP annotations of all unique variants
#  pos.ref.alt.indel_stack.txt: list of sites subsequently filtered out as indel stacks
#  chrM.pos2annot.txt: tab delimited file of all chrM positions and annotations including Type
#
# Output files:
# plots/Fig3E: stacked bar x variant type (CR, intergenic, tRNA, rRNA, coding) 
############################
library(ggplot2)

############################
# input files
############################
pass=read.delim("all.sample2var.passplus.txt",stringsAsFactors=FALSE)
pass$count=1

######################################################################################
# 3E: stacked bar x variant type (CR, intergenic, tRNA, rRNA, coding) file=var_by_type.combo.pdf
# make normalized stacked bar chart by variant type (CR, intergenic, tRNA, rRNA, coding LOW,MODERATE,HIGH)#
######################################################################################

# step 1: read in annotations of all chrM positions, and get 2 columns: varType and subset
annot=read.delim("chrM.pos2annot.txt",stringsAsFactors=FALSE)
annot[annot$Type=="coding","Type"]="protein-coding (all)"
annot$varType=annot$Type
annot$subset="mtDNA\npositions\nn=16569"

# step 2: get all unique variants (data frame mu) that are PASS, in release, and exclude indel_stack
# assign these varType and subset too
#pass=read.delim("all.sample2var.passplus.txt",stringsAsFactors=FALSE)
m=pass[(pass$HL>=0.10) & (pass$release3.1.1 == "true"),]
indelstack=read.delim("pos.ref.alt.indel_stack.txt",stringsAsFactors=FALSE)
indelstack=indelstack[,1]
m=m[!(m$POS.REF.ALT %in% indelstack),]

# subset = homoplasmic-haplogroup, homoplasmic-non-haplogroup, heteroplasmic-only (add n counts)
m$subset="none"
homMarker=unique(m[(m$HL>=0.95) & (m$hap_defining_variant=="hap_defining_variant"),"POS.REF.ALT"])
homNoMarker=unique(m[(m$HL>=0.95) & (m$hap_defining_variant=="0"),"POS.REF.ALT"])
hetOnly=unique(m[!(m$POS.REF.ALT %in% c(homMarker,homNoMarker)),"POS.REF.ALT"])
m[m$POS.REF.ALT %in% homMarker,"subset"]=paste("homoplasmic\n(haplogroup)\nn=",length(homMarker),sep="")
m[m$POS.REF.ALT %in% homNoMarker,"subset"]=paste("homoplasmic\n(non-haplogroup)\nn=",length(homNoMarker),sep="")
m[m$POS.REF.ALT %in% hetOnly,"subset"]=paste("heteroplasmic\nonly\nn=",length(hetOnly),sep="")

# assign mu to be just unique variants
mu=unique(m[,c("POS.REF.ALT","subset","count")])

# pull out the VEP annotations and assign to categories
vep=read.delim("summary.unfiltered.var.vep.txt",stringsAsFactors=FALSE)
vep$POS=gsub("[.ATCG]","",vep$POS.REF.ALT,perl=TRUE,fixed=FALSE, useBytes=FALSE)
mu$varType="none"
mu[mu$POS.REF.ALT %in% vep[vep$IMPACT %in% c("HIGH"),"POS.REF.ALT"],"varType"]="non-synonymous"
mu[mu$POS.REF.ALT %in% vep[vep$IMPACT %in% c("MODERATE"),"POS.REF.ALT"],"varType"]="non-synonymous"
mu[mu$POS.REF.ALT %in% vep[vep$IMPACT %in% c("LOW"),"POS.REF.ALT"],"varType"]="synonymous"
mu[mu$POS.REF.ALT %in% vep[vep$BIOTYPE %in% c("Mt_tRNA"),"POS.REF.ALT"],"varType"]="tRNA"
mu[mu$POS.REF.ALT %in% vep[vep$BIOTYPE %in% c("Mt_rRNA"),"POS.REF.ALT"],"varType"]="rRNA"
mu[mu$POS.REF.ALT %in% vep[(vep$BIOTYPE=="") & (vep$POS >= 577) & (vep$POS < 16024),"POS.REF.ALT"],"varType"]="intergenic"
mu[mu$POS.REF.ALT %in% vep[(vep$IMPACT=="MODIFIER") & (vep$BIOTYPE=="protein_coding"),"POS.REF.ALT"],"varType"]="intergenic"
mu[mu$POS.REF.ALT %in% vep[(vep$BIOTYPE=="") & !((vep$POS >= 577) & (vep$POS < 16024)),"POS.REF.ALT"],"varType"]="CR"

# step 3: combine the mtDNA annotations plus the unique variant annotations
tmp1=annot[,c("subset","varType")]
tmp2=mu[,c("subset","varType")]
tmp3=rbind(tmp1,tmp2)
# create factors
tmp3$varType=factor(tmp3$varType,levels=c("protein-coding (all)","non-synonymous","synonymous","rRNA","tRNA","intergenic","CR"))
tmp3$subset=factor(tmp3$subset,levels=rev(c("mtDNA\npositions\nn=16569","homoplasmic\n(haplogroup)\nn=4221","homoplasmic\n(non-haplogroup)\nn=4988","heteroplasmic\nonly\nn=1641")))

# now plot
mycolors=c("olivedrab","olivedrab3","olivedrab2","lightsteelblue","lightsteelblue1","gray85","black")
pdf("plots/Fig3E.pdf", width=8, height=4,useDingbats=FALSE)
ggplot(tmp3,aes(x=subset,fill=varType)) +
  geom_bar(position="fill") +
  scale_fill_manual(values=mycolors) +
  theme_classic() +
  coord_flip()
dev.off()

