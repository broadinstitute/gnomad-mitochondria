#!/broad/software/free/Linux/redhat_7_x86_64/pkgs/r_3.5.0-bioconductor/bin/Rscript
############################
## This script generates paper statistics and figure panels
##
# Input files:
#  unfiltered_sample_annotations.annot.txt: annotations per sample
#  all.sample2var.unfiltered.annotNUMTFP.txt: unfiltered variants
#  all.sample2var.passplus.annotNUMFP.txt: variants passing filters
#  TableS1.common_heteroplasmies.txt: common heteroplasmies
#
# Output files:
#  output.stats.figS2.txt
# Figure list (see plots/ directory)
# S2A: top heteroplasmies colored by validated NUMT 
# S2C: 25 NUMT-FP stacked bar by # samples 
# S2D: 12684 heteroplasmy vs theoretical 
# S2E: numtA heteroplasmy vs % NUMT-FP (color by mtCN)
# S2F: numtB heteroplasmy vs % NUMT-FP (color by mtCN)
############################
library(ggplot2)
library(RColorBrewer)

############################
# input files
############################
allm=read.delim("all.sample2var.unfiltered.annotNUMTFP.txt",stringsAsFactors=FALSE)
pass=read.delim("all.sample2var.passplus.annotNUMTFP.txt",stringsAsFactors=FALSE)
allm$count=1
pass$count=1
passhets=pass[pass$HL<0.95,]
# for passhets table, set mtcnbin and mtcnbinlabel (up to 500) and vafbin

# sample annotation file
samples=read.delim("unfiltered_sample_annotations.annot.txt",stringsAsFactors=FALSE)
samples$count=1

commonhets=read.delim("TableS1.common_heteroplasmies.txt",stringsAsFactors=FALSE)
numts=commonhets[commonhets$SIGNIFICANT,"POS.REF.ALT"]


############################
# Statistics output file
############################
outf="output.stats.figS2.txt"
line="Statistics for gnomAD"
write(line,file=outf)


############################
# S2E: numtA heteroplasmy vs % NUMT-FP (color by mtCN) 
# S2F: numtB heteroplasmy vs % NUMT-FP (color by mtCN) 
############################
m=passhets
m=m[(m$HL < 0.1599) & (m$HL >= 0.01),]
m=m[!is.na(m$mtCN),]
m$vafbin=floor(m$HL*100)/100
m$mtcnbin=floor(m$mtCN / 25)*25
m[(m$mtcnbin >= 225) & (m$mtcnbin <=500),"mtcnbin"]=225
m[m$mtcnbin > 500,"mtcnbin"]=500
m$mtcnbinlabel=paste(m$mtcnbin,m$mtcnbin+25,sep="-")
m[m$mtcnbinlabel=="225-250","mtcnbinlabel"]="225-500"
m[m$mtcnbinlabel=="500-525","mtcnbinlabel"]="500+"
m$mtcnbinlabel=factor(m$mtcnbinlabel,levels=unique(m$mtcnbinlabel)[order(unique(m$mtcnbin))])
m$mtcnbin=factor(m$mtcnbin,levels=sort(unique(m$mtcnbin)))
#aggregate by vafbin and mtcnbinlabel
tmp1=aggregate(m$count, by=list(vafbin=m$vafbin,mtcnbinlabel=m$mtcnbinlabel), FUN=sum)
colnames(tmp1)[3]="n"
mnumt=m[m$linked.NUMT.FP ==1,]
tmp2=aggregate(mnumt$count, by=list(vafbin=mnumt$vafbin,mtcnbinlabel=mnumt$mtcnbinlabel), FUN=sum)
colnames(tmp2)[3]="n.linkedNUMTFP"
tmp3=merge(tmp1,tmp2,all=TRUE,by=c("vafbin","mtcnbinlabel"))

mnumt=m[(m$linked.NUMT.FP ==1) & (m$numtID == "numtA"),]
tmp2=aggregate(mnumt$count, by=list(vafbin=mnumt$vafbin,mtcnbinlabel=mnumt$mtcnbinlabel), FUN=sum)
colnames(tmp2)[3]="n.numtA"
tmp3=merge(tmp3,tmp2,all=TRUE,by=c("vafbin","mtcnbinlabel"))

mnumt=m[(m$linked.NUMT.FP ==1) & (m$numtID == "numtB"),]
tmp2=aggregate(mnumt$count, by=list(vafbin=mnumt$vafbin,mtcnbinlabel=mnumt$mtcnbinlabel), FUN=sum)
colnames(tmp2)[3]="n.numtB"
tmp3=merge(tmp3,tmp2,all=TRUE,by=c("vafbin","mtcnbinlabel"))

tmp3=tmp3[tmp3$n>0,]
tmp3$percA=tmp3$n.numtA/tmp3$n
tmp3$percB=tmp3$n.numtB/tmp3$n

tmp3[is.na(tmp3)]=0
xlabs=paste(c(1:15)/100,c(2:16)/100,sep="-")
pdf("plots/FigS2E.pdf", width=4, height=4)
ggplot(tmp3,aes(y=percA,x=vafbin,color=mtcnbinlabel)) +
  geom_line() +
  geom_point(size=1) +
  ggtitle("% NUMT by VAF and mtCN") +
  geom_vline(xintercept=0.05, linetype="dashed", color = "black") +
  geom_vline(xintercept=0.15, linetype="dashed", color = "black") +
  geom_vline(xintercept=0.1, linetype="dashed", color = "black") +
  scale_x_continuous(name="Heteroplasmy bin", limits=c(0.01, 0.15),breaks=c(c(1:15)/100),labels=xlabs) +
  scale_y_continuous(name="% NUMT-FP from numtA") +
  labs(color = "mtCN") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90)) +
  scale_color_brewer(palette="Spectral")
dev.off()
pdf("plots/FigS2F.pdf", width=4, height=4)
ggplot(tmp3,aes(y=percB,x=vafbin,color=mtcnbinlabel)) +
  geom_line() +
  geom_point(size=1) +
  ggtitle("% NUMT by VAF and mtCN") +
  geom_vline(xintercept=0.05, linetype="dashed", color = "black") +
  geom_vline(xintercept=0.15, linetype="dashed", color = "black") +
  geom_vline(xintercept=0.1, linetype="dashed", color = "black") +
  scale_x_continuous(name="Heteroplasmy bin", limits=c(0.01, 0.15),breaks=c(c(1:15)/100),labels=xlabs) +
  scale_y_continuous(name="% NUMT-FP from numtB") +
  labs(color = "mtCN") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90)) +
  scale_color_brewer(palette="Spectral")
dev.off()

########################################################
# Plot HL vs 1/(1+mtCN) for 12684.G.A, and 16293.A.C
# FigS2D: 12684 heteroplasmy vs theoretical
########################################################
myhet="12684.G.A"
 subset=passhets[(passhets$POS.REF.ALT == myhet) & (passhets$HL < 0.50) & (passhets$HL >= 0.01),]
 subset$expected=1/(1+subset$mtCN)
 subset$mtcncolor=subset$mtCN
 subset[subset$mtcncolor > 100,"mtcncolor"]=100
 n=dim(subset)[[1]]
 pdf(paste("plots/FigS2D.pdf",sep=""), width=2.5, height=2,useDingbats=FALSE)
 ggplot(subset,aes(y=HL,x=expected,color=mtcncolor)) +
  geom_point(size=0.5) +
  theme_classic() +
  xlab("Expected VAF for heterozygous NUMT\n1/(1+mtCN)") +
  ylab("Observed VAF") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))
 dev.off()
line=paste("\nStatistics for NUMT-FP m.12684G>A:\nn=",
  n,"; Spearman r=",signif(cor(x=subset$HL,y=subset$expected,use="pairwise.complete.obs")[1],2),
  "; p=",cor.test(x=subset$HL,y=subset$expected,method="spearman",alternative="two.sided",exact=FALSE)$p.value,
  "; Note p=0 indicates p<2.2e-16",
  sep="")
write(line,file=outf,append=TRUE)


########################################################
# FigS2A: top heteroplasmies colored by validated NUMT 
########################################################
m=pass[(pass$HL >= 0.01) & (pass$HL < 0.50),]
tmp1=aggregate(m$count, by=list(POS.REF.ALT=m$POS.REF.ALT,numtID=m$numtID), FUN=sum)
tmp1=tmp1[order(tmp1$x,decreasing=TRUE),]
tmp1=tmp1[1:52,]
mylabels=c("Heteroplasmy correlates with 1/(1+mtCN), validated via PacBio",
  "Heteroplasmy correlates with 1/(1+mtCN), no validation data available",
  "Heteroplasmy does not correlate with 1/(1+mtCN)")
tmp1$label=mylabels[3]
tmp1[tmp1$POS.REF.ALT %in% numts,"label"]=mylabels[2]
tmp1[(tmp1$POS.REF.ALT %in% numts) & (tmp1$numtID != ""),"label"]=mylabels[1]
tmp1$label=factor(tmp1$label,levels=mylabels)
tmp1$POS.REF.ALT=factor(tmp1$POS.REF.ALT,levels=tmp1$POS.REF.ALT)
tmp1$letter=""
tmp1[tmp1$numtID=="numtA","letter"]="A"
tmp1[tmp1$numtID=="numtB","letter"]="B"

pdf("plots/FigS2A.pdf", width=15, height=3)
ggplot(tmp1, aes(fill=label, y=x, x=POS.REF.ALT)) +
  geom_hline(yintercept=5000, color = "gray") +
  geom_hline(yintercept=10000, color = "gray") +
  geom_hline(yintercept=15000, color = "gray") +
  geom_bar(position="stack", stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90)) +
  xlab("Variants") +
  ylab("# Samples with variant\n(heteroplasmy 1-50%)") +
  scale_fill_manual(values=c("red","orange","gray")) +
  geom_text(data=tmp1,aes(x=POS.REF.ALT,y=x,label=letter),vjust=-0.1)
dev.off()

########################################################
# Fig S2C: 25 NUMT-FP stacked bar by # samples 
########################################################
numtfp=allm[(allm$numtID != ""),]
numtfp$PASS="fail"
numtfp[numtfp$SAMPLE.VAR %in% pass[(pass$numtID != ""),"SAMPLE.VAR"],"PASS"]="PASS"
numtfp$vartype="no variant"
numtfp[numtfp$HL>=0.95,"vartype"]="homoplasmy"
numtfp[(numtfp$HL<0.95) & (numtfp$PASS == "PASS"),"vartype"]="PASS heteroplasmy"
numtfp[numtfp$HL<0.95 & (numtfp$PASS == "fail"),"vartype"]="FAIL heteroplasmy"
# aggregate by POS.REF.ALT
tmp1=aggregate(numtfp$count, by=list(POS.REF.ALT=numtfp$POS.REF.ALT,vartype=numtfp$vartype), FUN=sum)
tmp2=aggregate(numtfp$count, by=list(POS.REF.ALT=numtfp$POS.REF.ALT), FUN=sum)
colnames(tmp2)[2]="nWithVariant"
tmp2$nNoVariant=nrow(samples)-tmp2$nWithVariant
colnames(tmp2)=c("POS.REF.ALT","vartype","x")
tmp2$vartype="no variant"
# concatentate the known variants (homoplasmy, PASS heteroplasmy, FAIL heteroplasmy) with the no-variant
tmp3=rbind(tmp1,tmp2)
tmp3$fraction=tmp3$x/nrow(samples)
tmp3$vartype=factor(tmp3$vartype,levels=c("no variant","homoplasmy","FAIL heteroplasmy","PASS heteroplasmy"))
pdf("plots/FigS2C.pdf", width=5, height=3)
ggplot(tmp3, aes(fill=vartype, y=fraction, x=POS.REF.ALT)) +
  geom_bar(position="stack", stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90)) +
  labs(fill = "Variant Type") +
  xlab("NUMT-FP variants") +
  ylab("# variants") +
  scale_fill_manual(values=c("gray90","gray70","gray50","black")) +
  geom_hline(yintercept=0.40, linetype="dashed", color = "gray")
dev.off()
