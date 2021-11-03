#!/broad/software/free/Linux/redhat_7_x86_64/pkgs/r_3.5.0-bioconductor/bin/Rscript
############################
## This script generates paper statistics and figure panels
##
# Input files:
#  unfiltered_sample_annotations.annot.txt: annotations per sample
#  all.sample2var.passplus.annotFP.txt: variants passing filters 
#
# Figure list (see plots/ directory)
# FigS3A: mtCN bin x # hets/sample by insert size 
# FigS3B: mtCN bin x # NUMT-FP hets/sample by insert size
############################
library(ggplot2)
library(RColorBrewer)

############################
# input files
############################
pass=read.delim("all.sample2var.passplus.annotNUMTFP.txt",stringsAsFactors=FALSE)
pass$count=1
passhets=pass[pass$HL<0.95,]
# for passhets table, set mtcnbin and mtcnbinlabel (up to 500) and vafbin
passhets=passhets[!is.na(passhets$mtCN),]
passhets$vafbin=floor(passhets$HL*100)/100
passhets$mtcnbin=floor(passhets$mtCN / 25)*25
passhets[passhets$mtcnbin > 500,"mtcnbin"]=500
passhets$mtcnbinlabel=paste(passhets$mtcnbin,passhets$mtcnbin+25,sep="-")
passhets[passhets$mtcnbinlabel=="500-525","mtcnbinlabel"]="500+"
passhets$mtcnbinlabel=factor(passhets$mtcnbinlabel,levels=unique(passhets$mtcnbinlabel)[order(unique(passhets$mtcnbin))])
passhets$mtcnbin=factor(passhets$mtcnbin,levels=sort(unique(passhets$mtcnbin)))

# sample annotation file
samples=read.delim("unfiltered_sample_annotations.annot.txt",stringsAsFactors=FALSE)
samples$count=1
# annotate mtcnbin and mtcnbinlabel
samples$mtcnbin=floor(samples$mtCN / 25)*25
samples[samples$mtcnbin > 500,"mtcnbin"]=500
samples$mtcnbinlabel=paste(samples$mtcnbin,samples$mtcnbin+25,sep="-")
samples[samples$mtcnbinlabel=="500-525","mtcnbinlabel"]="500+"
samples$mtcnbinlabel=factor(samples$mtcnbinlabel,levels=unique(samples$mtcnbinlabel)[order(unique(samples$mtcnbin))])
samples$mtcnbin=factor(samples$mtcnbin,levels=sort(unique(samples$mtcnbin)))
# annotate insertsizebin and insertsizelabel
samples[is.na(samples$median_insert_size),"median_insert_size"]=357
samples$insertsizebin=floor(samples$median_insert_size / 50)*50
samples[samples$insertsizebin > 450,"insertsizebin"]=450
samples[samples$insertsizebin < 300,"insertsizebin"]=300
samples$insertsizebinlabel=paste(samples$insertsizebin,samples$insertsizebin+50,sep="-")
samples[samples$insertsizebinlabel=="450-500","insertsizebinlabel"]="450+"
samples$insertsizebinlabel=factor(samples$insertsizebinlabel,levels=unique(samples$insertsizebinlabel)[order(unique(samples$insertsizebin))])
samples$insertsizebin=factor(samples$insertsizebin,levels=sort(unique(samples$insertsizebin)))


############################
# insert size panels
# S3A: mtCN bin x # hets/sample by insert size 
# S3B: mtCN bin x # NUMT-FP hets/sample by insert size
############################
# add in insertsize
m=passhets
m=m[(m$HL < 0.50) & (m$HL >= 0.01),]
tmp=samples[,c("participant_id","insertsizebinlabel")]
m=merge(m,tmp)

# only show mtCN<300
m=m[m$mtCN <= 300,]

# calculate mean # hets/sample by mtCN bin
tmp1=aggregate(samples$count, by=list(mtcnbinlabel=samples$mtcnbinlabel,insertsizebinlabel=samples$insertsizebinlabel), FUN=sum)
colnames(tmp1)[3]="n"
tmp2=aggregate(m$count, by=list(mtcnbinlabel=m$mtcnbinlabel,insertsizebinlabel=m$insertsizebinlabel), FUN=sum)
colnames(tmp2)[3]="n.hets1to50"
tmp3=merge(tmp2,tmp1,by=c("mtcnbinlabel","insertsizebinlabel"))
tmp3$mean.hets1to50=tmp3$n.hets1to50/tmp3$n
pdf("plots/FigS3A.pdf", width=5, height=4,useDingbats=FALSE)
ggplot(tmp3,aes(y=mean.hets1to50,x=mtcnbinlabel,color=insertsizebinlabel)) +
  geom_point(size=1) +
  ggtitle("Mean # hets/sample by mtCN and insert-size") +
  xlab("mtDNA copy number (mtCN)") +
  ylab("Mean # variants/sample\n(heteroplasmy 1-50%)") +
  labs(color = "insert-size") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90)) +
  scale_color_manual(values=brewer.pal(n = 6, name = "Blues")[2:5]) +
  scale_y_continuous(expand = c(0, 0),limits=c(0,28))
dev.off()

# do same just for linked NUMT-FPs
m=m[m$linked.NUMT.FP==1,]
tmp1=aggregate(samples$count, by=list(mtcnbinlabel=samples$mtcnbinlabel,insertsizebinlabel=samples$insertsizebinlabel), FUN=sum)
colnames(tmp1)[3]="n"
tmp2=aggregate(m$count, by=list(mtcnbinlabel=m$mtcnbinlabel,insertsizebinlabel=m$insertsizebinlabel), FUN=sum)
colnames(tmp2)[3]="NUMTFP.hets1to50"
tmp3=merge(tmp2,tmp1,by=c("mtcnbinlabel","insertsizebinlabel"))
tmp3$meanNUMTFP.hets1to50=tmp3$NUMTFP.hets1to50/tmp3$n
#exclude mtCN>300
pdf("plots/FigS3B.pdf", width=5, height=4,useDingbats=FALSE)
ggplot(tmp3,aes(y=meanNUMTFP.hets1to50,x=mtcnbinlabel,color=insertsizebinlabel)) +
  geom_point(size=1) +
  ggtitle("Mean # NUMT-FP/sample by mtCN and insert-size") +
  xlab("mtDNA copy number (mtCN)") +
  ylab("Mean # NUMT-FP/sample\n(heteroplasmy 1-50%)") +
  labs(color = "insert-size") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90)) +
  scale_color_manual(values=brewer.pal(n = 6, name = "Blues")[2:5]) +
  scale_y_continuous(expand = c(0, 0),limits=c(0,28))
dev.off()
