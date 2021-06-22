#!/broad/software/free/Linux/redhat_7_x86_64/pkgs/r_3.5.0-bioconductor/bin/Rscript
############################
## This script generates figure 2, NUMT analyses, and related statistics
##
# Input files:
#  unfiltered_sample_annotations.annot.txt: annotations per sample
#  all.sample2var.unfiltered.txt: unfiltered variants
#  all.sample2var.passplus.txt: variants passing filters 
#  summary.unfiltered.var.vep.txt: VEP annotations of all unique variants
#  pos.ref.alt.indel_stack.txt: list of sites subsequently filtered out as indel stacks
#
# Output files:
#  output.stats.fig2.txt
#  TableS1.common_heteroplasmies.txt
#  TableS1.common_heteroplasmies.correlations.txt
#  all.sample2var.passplus.hets.annotNUMTFP.txt
#  all.sample2var.unfiltered.annotNUMTFP.txt
#
## Figure list (see plots/ directory)
## 2E: mtCN x # heteroplasmies/sample (red NUMT-FP) file plots/stackedbar.meanhets.by.mtCNbin.pdf
## 2F: 16293 heteroplasmy vs theoretical [file=plots/16293.A.C.observed.v.expected.pdf]
## 2G: heteroplasmy x % NUMT-FP (color by mtCN) file=plots/vafbin.vs.percNUMTFP.vs.mtCNbin.pdf
## 2H: density plot known cell lines vs other [file=plots/cell.density.nsamples.pdf]
## 2I: #mutations/sample in known cell lines vs mtCN50-500, group by VEP impact [file=plots/cell.stackedbar.meanhets.by.cellline.varType.pdf]
## 2J: overlaid histogram of VAF bin x % NUMT-FP (red for all samples, orange post-filtering [file=plots/hist.perc_numtfp]
## 2K: release samples only: heteroplasmy bin x # variants/sample (red NUMT-FP) [file=plots/stackedbar.release.hl_v_nvar.noHom.pdf]
############################
library(ggplot2)
library(reshape2)
library(RColorBrewer)

############################
# input files
############################
allm=read.delim("all.sample2var.unfiltered.txt",stringsAsFactors=FALSE)
pass=read.delim("all.sample2var.passplus.txt",stringsAsFactors=FALSE)
allm$count=1
pass$count=1

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

############################
# Statistics output file
############################
outf="output.stats.fig2.txt"
line="Statistics for gnomAD related to fig2 and pipeline"
write(line,file=outf)


############################
# identify candidate NUMTs, and linked NUMT-FPs from 2 validated NUMTs
#
# Output files:
#  TableS1.common_heteroplasmies.txt
#  TableS1.common_heteroplasmies.correlations.txt
#  all.sample2var.passplus.hets.annotNUMTFP.txt
#  all.sample2var.unfiltered.annotNUMTFP.txt
############################

getVar2Count <- function(mat,label) {
  var2count=aggregate(mat$count, by=list(POS.REF.ALT=mat$POS.REF.ALT), FUN=sum)
  colnames(var2count)[2]=label
  var2count
}

getVar2CountAndMerge <- function(mat,label,sofar) {
  tmp=getVar2Count(mat,label)
  merge(sofar,tmp,all=TRUE)
}

############################
# step 1: get common heteroplasmies: UNFILTERED and het 0-0.5 and min mincommon 
mincommon=1000
maxhet=0.5
myhets=allm[allm$HL < maxhet,]
tmp1=aggregate(myhets$count, by=list(POS.REF.ALT=myhets$POS.REF.ALT,numtID=myhets$numtID), FUN=sum)
colnames(tmp1)[3]="ALL_HETS_0to50"
common=tmp1[tmp1$ALL_HETS_0to50 > mincommon,]
ncommon=dim(common)[[1]]
# 122 common hets, with unfiltered variants 0-50% heteroplasmy present in 1000 samples
line=paste("\nStatistics for calculating candidate NUMTS\n",ncommon," candidate NUMTS with unfiltered variants (heteroplasmy 0-",maxhet,") detected in at least ",mincommon," samples",sep="")
write(line,file=outf,append=TRUE)

# add in frequency of pass variants 0-50 and pass variants 1-50
tmp=pass[(pass$POS.REF.ALT %in% common$POS.REF.ALT) & (pass$HL < maxhet),]
common=getVar2CountAndMerge(tmp,"PASS_HETS_0to50",common)
tmp=pass[(pass$POS.REF.ALT %in% common$POS.REF.ALT) & (pass$HL < maxhet) & (pass$HL >= 0.01),]
common=getVar2CountAndMerge(tmp,"PASS_HETS_1to50",common)

############################
# step 2: find common heteroplasmies that correlate with mtCN
allcommonhets=allm[(allm$POS.REF.ALT %in% common$POS.REF.ALT) & (allm$HL < maxhet),]
common$SPEARMAN=0
common$SPEARMAN.PVAL=1
mycor <- function(x,y,method) {
  cor(x=x,y=y,method=method,use="pairwise.complete.obs")
}
mycor.test <- function(x,y,method,alternative,exact) {
  cor.test(x=x,y=y,method=method,alternative=alternative,exact=exact,use="pairwise.complete.obs")
}
# test correlation in allcommonhets (not just PASS)
for (i in 1:ncommon) {
  myvar=common$POS.REF.ALT[i]
  tmp=allcommonhets[allcommonhets$POS.REF.ALT == myvar,]
  common[i,"SPEARMAN"]=mycor(x=tmp$HL,y=1/(1+tmp$mtCN),method="spearman")
  common[i,"SPEARMAN.PVAL"]=mycor.test(x=tmp$HL,y=1/(1+tmp$mtCN),method="spearman",alternative="two.sided",exact=FALSE)$p.value
}
common$SIGNIFICANT = (common$SPEARMAN > 0.45) & (common$SPEARMAN.PVAL * ncommon < 0.00001)
numts=common[common$SIGNIFICANT,"POS.REF.ALT"]
write.table(common,file="TableS1.common_heteroplasmies.txt",quote=FALSE,sep="\t",row.names=TRUE,col.names=NA)

line=paste(length(numts),"/",ncommon," common hets are significantly correlated with 1/(1+mtCN) (r> 0.45, bonferroni corrected p< 1e-5)",sep="")
write(line,file=outf,append=TRUE)
#67/122 common hets are significantly correlated with 1/(1+mtCN) (r> 0.45, bonferroni corrected p< 1e-5)

############################
# step 3: Identify subset that co-occur 
# Note -- need to use full matrix, and then NA homoplasmic variants (because cannot tell if they are shared)

# pull subset of numts, then convert from Long to Wide form using dcast from reshape2 package
numtm=allm[allm$POS.REF.ALT %in% numts,]
commonLongForm=numtm[,c("participant_id","POS.REF.ALT","HL")]
commonWideForm <- dcast(commonLongForm, POS.REF.ALT ~ participant_id)
commonWideForm[is.na(commonWideForm)] = 0
rownames(commonWideForm)=commonWideForm[,"POS.REF.ALT"]
commonWideForm=commonWideForm[,2:(dim(commonWideForm)[[2]])]
# Null out high heteroplasmy variants (cannot assess low het variant presence)
commonWideForm[commonWideForm > 0.50] = NA
#dim(commonWideForm)
#[1]    67 69947
# do not add in blank columns for all other samples with no heteroplasmies (unneeded since we capture the vast majority: 69947/7-359 samples)
# get correlations
varcor=round(cor(t(commonWideForm),use="pairwise.complete.obs"),3)
varcor[is.na(varcor)]=0
varcor[varcor == 1] = 0
write.table(varcor,file="TableS1.common_heteroplasmies.correlations.txt",quote=FALSE,sep="\t",row.names=TRUE,col.names=NA)

############################
# step 4: identify linked NUMT-FP variants from 2 validated numts (numtA and numtB)
# 2 heteroplasmies (0-50%) in same NUMT, AND, HL within 0.05 of median for that NUMT

mycols=colnames(allm)
# add new columns concatenating other columns to be used as IDs 
allm$SAMPLE.NUMT=paste(allm$participant_id,allm$numtID,sep=".")
allm$SAMPLE.VAR=paste(allm$participant_id,allm$POS.REF.ALT,sep=".")
pass$SAMPLE.VAR=paste(pass$participant_id,pass$POS.REF.ALT,sep=".")

numtfp=allm[(allm$numtID != "") & (allm$HL < 0.50),]
# aggregate by SAMPLE.NUMT
tmp1=aggregate(numtfp$count, by=list(SAMPLE.NUMT=numtfp$SAMPLE.NUMT), FUN=sum)
colnames(tmp1)[2]="unfiltered_NUMTFP"
tmp1=tmp1[tmp1$unfiltered_NUMTFP > 1,]
tmp2=aggregate(numtfp$HL, by=list(SAMPLE.NUMT=numtfp$SAMPLE.NUMT), FUN=median)
colnames(tmp2)[2]="medianHL"
# tmp3 is the set of variants linked to another NUMT-FP (but any heteroplasmy)
tmp3=allm[,c("SAMPLE.VAR","SAMPLE.NUMT","HL")]
tmp3=tmp3[(tmp3$SAMPLE.NUMT %in% tmp1$SAMPLE.NUMT) & (tmp3$HL < 0.5),]
# tmp4 now joins in medianHL so we can find outliers that are too far from the median HL
tmp4=merge(tmp3,tmp2,by="SAMPLE.NUMT")
tmp4$diff=abs(tmp4$HL - tmp4$medianHL)
# linkedNUMTFP is the set of sample-variants that are linked to another NUMT heteroplasmy (0-50%) AND within 0.05 VAF of median for that NUMT
linkedNUMTFP=tmp4[tmp4$diff < 0.05,"SAMPLE.VAR"]

allm$linked.NUMT.FP=0
allm[(allm$SAMPLE.VAR %in% linkedNUMTFP) & (allm$HL < 0.5),"linked.NUMT.FP"]=1

# now add linked.NUMT.FP to pass
pass$SAMPLE.VAR=paste(pass$participant_id,pass$POS.REF.ALT,sep=".")
pass$linked.NUMT.FP=0
pass[(pass$SAMPLE.VAR %in% linkedNUMTFP) & (pass$HL < 0.5),"linked.NUMT.FP"]=1

# pull out just relevant columns and write table, now annotated with linked.NUMT.FP
allm=allm[,c(mycols,"linked.NUMT.FP")]
write.table(allm,file="all.sample2var.unfiltered.annotNUMTFP.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
pass=pass[,c(mycols,"linked.NUMT.FP")]
write.table(pass,file="all.sample2var.passplus.annotNUMTFP.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)


############################
# Now generate figures related to NUMTs
############################

#pass=read.delim("all.sample2var.passplus.annotNUMTFP.txt",stringsAsFactors=FALSE)
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

############################
# stats for paper
# % of all pass hets 1-50% that are linked-FPs
############################
line="\nStats on heteroplasmic PASS variants that are NUMT-FP (70K samples):"
write(line,file=outf,append=TRUE)

# PASS variants 1-50%:
mycalc<-function(minthresh,maxthresh,release) {
  if (release) {
    val1=sum(passhets[(passhets$HL >= minthresh) & (passhets$HL < maxthresh) & (passhets$release3.1.1 == "true"),"linked.NUMT.FP"])
    val2=length(passhets[(passhets$HL >= minthresh) & (passhets$HL < maxthresh) & (passhets$release3.1.1 == "true"),"linked.NUMT.FP"])
  } else {
    val1=sum(passhets[(passhets$HL >= minthresh) & (passhets$HL < maxthresh),"linked.NUMT.FP"])
    val2=length(passhets[(passhets$HL >= minthresh) & (passhets$HL < maxthresh),"linked.NUMT.FP"])
  }
  val3=round(val1/val2,2)
  if (release) {
    line=paste(val1,"/",val2,"(",val3,") PASS variants IN RELEASE 3.1.1 (heteroplasmy ",minthresh,"-",maxthresh,") are linked.NUMT.FP",sep="")
  } else {
    line=paste(val1,"/",val2,"(",val3,") PASS variants (heteroplasmy ",minthresh,"-",maxthresh,") are linked.NUMT.FP",sep="")
  }
 write(line,file=outf,append=TRUE)
}
mycalc(0.01,0.5,0)
mycalc(0.01,0.10,0)
mycalc(0.01,0.05,0)
# do all bins 0.05-0.10; 0.10-0.15, etc to 0.90-0.95
for (i in 1:18) {
  mycalc(i*0.05,i*0.05 + 0.05,0)
}

# do the same but only for the release samples (remove mtCN and contaminants)
line="\n\nStats on heteroplasmic PASS variants IN RELEASE SAMPLE that are NUMT-FP (56434 samples):"
write(line,file=outf,append=TRUE)
mycalc(0.01,0.5,1)
mycalc(0.01,0.10,1)
mycalc(0.01,0.05,1)
# do all bins 0.05-0.10; 0.10-0.15, etc to 0.90-0.95
for (i in 1:18) {
  mycalc(i*0.05,i*0.05 + 0.05,1)
}

############################
# Fig2D: mtCN x # heteroplasmies/sample (red NUMT-FP) 
############################
m=passhets
m=m[(m$HL < 0.50) & (m$HL >= 0.01),]
# get mean # heteroplasmies/sample per mtcn bin
tmp1=aggregate(samples$count, by=list(mtcnbinlabel=samples$mtcnbinlabel), FUN=sum)
colnames(tmp1)[2]="n"
tmp2=aggregate(m$count, by=list(mtcnbinlabel=m$mtcnbinlabel,linked.NUMT.FP=m$linked.NUMT.FP), FUN=sum)
colnames(tmp2)[3]="n.hets1to50"
tmp3=merge(tmp2,tmp1,by=c("mtcnbinlabel"))
tmp3$mean.hets1to50=tmp3$n.hets1to50/tmp3$n
tmp3$linked.NUMT.FP=factor(tmp3$linked.NUMT.FP,levels=c(1,0))
# exclude samples mtCN >300
tmp3=tmp3[!(tmp3$mtcnbinlabel %in% c("300-325","325-350","350-375","375-400","400-425","425-450","450-475","475-500","500+")),]
pdf("plots/Fig2D.pdf", width=5, height=3)
ggplot(tmp3, aes(fill=linked.NUMT.FP, y=mean.hets1to50, x=mtcnbinlabel)) + geom_bar(position="stack", stat="identity")+theme_classic() + theme(axis.text.x = element_text(angle=90))+ labs(fill = "NUMT-FP")+ xlab("mtDNA copy number (mtCN)") + ylab("Mean # variants/sample\n(heteroplasmy 1-50%)")+ scale_fill_manual(values=c("red","black"))
dev.off()

############################
# Fig2F: heteroplasmy x % NUMT-FP (color by mtCN)
############################
m=passhets
m=m[(m$HL < 0.1599) & (m$HL >= 0.01),]
m=m[!is.na(m$mtCN),]
m$vafbin=floor(m$HL*100)/100
m$mtcnbin=floor(m$mtCN / 25)*25
m[m$mtcnbin > 500,"mtcnbin"]=500
m$mtcnbinlabel=paste(m$mtcnbin,m$mtcnbin+25,sep="-")
m[m$mtcnbinlabel=="500-525","mtcnbinlabel"]="500+"
m$mtcnbinlabel=factor(m$mtcnbinlabel,levels=unique(m$mtcnbinlabel)[order(unique(m$mtcnbin))])
m$mtcnbin=factor(m$mtcnbin,levels=sort(unique(m$mtcnbin)))
# aggregate by vaf and mtCN
tmp1=aggregate(m$count, by=list(vafbin=m$vafbin,mtcnbinlabel=m$mtcnbinlabel), FUN=sum)
colnames(tmp1)[3]="n"
mnumt=m[(m$linked.NUMT.FP ==1),]
tmp2=aggregate(mnumt$count, by=list(vafbin=mnumt$vafbin,mtcnbinlabel=mnumt$mtcnbinlabel), FUN=sum)
colnames(tmp2)[3]="n.linkedNUMTFP"
tmp3=merge(tmp1,tmp2,all=TRUE,by=c("vafbin","mtcnbinlabel"))
tmp3$percNUMTFP=tmp3$n.linkedNUMTFP/tmp3$n
tmp3=tmp3[!is.na(tmp3$percNUMTFP),]
# exclude mtCN 250-500 for readability (and because these bins have few samples)
tmp3=tmp3[!(tmp3$mtcnbinlabel %in% c("250-275","275-300","300-325","325-350","350-375","375-400","400-425","425-450","450-475","475-500")),]
xlabs=paste(c(1:15)/100,c(2:16)/100,sep="-")
pdf("plots/Fig2F.pdf", width=4, height=4)
ggplot(tmp3,aes(y=percNUMTFP,x=vafbin,color=mtcnbinlabel))+geom_line() + geom_point(size=1) + ggtitle("% NUMT by VAF and mtCN")+ geom_vline(xintercept=0.05, linetype="dashed", color = "black")+ geom_vline(xintercept=0.15, linetype="dashed", color = "black")+ geom_vline(xintercept=0.1, linetype="dashed", color = "black") + scale_x_continuous(name="Heteroplasmy bin", limits=c(0.01, 0.15),breaks=c(c(1:15)/100),labels=xlabs) + scale_y_continuous(name="% NUMT-FP")+ labs(color = "mtCN")+theme_classic() + theme(axis.text.x = element_text(angle=90))+ scale_color_brewer(palette="Spectral")
dev.off()

########################################################
# Fig2E: 16293 heteroplasmy vs theoretical [file=plots/16293.A.C.observed.v.expected.pdf]
########################################################
myhet="16293.A.C"
 subset=passhets[(passhets$POS.REF.ALT == myhet) & (passhets$HL < 0.50) & (passhets$HL >= 0.01),]
 subset$expected=1/(1+subset$mtCN)
 subset$mtcncolor=subset$mtCN
 subset[subset$mtcncolor > 100,"mtcncolor"]=100
 n=dim(subset)[[1]]
 pdf(paste("plots/Fig2E.pdf",sep=""), width=2, height=2,useDingbats=FALSE)
 ggplot(subset,aes(y=HL,x=expected))+geom_point(size=0.5,color="red") + theme_classic()+xlab("Expected VAF for heterozygous NUMT\n1/(1+mtCN)")+ylab("Observed VAF")+ scale_x_continuous(expand = c(0, 0),limits=c(0,0.25),breaks=c(0.1,0.2)) + scale_y_continuous(expand = c(0, 0),limits=c(0,0.25),breaks=c(0.1,0.2))
 dev.off()
line=paste("\nStatistics for NUMT-FP m.16293A>C:\nn=",n,"; Spearman r=",signif(cor(x=subset$HL,y=subset$expected,use="pairwise.complete.obs")[1],2),"; p=",cor.test(x=subset$HL,y=subset$expected,method="spearman",alternative="two.sided",exact=FALSE)$p.value,
  "; Note p=0 indicates p<2.2e-16",sep="")
write(line,file=outf,append=TRUE)


########################################################
# Fig2L: release samples only: heteroplasmy bin x # variants/sample (red NUMT-FP) [file=plots/stackedbar.release.hl_v_nvar.noHom.pdf]
########################################################
#pass=read.delim("all.sample2var.passplus.annotNUMTFP.txt",stringsAsFactors=FALSE)

# lets look only at release samples, and only HL>=0.01
m=pass
m=m[(m$HL >=0.01) & (m$release3.1.1=="true"),]
m[m$HL ==1,"HL"]=0.999
m$vafbin=floor(m$HL*20)/20
m$vafbinlabel=paste(m$vafbin,m$vafbin+0.05,sep="-")
m$vafbinlabel=factor(m$vafbinlabel,levels=unique(m$vafbinlabel)[order(unique(m$vafbin))])
m$vafbin=factor(m$vafbin,levels=sort(unique(m$vafbin)))

# get mean # variants / heteroplasmy bin
tmp1=aggregate(m$count, by=list(vafbinlabel=m$vafbinlabel,linked.NUMT.FP=m$linked.NUMT.FP), FUN=sum)
colnames(tmp1)[3]="n"
tmp1$linked.NUMT.FP=factor(tmp1$linked.NUMT.FP,levels=c(1,0))
# exclude homoplasmic because it causes Y-axis to be too large
tmp1=tmp1[tmp1$vafbinlabel!="0.95-1",]
pdf("plots/Fig2L.pdf", width=5, height=3)
ggplot(tmp1, aes(fill=linked.NUMT.FP, y=n, x=vafbinlabel)) + geom_bar(position="stack", stat="identity")+theme_classic() + theme(axis.text.x = element_text(angle=90))+ labs(fill = "NUMT-FP")+ xlab("Heteroplasmy bin") + ylab("# variants")+ scale_fill_manual(values=c("red","black"))+ geom_vline(xintercept=2.5, linetype="dashed", color = "black")
dev.off()


########################################################
# Fig2K: overlaid histogram of VAF bin x % NUMT-FP (red for all samples, orange post-filtering
########################################################
#pass=read.delim("all.sample2var.passplus.annotNUMTFP.txt",stringsAsFactors=FALSE)
# lets look only at HL 0.01-0.50
m=pass
m=m[(m$HL >=0.01) & (m$HL<=0.50),]
m$vafbin=floor(m$HL*100)/100
m$vafbin=factor(m$vafbin,levels=sort(unique(m$vafbin)))

# get counts of numt-fp all data
subm=m[m$linked.NUMT.FP==0,]
tmp1=aggregate(subm$count, by=list(vafbin=subm$vafbin), FUN=sum)
colnames(tmp1)[2]="all.nonumt"
subm=m[m$linked.NUMT.FP==1,]
tmp2=aggregate(subm$count, by=list(vafbin=subm$vafbin), FUN=sum)
colnames(tmp2)[2]="all.numt"
tmp1=merge(tmp1,tmp2,all=TRUE)
subm=m[(m$linked.NUMT.FP==1) & (m$release3.1.1=="true"),]
tmp2=aggregate(subm$count, by=list(vafbin=subm$vafbin), FUN=sum)
colnames(tmp2)[2]="release.numt"
tmp1=merge(tmp1,tmp2,all=TRUE)
subm=m[(m$linked.NUMT.FP==0) & (m$release3.1.1=="true"),]
tmp2=aggregate(subm$count, by=list(vafbin=subm$vafbin), FUN=sum)
colnames(tmp2)[2]="release.nonumt"
tmp1=merge(tmp1,tmp2,all=TRUE)
tmp1[is.na(tmp1)]=0
tmp1$all.p.numt=tmp1$all.numt/(tmp1$all.numt+tmp1$all.nonumt)
tmp1$release.p.numt=tmp1$release.numt/(tmp1$release.numt+tmp1$release.nonumt)
tmp1$nonrelease.p.numt=tmp1$all.p.numt - tmp1$release.p.numt
tmp1[tmp1$nonrelease.p.numt<0,"nonrelease.p.numt"]=0

# for display, need to make long form of table
tmp2=tmp1[,c("vafbin","release.p.numt")]
tmp3=tmp1[,c("vafbin","nonrelease.p.numt")]
colnames(tmp2)[2]="p"
colnames(tmp3)[2]="p"
tmp2$label="release"
tmp3$label="filtered samples"
tmp4=rbind(tmp2,tmp3)
#tmp4$label=factor(tmp4$label,levels=c("release","filtered samples"))
tmp4$label=factor(tmp4$label,levels=c("filtered samples","release"))
tmp4$vafbin=as.numeric(as.character(tmp4$vafbin))
mybreaks=c(0,5,10,15,20,25,30,35,40,45,50)/100

pdf("plots/Fig2K.pdf", width=5, height=5)
ggplot(tmp4, aes(fill=label, y=p, x=vafbin)) + geom_bar(position="stack", stat="identity")+theme_classic() + scale_x_continuous(breaks=mybreaks, labels=mybreaks) + theme(axis.text.x = element_text(angle=90))+ xlab("VAF bin") + ylab("Fraction variants at 25 NUMT-FP sites")+ scale_fill_manual(values=c("red","orange"))+ geom_vline(xintercept=0.10, linetype="dashed", color = "black")+ geom_vline(xintercept=0.05, linetype="dashed", color = "gray")
dev.off()


############################
# Section D: figures related to cell lines
# Fig2G: density plot known cell lines vs other
# Fig2H: #mutations/sample in known cell lines vs mtCN50-500
############################
#pass=read.delim("all.sample2var.passplus.annotNUMTFP.txt",stringsAsFactors=FALSE)
#pass$count=1
#passhets=pass[pass$HL < 0.95,]

# assign cell line status based on participant id prefix in NA, HG, LP, JP
samples$pref=substr(samples$participant_id,1,2)
samples$is.cell="unknown"
samples[samples$pref %in% c("NA","HG","LP","JP","v3"),"is.cell"]="cell line"
pass$pref=substr(pass$participant_id,1,2)
pass$is.cell="unknown"
pass[pass$pref %in% c("NA","HG","LP","JP","v3"),"is.cell"]="cell line"
passhets=pass[pass$HL < 0.95,]
#stats
line=paste("\nCell Line Stats\n",sum(samples$is.cell== "cell line")," known cell lines",sep="")
write(line,file=outf,append=TRUE)
# 3523 known cell lines

############################
# Fig2G: density plot known cell lines vs other [file=plots/cell.density.nsamples.pdf]
############################
pdf("plots/Fig2G.pdf", width=4, height=2)
ggplot(samples, aes(x=mtCN, color=is.cell)) + geom_density()+theme_classic() + ylab("Density") + scale_x_continuous(name="mtDNA copy number (mtCN)", limits=c(0,2000)) + scale_color_manual(values=c("blue", "black"))
dev.off()
##############################
## 2H: #mutations/sample in known cell lines vs mtCN50-500, group by VEP impact [file=plots/cell.stackedbar.meanhets.by.cellline.varType.pdf]
##############################
# get variants heteroplasmy 10-95%
m=pass
m=m[(m$HL < 0.95)&(m$HL >= 0.10),]

# variant annotations: read in and classify varType as lof, nonsyn, syn, RNA, missense
vep=read.delim("summary.unfiltered.var.vep.txt",stringsAsFactors=FALSE)
m$varType="non-coding"
m[m$POS.REF.ALT %in% vep[vep$IMPACT %in% c("LOW"),"POS.REF.ALT"],"varType"]="synonymous"
m[m$POS.REF.ALT %in% vep[vep$Consequence == "missense_variant","POS.REF.ALT"],"varType"]="missense"
m[m$POS.REF.ALT %in% vep[vep$Consequence %in% c("stop_gained&frameshift_variant","stop_gained","frameshift_variant"),"POS.REF.ALT"],"varType"]="pLOF"
m[m$POS.REF.ALT %in% vep[vep$BIOTYPE %in% c("Mt_tRNA","Mt_rRNA"),"POS.REF.ALT"],"varType"]="rRNA/tRNA"
m$varType=factor(m$varType,levels=c("pLOF","missense","synonymous","rRNA/tRNA","non-coding"))
m$count=1

#cell2 is either known cell line ofr mtCN 50-500
samples$cell2=samples$is.cell
samples[(samples$mtCN >= 50) & (samples$mtCN <= 500) & (samples$is.cell == "unknown"),"cell2"]="mtCN 50-500"
samples$cell2=factor(samples$cell2,levels=c("unknown","mtCN 50-500","cell line"))
m$cell2=m$is.cell
m[(m$mtCN >= 50) & (m$mtCN <= 500) & (m$is.cell == "unknown"),"cell2"]="mtCN 50-500"
m$cell2=factor(m$cell2,levels=c("unknown","mtCN 50-500","cell line"))
#stats
line=paste(sum(samples$cell2=="cell line")," known cell line samples",sep="")
write(line,file=outf,append=TRUE)
line=paste(sum(samples$cell2=="mtCN 50-500")," mtCN 50-500 samples",sep="")
write(line,file=outf,append=TRUE)

# now get mean # variants of each varType in cell lines or mtCN50-500
tmp1=aggregate(samples$count, by=list(cell2=samples$cell2), FUN=sum)
colnames(tmp1)[2]="n"
tmp2=aggregate(m$count, by=list(cell2=m$cell2,varType=m$varType), FUN=sum)
colnames(tmp2)[3]="n.hets10to95"
tmp3=merge(tmp2,tmp1,by=c("cell2"))
tmp3$mean.hets10to95=tmp3$n.hets10to95/tmp3$n
tmp3$varType=factor(tmp3$varType,levels=c("pLOF","missense","rRNA/tRNA","synonymous","non-coding"))
tmp3=tmp3[tmp3$cell2 != "unknown",]
pdf("plots/Fig2H.pdf", width=5, height=3)
ggplot(tmp3, aes(fill=cell2, y=mean.hets10to95, x=varType)) + geom_bar(position="dodge", stat="identity")+theme_classic() + theme(axis.text.x = element_text(angle=90))+ xlab("mtDNA copy number (mtCN)") + ylab("Mean # variants/sample\n(heteroplasmy 10-95%)")+ scale_fill_manual(values=c("gray","blue"))
dev.off()

# print enrichment and significance (Fisher's Exact) for each level
line=paste("Cell line enrichment over mtCN50-500")
write(line,file=outf,append=TRUE)
for (mylevel in levels(tmp3$varType)) {
 #print(signif(fisher.test(x=tmp3[tmp3$varType==mylevel,c("n.hets10to95","n")])$p.value,2))
 line=paste(mylevel," enrichment=",signif(tmp3[(tmp3$varType == mylevel)&(tmp3$cell2=="cell line"),"mean.hets10to95"]/tmp3[(tmp3$varType == mylevel)&(tmp3$cell2=="mtCN 50-500"),"mean.hets10to95"],2)," p=",signif(fisher.test(x=tmp3[tmp3$varType==mylevel,c("n.hets10to95","n")])$p.value,2),sep="")
 write(line,file=outf,append=TRUE)
} 

##################
# Fig 2J: statistics breakdown of all unique variants into 6 mutually exclusive bins:
# note: exclude indel stack variants -- which are subsequently excluded by HAIL scripts
#                                     RELEASE-samples    FILTERED-samples-only
# VAF 0.95-1.00
# VAF 0.10-0.95 (excl 0.95-1.00)
# VAF 0.01-0.10 (excl 0.10-1.00)
##################
#pass=read.delim("all.sample2var.passplus.annotNUMTFP.txt",stringsAsFactors=FALSE)
#  pos.ref.alt.indel_stack.txt: list of sites subsequently filtered out as indel stacks
indelstack=read.delim("pos.ref.alt.indel_stack.txt",stringsAsFactors=FALSE)
indelstack=indelstack[,1]
m=pass[pass$HL>=0.01,]
m=m[!(m$POS.REF.ALT %in% indelstack),]
rel95=unique(m[(m$HL>=0.95) & (m$release3.1.1=="true"),"POS.REF.ALT"])
sofar=rel95
rel10=unique(m[(m$HL>=0.10) & (m$release3.1.1=="true") & !(m$POS.REF.ALT %in% sofar),"POS.REF.ALT"])
sofar=c(sofar,rel10)
rel01=unique(m[(m$HL>=0.01) & (m$release3.1.1=="true") & !(m$POS.REF.ALT %in% sofar),"POS.REF.ALT"])
sofar=c(sofar,rel01)
exc95=unique(m[(m$HL>=0.95) & (m$release3.1.1 != "true") & !(m$POS.REF.ALT %in% sofar),"POS.REF.ALT"])
sofar=c(sofar,exc95)
exc10=unique(m[(m$HL>=0.10) & (m$release3.1.1 != "true") & !(m$POS.REF.ALT %in% sofar),"POS.REF.ALT"])
sofar=c(sofar,exc10)
exc01=unique(m[(m$HL>=0.01) & (m$release3.1.1 != "true") & !(m$POS.REF.ALT %in% sofar),"POS.REF.ALT"])
sofar=c(sofar,exc01)

line=paste("\nBreakdown of all PASS unique variants into 6 mutually-exclusive categories:\n",
  length(rel95),"\t VAF 0.95-1.00 in RELEASE samples\n",
  length(rel10),"\t VAF 0.10-0.95 in RELEASE samples, excluding all above variants\n",
  length(rel01),"\t VAF 0.01-0.10 in RELEASE samples, excluding all above variants\n",
  length(exc95),"\t VAF 0.95-1.00 in FILTERED OUT samples, excluding all above variants\n",
  length(exc10),"\t VAF 0.10-0.95 in FILTERED OUT samples, excluding all above variants\n",
  length(exc01),"\t VAF 0.01-0.10 in RELEASE samples, excluding all above variants\n",
  length(sofar),"\t total variants (sum of 6 previous categories)\n",
  sep="")
write(line,file=outf,append=TRUE)


