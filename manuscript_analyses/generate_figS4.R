#!/broad/software/free/Linux/redhat_7_x86_64/pkgs/r_3.5.0-bioconductor/bin/Rscript
############################
# This script generates figure S4
#
# Input files:
#  all.sample2var.passplus.txt: variants passing filters
#  pos.ref.alt.indel_stack.txt: variants in indel stacks (excluded from analysis)
#
## Figure list (see plots/ directory);  All figures concern variants in 56434 PASS samples
## A: Number variant calls per VAFbin (0.01-0.95), colored by label {obs. homoplasmic, heteroplasmic-only}
## B: Fraction variant calls per VAFbin (0.01-1.00), colored by label {obs. homoplasmic, heteroplasmic-only}
## C: Fraction variant calls per VAFbin (0.01-1.00), colored by # variants/SAMPLE in this VAFbin -- excluding multiallelic
## D: Number variant calls per VAFbin (0.90-1.00), colored by label {obs. homoplasmic, heteroplasmic-only}
## E: Fraction variant calls per VAFbin (0.90-1.00), colored by Max Observed Heteroplasmy {<0.95, 0.95-0.99, 1.00}
## F: Fraction variant calls per VAFbin (0.90-1.00), colored by # variants/SAMPLE in this VAFbin -- excluding multiallelic
## G: Number variant calls per unique variant, colored by number distinct haplogroups
############################
library(ggplot2)

########################################################
# Variants in 56434 samples, including VAF 0.01-0.10
########################################################
m=read.delim("all.sample2var.passplus.txt",stringsAsFactors=FALSE)
m=m[(m$HL >=0.01) & (m$release3.1.1=="true"),]
m$count=1

indelstack=read.delim("pos.ref.alt.indel_stack.txt",stringsAsFactors=FALSE)
indelstack=indelstack[,1]
m=m[!(m$POS.REF.ALT %in% indelstack),]


# annotate with obshom label: obs. homoplasmic; heterplasmic only
obshom=unique(m[m$HL >= 0.95,"POS.REF.ALT"])
m$obshom="heteroplasmic only"
m[m$POS.REF.ALT %in% obshom,"obshom"]="obs. homoplasmic"
m$obshom=factor(m$obshom,levels=c("heteroplasmic only","obs. homoplasmic"))

# annotate with MaxObservedHeteroplasmy (MOH) in: <0.95, 0.95-0.99, 1.00
obs1.00=unique(m[m$HL == 1.00,"POS.REF.ALT"])
m$MOH="<0.95"
m[m$POS.REF.ALT %in% obshom,"MOH"]="0.95-0.99"
m[m$POS.REF.ALT %in% obs1.00,"MOH"]="1.00"
m$MOH=factor(m$MOH,levels=c("<0.95","0.95-0.99","1.00"))

# annotate VAF bin in bins of 5, make top bin 0.95-1.00 (instead of 1.00-1.05 which does not make sense)
m$vafbin=floor(m$HL*20)/20
m[m$vafbin==1,"vafbin"]=0.95
m$vafbinlabel=paste(m$vafbin,m$vafbin+0.05,sep="-")
m$vafbinlabel=factor(m$vafbinlabel,levels=unique(m$vafbinlabel)[order(unique(m$vafbin))])
m$vafbin=factor(m$vafbin,levels=sort(unique(m$vafbin)))

# annotate VAF bin in bins of 1 (eg 0.90, 0.91, 0.92...)
m$vaf1bin=floor(m$HL*100)/100
m$vaf1bin=factor(m$vaf1bin,levels=sort(unique(m$vaf1bin)))

# annotate which variant calls are multiallelic (so we can exclude these for certain analyses)
### step 1: convert POS.REF.ALT into SAMPLE.POS 
### getPos function splits a string POS.REF.ALT and returns just the numeric POS
getPos <-function(str) {
  as.numeric(strsplit(str,split="\\.")[[1]][1])
}
m$POS=lapply(m$POS.REF.ALT,getPos)
m$SAMPLE.POS=paste(m$participant_id,m$POS,sep=".")
### step 2: aggregate SAMPLE.POS and assign those with >1 count as "isMultiAllelic"
countMultiAlleles=aggregate(m$count, by=list(SAMPLE.POS=m$SAMPLE.POS), FUN=sum)
multiallelic=countMultiAlleles[countMultiAlleles$x>1,"SAMPLE.POS"]
m$isMultiAllelic=0
m[m$SAMPLE.POS %in% multiallelic ,"isMultiAllelic"]=1

# get subset of dataframe just with variants VAF 0.90-1.00
m90=m[m$HL >= 0.90,]

########################################################
# Fig S4A: number of variants per VAF bin
########################################################

# get total variants / heteroplasmy bin
vafbin2count=aggregate(m$count, by=list(vafbinlabel=m$vafbinlabel,obshom=m$obshom), FUN=sum)
colnames(vafbin2count)[3]="n"
# exclude homoplasmic because it causes Y-axis to be too large
tmp2=vafbin2count[vafbin2count$vafbinlabel!="0.95-1",]
pdf("plots/figS4A.pdf", width=5, height=3)
ggplot(tmp2, aes(fill=obshom, y=n, x=vafbinlabel)) +
  geom_bar(position="stack", stat="identity")+
  theme_classic() +
  theme(axis.text.x = element_text(angle=90))+
  labs(fill = "Label")+
  xlab("VAF bin") +
  ylab("# Variant calls")+
  scale_fill_manual(values=c("gray","black"))+
  geom_vline(xintercept=2.5, linetype="dashed", color = "black")
dev.off()

########################################################
# Fig S4B: fraction of variants per VAF bin
########################################################
pdf("plots/figS4B.pdf", width=5, height=3)
ggplot(vafbin2count, aes(fill=obshom, y=n, x=vafbinlabel)) +
  geom_bar(position="fill", stat="identity")+
  theme_classic() +
  theme(axis.text.x = element_text(angle=90))+
  labs(fill = "Label")+
  xlab("VAF bin") +
  ylab("Fraction variant calls")+
  scale_fill_manual(values=c("gray","black"))+
  geom_vline(xintercept=2.5, linetype="dashed", color = "black")
dev.off()

########################################################
# Fig S4C: Fraction variant calls per VAFbin (0.01-1.00),
# colored by # variants/SAMPLE in this VAFbin -- excluding multiallelic
########################################################

# get dataframe excluding multiallelic sites
mNoMultiAllelic=m[m$isMultiAllelic==0,]

# now count # heteroplasmies/sample by VAFbin
tmp=aggregate(mNoMultiAllelic$count, by=list(participant_id=mNoMultiAllelic$participant_id,vafbin=mNoMultiAllelic$vafbin), FUN=sum)
colnames(tmp)[3]="sample.hets.per.vafbin"
mat=merge(mNoMultiAllelic,tmp,all=TRUE)

# bin these into categories, max 3+
mat[mat$sample.hets.per.vafbin >= 3,"sample.hets.per.vafbin"]=3
mat[mat$sample.hets.per.vafbin == 3,"sample.hets.per.vafbin"]="3+"
mat$sample.hets.per.vafbin=factor(mat$sample.hets.per.vafbin,levels=c("1","2","3+"))

#now aggregate by categories
byVafBin=aggregate(mat$count, by=list(vafbinlabel=mat$vafbinlabel,sample.hets.per.vafbin=mat$sample.hets.per.vafbin), FUN=sum)
colnames(byVafBin)[3]="n"
pdf("plots/figS4C.pdf", width=5, height=3)
ggplot(byVafBin, aes(fill=sample.hets.per.vafbin, y=n, x=vafbinlabel)) +
  geom_bar(position="fill", stat="identity")+
  theme_classic() +
  theme(axis.text.x = element_text(angle=90))+
  labs(fill = "# variants/sample in this VAF bin")+
  xlab("VAF bin") +
  ylab("Fraction variant calls")+
  scale_fill_manual(values=c("steelblue1","steelblue3","steelblue4"))+
  geom_vline(xintercept=2.5, linetype="dashed", color = "black")
dev.off()

########################################################
# Fig S4D: number of variants per VAF bin (0.90-1.00)
########################################################

# get total variants / heteroplasmy bin
vafbin2count=aggregate(m90$count, by=list(vaf1bin=m90$vaf1bin,obshom=m90$obshom), FUN=sum)
colnames(vafbin2count)[3]="n"
pdf("plots/figS4D.pdf", width=5, height=3)
ggplot(vafbin2count, aes(fill=obshom, y=n, x=vaf1bin)) +
  geom_bar(position="stack", stat="identity")+
  theme_classic() +
  theme(axis.text.x = element_text(angle=90))+
  labs(fill = "Label")+
  xlab("VAF bin") +
  ylab("# Variant calls")+
  scale_fill_manual(values=c("gray","black"))+
  geom_vline(xintercept=5.5, linetype="dashed", color = "black")
dev.off()

########################################################
# Fig S4E: Fraction variant calls per VAFbin (0.90-1.00), colored by Max Observed Heteroplasmy {<0.95, 0.95-0.99, 1.00}
########################################################
vafbin2count=aggregate(m90$count, by=list(vaf1bin=m90$vaf1bin,MOH=m90$MOH), FUN=sum)
colnames(vafbin2count)[3]="n"
pdf("plots/figS4E.pdf", width=5, height=3)
ggplot(vafbin2count, aes(fill=MOH, y=n, x=vaf1bin)) +
  geom_bar(position="fill", stat="identity")+
  theme_classic() +
  theme(axis.text.x = element_text(angle=90))+
  labs(fill = "Max Obs. VAF")+
  xlab("VAF bin") +
  ylab("Fraction variant calls")+
  scale_fill_manual(values=c("gray","gray50","black"))+
  geom_vline(xintercept=5.5, linetype="dashed", color = "black")
dev.off()


########################################################
# Fig S4F: Fraction variant calls per VAFbin (0.90-1.00),
# colored by # variants/SAMPLE in this VAFbin -- excluding multiallelic
########################################################

# get dataframe excluding multiallelic sites
mNoMultiAllelic=m90[m90$isMultiAllelic==0,]

# now count # heteroplasmies/sample by VAFbin
tmp=aggregate(mNoMultiAllelic$count, by=list(participant_id=mNoMultiAllelic$participant_id,vaf1bin=mNoMultiAllelic$vaf1bin), FUN=sum)
colnames(tmp)[3]="sample.hets.per.vafbin"

mat=merge(mNoMultiAllelic,tmp,all=TRUE)

# bin these into categories, max 3+
mat[mat$sample.hets.per.vafbin >= 3,"sample.hets.per.vafbin"]=3
mat[mat$sample.hets.per.vafbin == 3,"sample.hets.per.vafbin"]="3+"
mat$sample.hets.per.vafbin=factor(mat$sample.hets.per.vafbin,levels=c("1","2","3+"))

#now aggregate by categories
byVaf1bin=aggregate(mat$count, by=list(vaf1bin=mat$vaf1bin,sample.hets.per.vafbin=mat$sample.hets.per.vafbin), FUN=sum)
colnames(byVaf1bin)[3]="n"
pdf("plots/figS4F.pdf", width=5, height=3)
ggplot(byVaf1bin, aes(fill=sample.hets.per.vafbin, y=n, x=vaf1bin)) +
  geom_bar(position="fill", stat="identity")+
  theme_classic() +
  theme(axis.text.x = element_text(angle=90))+
  labs(fill = "# variants/sample in this VAF bin")+
  xlab("VAF bin") +
  ylab("Fraction variant calls")+
  scale_fill_manual(values=c("steelblue1","steelblue3","steelblue4"))+
  geom_vline(xintercept=5.5, linetype="dashed", color = "black")
dev.off()

########################################################
# Fig S4G: # variant calls per unique variant, colored by number distinct haplogroups
########################################################

# restrict to VAF >=0.10
mat=m[m$HL >= 0.10,]

# annotate category cat as het-only, het-hom, or hom-only
obshet1=unique(m[(m$HL >=0.10 ) & (m$HL < 0.95),"POS.REF.ALT"])
obshom1=unique(m[m$HL >=0.95,"POS.REF.ALT"])
mat$cat="none"
mat[(mat$POS.REF.ALT %in% obshet1) & (mat$POS.REF.ALT %in% obshom1),"cat"]="het-hom"
mat[(mat$POS.REF.ALT %in% obshet1) & !(mat$POS.REF.ALT %in% obshom1),"cat"]="het-only"
mat[!(mat$POS.REF.ALT %in% obshet1) & (mat$POS.REF.ALT %in% obshom1),"cat"]="hom-only"
mat$cat=factor(mat$cat,levels=c("het-only","het-hom","hom-only"))

# annotate each sample's haplogroup2 (the first 2 letters of the haplogroup) into the matrix
sample=read.delim("unfiltered_sample_annotations.annot.txt",stringsAsFactors=FALSE)
sample$haplogroup2=substr(sample$major_haplogroup,0,2)
mat=merge(mat,sample[,c("participant_id","haplogroup2")])

# want to concentrate now on the heteroplasmies
hetsOnly=mat[mat$HL < 0.95,]

# get # hets per POS.REF.ALT
nhetsPerVar=aggregate(hetsOnly$count, by=list(POS.REF.ALT=hetsOnly$POS.REF.ALT), FUN=sum)
colnames(nhetsPerVar)[2]="nhetsPerVar"
hetsOnly=merge(hetsOnly,nhetsPerVar)

# now subset to just to 2+ hets per var
hets2Only=hetsOnly[hetsOnly$nhetsPerVar > 1,]

# get # per haplogroup2
byVarHap=aggregate(hets2Only$count, by=list(POS.REF.ALT=hets2Only$POS.REF.ALT,haplogroup2=hets2Only$haplogroup2), FUN=sum)
colnames(byVarHap)[3]="n"
byVarHap=byVarHap[order(byVarHap$POS.REF.ALT,byVarHap$n,decreasing=TRUE),]
byVarHap$count=1

# get # haplogroups per POS.REF.ALT
nhaplogroupsPerVar=aggregate(byVarHap$count, by=list(POS.REF.ALT=byVarHap$POS.REF.ALT), FUN=sum)
colnames(nhaplogroupsPerVar)[2]="nhaplogroupsPerVar"
hets2Only=merge(hets2Only,nhaplogroupsPerVar)

# now rank order variants by cat then by # haplogroups
byVar=aggregate(hets2Only$count, by=list(POS.REF.ALT=hets2Only$POS.REF.ALT,
                                   cat=hets2Only$cat,
                                   nhaplogroupsPerVar=hets2Only$nhaplogroupsPerVar), FUN=sum)
colnames(byVar)[4]="n"
byVar=byVar[order(byVar$cat,byVar$n,byVar$nhaplogroupsPerVar,decreasing=TRUE),]
byVar$POS.REF.ALT.rank=c(1:(dim(byVar)[[1]]))
hets2Only=merge(hets2Only,byVar[,c("POS.REF.ALT","POS.REF.ALT.rank")],all=TRUE)

# now max out nhaplogroupsPerVar at 3
hets2Only[hets2Only$nhaplogroupsPerVar >= 3 ,"nhaplogroupsPerVar"]="3+"
hets2Only$nhaplogroupsPerVar=factor(hets2Only$nhaplogroupsPerVar,levels=c("1","2","3+"))

hets2Only$count=1
byPos=aggregate(hets2Only$count,by=list(POS.REF.ALT=hets2Only$POS.REF.ALT,
                                  POS.REF.ALT.rank=hets2Only$POS.REF.ALT.rank,
                                  nhaplogroupsPerVar=hets2Only$nhaplogroupsPerVar),FUN=sum)
byPos$log2n=log(byPos$x,2)

pdf("plots/figS4G.pdf", width=22, height=6)
ggplot(byPos, aes(x=POS.REF.ALT.rank,y=x,fill=nhaplogroupsPerVar)) +
  geom_bar(stat="identity")+
  theme_classic() +
  xlab("Rank: unique variants with 2+ heteroplasmies, ranked by (i) het-hom vs het-only, (ii) # heteroplasmies/var, (iii) # haplogroups/var") +
  ylab("# Variant Calls") + scale_y_continuous(trans="log10")+ annotation_logticks()+
  scale_fill_manual(values=c("gray","goldenrod1","goldenrod3"))
dev.off()
