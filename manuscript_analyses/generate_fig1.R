#!/broad/software/free/Linux/redhat_7_x86_64/pkgs/r_3.5.0-bioconductor/bin/Rscript
############################
# This script generates fig 1 panels and related statistics
#
# Input files:
#  all.sample2var.unfiltered.txt: unfiltered variants
#
# Output files:
#  output.stats.fig1.txt
#
# Figure list (see plots/ directory)
# Fig1B: nDNA_cov x mtDNA_cov by cohort (3 cohorts colored)
# Fig1C: mtDNA_cov hist (3 selected cohorts colored) 
# Fig1D: nDNA_cov hist (3 selected cohorts colored) 
# Fig1E: mtCN hist (3 selected cohorts colored) 
############################
library(ggplot2)

############################
# input files: sample annotation file
############################
samples=read.delim("unfiltered_sample_annotations.annot.txt",stringsAsFactors=FALSE)
samples$count=1

############################
# Statistics output file
############################
outf="output.stats.fig1.txt"
line="Statistics for gnomAD related to fig1"
write(line,file=outf)

# get 3 selected cohorts
# note this mistakenly excludes 67 1KG/HGDP samples that have prefix v3.1::HG or v3.1::NA (only 4 of which pass sample filters);  however this does not substantially change results
samples$pref=substr(samples$participant_id,1,2)
samples[samples$pref %in% c("NA","HG"),"research_project"]="1KG/HGDP"
samples[is.na(samples$research_project),"research_project"]=0
samples$selectedcohort=0
samples[samples$research_project=="NHLBI_WholeGenome_Sequencing","selectedcohort"]="NHLBI"
samples[samples$research_project=="1KG/HGDP","selectedcohort"]="1KG/HGDP"
samples[samples$research_project=="TOPMED_COPD","selectedcohort"]="TOPMED_COPD"

# print statistics per selected cohort: # Samples, median mtCN
line=paste("\nCohort statistics:\nCOHORT NUM_SAMPLES median_mtCN")
write(line,file=outf,append=TRUE)
for (mycohort in c("NHLBI","TOPMED_COPD","1KG/HGDP")) {
 line=paste(mycohort,"\t",sum(samples$selectedcohort==mycohort),"\t",median(samples[samples$selectedcohort==mycohort,"mtCN"]),sep="")
 write(line,file=outf,append=TRUE)
} 

# create subset dataframe with just the columns of interest
# for overlaid histogram, dataframe will concat all samples to selected cohorts
tmp1=samples[samples$selectedcohort != "0",c("selectedcohort","wgs_median_coverage","mt_mean_coverage","mtCN")]
colnames(tmp1)[1]="cohort"
tmp2=samples[,c("count","wgs_median_coverage","mt_mean_coverage","mtCN")]
colnames(tmp2)[1]="cohort"
tmp2$cohort="all"
tmp=rbind(tmp2,tmp1)
tmp$cohort=factor(tmp$cohort,levels=c("all","NHLBI","TOPMED_COPD","1KG/HGDP"))
mycolors=c("gray","peachpuff","sandybrown","salmon4")

############################
# Fig1C: mtDNA_cov hist (3 selected cohorts colored) 
############################
maxX=15000
pdf("plots/Fig1C.pdf", width=6, height=4)
ggplot(tmp, aes(x=mt_mean_coverage,color=cohort,fill=cohort)) +
  geom_histogram(position="identity",binwidth=25) +
  theme_classic() + ylab("# Samples") +
  scale_x_continuous(name="mtDNA mean coverage", limits=c(0,maxX),expand=c(0,0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values=mycolors) +
  scale_color_manual(values=mycolors)
dev.off()
line=paste("\nCoverage histograms: missing data not shown:\nmtDNA coverage: ",sum(tmp$mt_mean_coverage > maxX)," samples with coverage ",maxX,"-",max(tmp$mt_mean_coverage),sep="")
write(line,file=outf,append=TRUE)
############################
# Fig1D: nDNA_cov hist (3 selected cohorts colored) 
############################
minX=10
maxX=60
pdf("plots/Fig1D.pdf", width=6, height=4)
ggplot(tmp, aes(x=wgs_median_coverage,color=cohort,fill=cohort)) +
  geom_histogram(position="identity",binwidth=1) +
  theme_classic() +
  ylab("# Samples") +
  scale_x_continuous(name="nDNA median coverage", limits=c(minX,maxX),expand=c(0,0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values=mycolors) +
  scale_color_manual(values=mycolors)
dev.off()
line=paste("nDNA coverage: ",sum(tmp$wgs_median_coverage > maxX)," samples with coverage ",maxX,"-",max(tmp$wgs_median_coverage),sep="")
write(line,file=outf,append=TRUE)
############################
# Fig1E: mtCN hist (3 selected cohorts colored)
############################
maxX=1250
pdf("plots/Fig1E.pdf", width=6, height=4)
ggplot(tmp, aes(x=mtCN,color=cohort,fill=cohort)) +
  geom_histogram(position="identity",binwidth=5) +
  theme_classic() +
  ylab("# Samples") +
  scale_x_continuous(name="mtDNA copy number (mtCN)", limits=c(0,maxX),breaks=c(250,500,750,1000,1250),expand=c(0,0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values=mycolors) +
  scale_color_manual(values=mycolors) +
  geom_vline(xintercept=50, linetype="dashed", color = "black") +
  geom_vline(xintercept=500, linetype="dashed", color = "black")
dev.off()
line=paste("mtCN hist: ",sum(tmp$mtCN > maxX)," samples with mtCN ",maxX,"-",max(tmp$mtCN),sep="")
write(line,file=outf,append=TRUE)

##########################
# Fig1B: nDNA_cov x mtDNA_cov by cohort (3 cohorts colored) file=plots/scatterplot.nDNA_vs_mtDNA.cohorts.pdf
##########################                                       
# step 1: get matrix by cohort with mean/sd for nDNA cov and mean/sd for mtDNA cov
# plus factors for cohort subset and size
tmp1=aggregate(samples$count, by=list(cohort=samples$research_project), FUN=sum)
colnames(tmp1)[2]="n"
tmp2=aggregate(samples$wgs_median_coverage, by=list(cohort=samples$research_project), FUN=mean)
colnames(tmp2)[2]="nDNA.mean"
tmp1=merge(tmp1,tmp2,by=c("cohort"))
tmp2=aggregate(samples$wgs_median_coverage, by=list(cohort=samples$research_project), FUN=sd)
colnames(tmp2)[2]="nDNA.sd"
tmp1=merge(tmp1,tmp2,by=c("cohort"))
tmp2=aggregate(samples$mt_mean_coverage, by=list(cohort=samples$research_project), FUN=mean)
colnames(tmp2)[2]="mtDNA.mean"
tmp1=merge(tmp1,tmp2,by=c("cohort"))
tmp2=aggregate(samples$mt_mean_coverage, by=list(cohort=samples$research_project), FUN=sd)
colnames(tmp2)[2]="mtDNA.sd"
tmp1=merge(tmp1,tmp2,by=c("cohort"))
# get sd min/max
tmp1$ndnamin=tmp1$nDNA.mean-tmp1$nDNA.sd
tmp1$ndnamax=tmp1$nDNA.mean+tmp1$nDNA.sd
tmp1$mtdnamin=tmp1$mtDNA.mean-tmp1$mtDNA.sd
tmp1$mtdnamax=tmp1$mtDNA.mean+tmp1$mtDNA.sd
# get NumSamples by bins (for display)
tmp1$NumSamples="10000+"
tmp1[tmp1$n < 10000,"NumSamples"]="1000-10000"
tmp1[tmp1$n < 1000,"NumSamples"]="<1000"
tmp1$NumSamples=factor(tmp1$NumSamples,levels=c("<1000","1000-10000","10000+"))
# get CohortSet so we can color 3 selected cohorts
tmp1$CohortSet="all"
tmp1[tmp1$cohort=="1KG/HGDP","CohortSet"]="1KG/HGDP"
tmp1[tmp1$cohort=="NHLBI_WholeGenome_Sequencing","CohortSet"]="NHLBI"
tmp1[tmp1$cohort=="TOPMED_COPD","CohortSet"]="TOPMED_COPD"
tmp1$CohortSet=factor(tmp1$CohortSet,levels=c("all","NHLBI","TOPMED_COPD","1KG/HGDP"))

pdf("plots/Fig1B.pdf", width=6, height=4,useDingbats=FALSE)
ggplot(data = tmp1,aes(x = nDNA.mean,y = mtDNA.mean,color=CohortSet)) +
  geom_errorbar(aes(ymin = mtdnamin,ymax = mtdnamax),color="gray90",size=0.5) +
  geom_errorbarh(aes(xmin = ndnamin,xmax = ndnamax),color="gray90",size=0.5) +
  theme_classic() +
  geom_point(aes(size=NumSamples)) +
  geom_text(aes(label=ifelse(CohortSet != "all",as.character(CohortSet),'')) ,vjust=-1,size=2) +
  scale_color_manual(values=mycolors) +
  coord_cartesian(ylim = c(0, 20000)) +
  xlab("mean nDNA coverage") +
  ylab("mean mtDNA coverage")
dev.off()

