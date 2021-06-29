#extract G>A synonymous variants not seen at homoplasmy in gnomAD

dir.create("tables")

Ts_list <- c("A>G","C>T","G>A","T>C")
annotated_syn_vcf <- read.delim(file='annotated_synthetic.vcf', header=TRUE, sep = "\t", na.strings = "")
annotated_syn_vcf$variant <- paste("m.",annotated_syn_vcf$POS,annotated_syn_vcf$REF,">",annotated_syn_vcf$ALT,sep="")
annotated_syn_vcf$mut <- ifelse(paste(annotated_syn_vcf$REF,">",annotated_syn_vcf$ALT,sep="") %in% Ts_list,as.character(paste(annotated_syn_vcf$REF,">",annotated_syn_vcf$ALT,sep="")),"Tv")
annotated_syn_vcf$mitomap_ac <- (annotated_syn_vcf$Mitomap_af * 51836) #to convert to AC

write.table(annotated_syn_vcf[grepl("synonymous",annotated_syn_vcf$Consequence) & annotated_syn_vcf$mut=="G>A" & annotated_syn_vcf$gnomAD_max_hl<0.95,c("variant","SYMBOL","Consequence","HGVSc","HGVSp","gnomAD_max_hl","Helix_max_hl","mitomap_ac")],
            file = "tables/TableS3.txt", append = FALSE, sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)




