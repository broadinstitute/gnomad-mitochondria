library(stringr)

dir.create("tables")

gnomad <- read.delim(file = 'reformated.vcf', header = TRUE, sep = "\t", na.strings = "")

gnomad$variant <- paste("m.", gnomad$POS, gnomad$REF, ">", gnomad$ALT, sep = "")
gnomad$HGVSp <- str_split(gnomad$HGVSp, ":", simplify = TRUE)[, 2]

write.table(gnomad[((gnomad$VARIANT_CLASS != "SNV" & grepl('frameshift', gnomad$Consequence)) | gnomad$Consequence == "stop_gained") & gnomad$max_hl >= 0.95, c("variant", "SYMBOL", "HGVSp", "Consequence", "max_hl", "AC_hom", "AC_het")],
            file = "tables/TableS3.txt", append = FALSE, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


