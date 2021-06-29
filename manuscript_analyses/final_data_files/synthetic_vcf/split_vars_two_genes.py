#!/usr/bin/python

import csv
import sys

vep_vcf = sys.argv[1]  # synthetic_vcf/NC_012920.1_synthetic_vep.vcf


def split(file1):

    with open(file1) as tsv_file:
        vep_vcf = csv.reader(tsv_file, delimiter="\t")

        file_name = file1.split(".vcf")[0]

        file_output_vep = open("%s_splitvarstwogenes.vcf" % file_name, "w")
        header_vep = "CHROM	POS	ID	REF	ALT	QUAL	FILTER	Allele	Consequence	IMPACT	SYMBOL	Gene	Feature_type	Feature	BIOTYPE	EXON	INTRON	HGVSc	HGVSp	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	DISTANCE	STRAND	FLAGS	VARIANT_CLASS	SYMBOL_SOURCE	HGNC_ID	HGVS_OFFSET"
        file_output_vep.write(header_vep + "\n")

        for row in vep_vcf:
            if not row[0].startswith("#"):
                variant_details = "\t".join(row[0:7])
                info = row[7]

                vep_append = info.replace("|", "\t")
                # for variants in two genes, split on ,
                if "," in vep_append:
                    vep_append1 = vep_append.split(",")[0]
                    vep_append2 = vep_append.split(",")[1]
                    file_output_vep.write(variant_details + "\t" + vep_append1 + "\n")
                    file_output_vep.write(variant_details + "\t" + vep_append2 + "\n")

                else:
                    file_output_vep.write(variant_details + "\t" + vep_append + "\n")


split(vep_vcf)
