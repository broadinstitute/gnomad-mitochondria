#command for vep file output
#input is vcf produced by make_syntheticvcf_mtDNA.sh
#Using vep installed locally
ensembl-vep/vep --cache -i synthetic_vcf/NC_012920.1_synthetic.vcf -o synthetic_vcf/NC_012920.1_synthetic_vep.vcf --distance 0 --biotype --symbol --hgvs --variant_class --force_overwrite --vcf

#synthetic_vcf/split_vars_two_genes.py is used to separate variants within two genes, with two variant consequences per VEP

#used python script check_synthetic_vcf.py to create the same synthetic vcf, as a check of make_syntheticvcf_mtDNA.sh
diff -y --suppress-common-lines NC_012920.1_synthetic.vcf alternate_synthetic.vcf
#this showed synthetic vcf files are the same, for both scripts
