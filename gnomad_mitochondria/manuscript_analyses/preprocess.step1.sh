# preprocess step1: processes output data on the 56434 samples included in the v3.1 release 
# Specifically, this script:
# a) converts VCF file into tab delimited files for manuscript figure preparation
#    (summary.var.txt has one line per variant and all.sample2var.txt has one line per sample-variant)
# b) calculates MNVs 
# Note I have always run this with lots of memory using uger shell command: ish -l h_vmem=30g

##############################################################
# Dependencies: copies raw data from Kristen's directory 
#  Variable SOURCE_DIR is location of raw files exported from Hail by Kristen
#  These are produced by running combine_vcfs.py and add_annotations.py
#  $SOURCE_DIR/combined_sites_only_gnomad.vcf.bgz
#  $SOURCE_DIR/sample_annotations_gnomad.vcf.bgz
#  $SOURCE_DIR/sample_annotations_gnomad.txt
#
# Scripts in this directory:
#  ./vcf_to_table.nosplit.pl
#  ./merge_sets.pl
#  ./count_mnv.pl
#  ./preprocess.step2.sh
#
# OUTPUTS
#  sample_annotations_gnomad.txt : one line per sample
#  summary.var.txt : one line per variant
#  all.sample2var.txt : one line per sample-variant
#  pos.ref.alt.indel_stack.txt : one line per variant (that is an indel_stack)
#  summary.var.mnv.txt : summary data for calculating MNVs
#  mtdna_mnvs.txt : one line per output mnv
##############################################################


##############################################################
# get data from Kristen's directory and unzip it
##############################################################
cp $SOURCE_DIR/combined_sites_only_gnomad.vcf.bgz .
cp $SOURCE_DIR/sample_annotations_gnomad.vcf.bgz .
cp $SOURCE_DIR/sample_annotations_gnomad.txt .
mv combined_sites_only_gnomad.vcf.bgz combined_sites_only_gnomad.vcf.gz  
mv sample_annotations_gnomad.vcf.bgz sample_annotations_gnomad.vcf.gz
gunzip combined_sites_only_gnomad.vcf.gz
gunzip sample_annotations_gnomad.vcf.gz

##############################################################
# create a text summary file of variants from combined_sites_only_gnomad.vcf INFO fields:
# output columns: AN, AC_hom, AC_het, hap_defining_variant
##############################################################
grep -v "^\#\#" combined_sites_only_gnomad.vcf | perl -p -e 's/^(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t/$2\.$4\.$5\t$7\t/g' | perl -p -e 's/\t[^\t]*hap_defining_variant;.*AN=(\d+);AC_hom=(\d+);AC_het=(\d+);.*/\thap_defining_variant\t$1\t$2\t$3/g' | perl -p -e 's/\t[^\t]*AN=(\d+);AC_hom=(\d+);AC_het=(\d+);.*/\t0\t$1\t$2\t$3/g' | perl -p -e 's/INFO/hap_defining_variant\tAN\tAC_hom\tAC_het/g' > summary.var.txt

##############################################################
# convert vcf to all.sample2var.txt (one line per person-variant)
##############################################################
cat sample_annotations_gnomad.vcf | ./vcf_to_table.nosplit.pl -pass 6 -m 0.0001 > all.sample2var.txt

mv all.sample2var.txt all.sample2var.with_spaces.txt
# now convert to partipant_id
cut -f 1 all.sample2var.with_spaces.txt > t1
cut -f 2-5 all.sample2var.with_spaces.txt > t2
./merge_sets.pl -header -byname -x sample_annotations_gnomad.txt=SAMPLE_ID,s,participant_id t1 | cut -f 2 > t1.remap
paste t1.remap t2 > all.sample2var.txt

# get indel stacks
echo "INDEL.STACK" > pos.ref.alt.indel_stack.txt
grep indel summary.var.txt | cut -f 1 >> pos.ref.alt.indel_stack.txt

############################################################
# find MNVs: multi-nucleotide variants
# define: homoplasmic variants that are present in the same
# sample with an adjacent homoplasmic variants in at least 
# 90% of the time this variant is observed
############################################################
# sort all homoplasmies by sample then position
head -n 1 all.sample2var.txt | perl -p -e 's/DP/DP\tPOS/g' > all.sample2var.homs.txt
cat all.sample2var.txt | awk -F "\t" ' $4 >= 0.95 { print } ' | perl -p -e 's/^(\S+)\t(\S+)\.(\S+)\.(\S+)\t(.*)$/$1\t$2\.$3\.$4\t$5\t$2/g' | sort -k1,1 -k2,6n >> all.sample2var.homs.txt

# get MNV count
cat all.sample2var.homs.txt | ./count_mnv.pl > summary.var.mnv.txt
# flag as MNV any variants that are MNV in >=90% of samples with the given variant
cat summary.var.mnv.txt | awk -F "\t" ' $4 > 0.90 { print } ' > mtdna_mnvs.txt

