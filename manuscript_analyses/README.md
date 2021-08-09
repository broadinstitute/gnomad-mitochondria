This directory contains scripts and data files needed to generate figures and analyses included in the gnomAD mtDNA manuscript.

############################################################################
# Step 1
# Run preprocess scripts (this will take several hours and require high memory)
# Note: 
#  preprocess.step1.sh generates processed data on the 56434 samples included in the release; 
#  preprocess.step2.sh generates processed data on the 70K samples before filtering
#############################################################################
ish -l h_vmem=30g (optional but recommended: will launch interactive shell with 30G memory)
sh preprocess.step1.sh
sh preprocess.step2.sh

#############################################################################
# Note these preprocess bash shell scripts call perl scripts for further data processing
# Run script with -h option for full details.
#
# vcf_to_table.nosplit.pl : perl script inputs a VCF file and outputs a tab delimited 
# output file with one line per variant per sample (with parameters specifying minimum 
# thresholds for output)
#
# merge_sets.pl : perl script performs a database join (similar to Excel vlookup function) 
# between two tab-delimited files that share a common key column
#
# count_mnv.pl : perl script to identify multi-nucleotide variants (MNVs);  
# Input is tab-delimited file (one line per variant per person) sorted by SAMPLE_ID then by POS;  
# Output is set of unique variants (POS.REF.ALT) with the allele count of that variant in the 
# dataset (AC), the allele count of that variant as part of a multi-nucleotide variant (AC_MNV) 
# and the percentage (AC_MNV/AC).  
# For example: adjacent variants 5185.G.A (present 1 sample) and 5186.A.T (present in 80 samples) 
# were observed together in 1 sample and thus two lines would be 
# output “5185.G.A 1 1 1” and “5186.A.T 80 1 0.0125”
#############################################################################

#############################################################################
# Step 2
# Call R scripts that generate statistics and figures (see plots/ directory); 
# note fig2.R takes >30min
#############################################################################
use R-3.5
mkdir -p plots
Rscript --vanilla --default-packages="utils,grDevices,graphics,stats" generate_fig1.R  
Rscript --vanilla --default-packages="utils,grDevices,graphics,stats" generate_fig2.R
Rscript --vanilla --default-packages="utils,grDevices,graphics,stats" generate_figS2.R
Rscript --vanilla --default-packages="utils,grDevices,graphics,stats" generate_figS3.R
Rscript --vanilla --default-packages="utils,grDevices,graphics,stats" generate_fig3E.R
