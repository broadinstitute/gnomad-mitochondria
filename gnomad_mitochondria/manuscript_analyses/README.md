# Manuscript analyses
This directory contains scripts and data files needed to generate figures and analyses included in the gnomAD mtDNA manuscript. Note that the gnomAD data files used (VCF with genotype and heteroplasmy level data for each sample, and txt with sample annotations) are not publicly available and in the code below refers to raw data present in the $SOURCE_DIR directory that is not publicly available. However these scripts detail how the raw data were processed and transformed into the manuscript figure panels and statistics. 


#######################################################################################################
## Notes for analyses and figures relating to NUMTs, coverage, insert sizes, and characteristics of heteroplasmic and homoplasmic variants (figures 1, 2, 3E, S2, S3, S4)
#######################################################################################################


Run preprocess scripts (this will take several hours and require high memory) that create tab delimited files required for R scripts
preprocess.step1.sh generates processed data on the 56,434 samples included in the release 
preprocess.step2.sh generates processed data on the 70K samples before filtering

ish -l h_vmem=30g (optional but recommended: will launch interactive shell with 30G memory)
```
sh preprocess.step1.sh
sh preprocess.step2.sh
```

The following preprocess bash shell scripts call several perl scripts for further data processing (run script with -h option for full details)
 vcf_to_table.nosplit.pl : perl script inputs a VCF file and outputs a tab delimited output file with one line per variant per sample (with parameters specifying minimum thresholds for output)
 merge_sets.pl : perl script performs a database join (similar to Excel vlookup function) between two tab-delimited files that share a common key column
 count_mnv.pl : perl script to identify multi-nucleotide variants (MNVs); Input is tab-delimited file (one line per variant per person) sorted by SAMPLE_ID then by POS; Output is set of unique variants (POS.REF.ALT) with the allele count of that variant in the dataset (AC), the allele count of that variant as part of a multi-nucleotide variant (AC_MNV) and the percentage (AC_MNV/AC). For example: adjacent variants 5185.G.A (present in 1 sample) and 5186.A.T (present in 80 samples) were observed together in 1 sample and thus two lines would be output “5185.G.A 1 1 1” and “5186.A.T 80 1 0.0125”


### Generate plots for figures 1, 2, S2, S3, 3E, S4

Call R scripts that generate statistics and figures (see plots/ directory); 
note fig2.R takes > 30min

```
use R-3.5
mkdir -p plots
Rscript --vanilla --default-packages="utils,grDevices,graphics,stats" generate_fig1.R
Rscript --vanilla --default-packages="utils,grDevices,graphics,stats" generate_fig2.R
Rscript --vanilla --default-packages="utils,grDevices,graphics,stats" generate_fig3E.R
Rscript --vanilla --default-packages="utils,grDevices,graphics,stats" generate_figS2.R
Rscript --vanilla --default-packages="utils,grDevices,graphics,stats" generate_figS3.R
Rscript --vanilla --default-packages="utils,grDevices,graphics,stats" generate_figS4.R
```



#######################################################################################################
## Notes for analyses and figures relating to variant statistics, sample distribution, and coverage post filtering
#######################################################################################################

`generate_fig1A_4_3BCD_S5ABC_S6.R` is used to plot variant statistics, sample distributions, and coverage metrics. This script takes four input files:
1. sample_annotations: The sample annotations file output by add_annotations.py for release samples
2. variant_data: Combined sites-only file (with txt extension) output by add_annotations.py
3. coverage: Mean coverages, tab-delimited file of 'locus', in format of chrM:position, and 'mean', mean coverage across all samples at that position (metrics can be obtained and extracted by running annotate_coverage.py)
4. plot_directory: Directory to which plots should be written

The output consists of plots for Figure 1A, Figure 4, Figure 3BCD, Supplemental Figure S5ABC, and Supplemental Figure S6.


#######################################################################################################
## Notes for analyses and figures relating to patterns of variation and pathogenic variation
#######################################################################################################

`reformat_vcf.py` is used to process and annotate the gnomAD data. The input files required for this script are located in `final_data_files`:
- In silico predictions, in directory `insilicos`: (a) APOGEE from MitImpact https://mitimpact.css-mendel.it/ and (b) HmtVar from https://www.hmtvar.uniba.it/. Note that you will need to unzip the HmtVar annotation file before using.
- Other database files, in directory `other_databases`: (a) MITOMAP disease and polymorphism variant lists, downloaded from https://www.mitomap.org/MITOMAP/resources, and (b) HelixMTdb dataset, downloaded from https://www.helix.com/pages/mitochondrial-variant-database.
- A synthetic vcf with all possible SNVs in the mtDNA reference sequence NC_012920.1, annotated with VEP, in the directory `synthetic_vcf`. Note that the files `NC_012920.1_synthetic_vep.vcf.gz` and `NC_012920.1_synthetic_vep_splitvarstwogenes.vcf.gz` will need to be unzipped before using.

Note that the gnomAD data files used (VCF with genotype and heteroplasmy level data for each sample, and txt with sample annotations) are not publicly available. The path to the directories with the required files (`insilicos`,`other_databases`,`synthetic_vcf`, and the input gnomAD data) can be specified via `reformat_vcf.py`, the default follows the directory structure of this repo. 

The script produces two output files: (1) `reformated.vcf`, which is used as the input file for figures 5, 6, S5D, S7, S8 and table S3, the input for `generate_tableS4.R` is generated via the second output file, (2) `annotated_synthetic.vcf`. The latter file is provided in this repo.

### Generate plots for figures 5, 6, S5D, S7, S8 and tables 3, 4
The figures were generated in R. The script `generate_all_figs5_6_S5d_S7_S8_tables_S3_S4.R` can be used to produce all specified figures and tables. Alternatively, the figures can be generated individually using the relevant scripts (named for the figure or table they produce). Note that `options(bitmapType = 'cairo', device = 'png')` can be commented/uncommented depending on if you are using Rstudio or the terminal to run the R scripts (relevant for collate_figure5.R, generate_figS7.R, generate_figS8.R, generate_fig6.R, generate_figS5D.R, generate_all_figs5_6_S5d_S7_S8_tables_S3_S4.R scripts). 

**Version details**
R v3.6.1: `cowplot 1.0.0` `dplyr 0.8.5` `forcats 0.5.0` `ggbeeswarm 0.6.0` `ggplot2 2_3.3.0 ` `ggpubr 0.3.0` `ggrepel 0.8.2` `gridExtra 2.3` `plyr 1.8.6` `tidyr 1.1.0` `scales 1.1.1` `stringr 1.4.0`




