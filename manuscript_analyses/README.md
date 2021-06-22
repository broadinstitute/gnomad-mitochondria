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

 