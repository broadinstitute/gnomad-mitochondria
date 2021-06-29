#Make a synthetic vcf

#Downloaded fasta file of mtDNA reference sequence from https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1?report=fasta

#Fasta file has a header line, and empty last line, so want to remove these
#1d deletes the first line (1 to only act on the first line, d to delete it)
#$d deletes the last line ($ to only act on the last line, d to delete it)
#Using sed replace every character with itself followed by a newline (doesn't include empty lines for linefeeds)
#This still maintains some empty lines, presumably from the line breaks. To remove any empty lines, use grep
#Then, use nl to number each line

sed -e '1d;$d' NC_012920.1.fasta | sed -e 's/\(.\)/\1\'$'\n/g' | grep -v "^$" | nl > temp_NC_012920.1.txt

# wc -l should be 16569, number of bp in the mtDNA human reference genome

awk 'OFS="\=" {
	if ($2 == "A")
	{print $1 "\t" $2 "\t" "T" "\n" $1 "\t" $2 "\t" "C" "\n" $1 "\t" $2 "\t" "G"}
	else if ($2 == "C")
	{print $1 "\t" $2 "\t" "T" "\n" $1 "\t" $2 "\t" "A" "\n" $1 "\t" $2 "\t" "G"}
	else if ($2 == "T")
	{print $1 "\t" $2 "\t" "C" "\n" $1 "\t" $2 "\t" "A" "\n" $1 "\t" $2 "\t" "G"}
	else # so if G (or N)
	{print $1 "\t" $2 "\t" "T" "\n" $1 "\t" $2 "\t" "A" "\n" $1 "\t" $2 "\t" "C"}
}' temp_NC_012920.1.txt > temp2_NC_012920.1.txt

# now annotate for vep input

cat temp2_NC_012920.1.txt | awk 'OFS="\=" {print $1 "\t" $1 "\t" $2 "/" $3}' | sed -e 's/^/MT	/' -e 's/$/	+/' > NC_012920.1_synthetic.vcf

echo "Ready for VEP!"
rm temp*
