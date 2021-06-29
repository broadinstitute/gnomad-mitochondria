#!/usr/bin/python

import csv 
import sys 
import itertools
import decimal
import argparse
csv.field_size_limit(sys.maxsize)

'''gather additional annotations for gnomAD mitochondrial variants, not in vcf:
MITOMAP: disease-associated variants, and database allele frequencies
APOGEE: in silico prediction for non-synonymous variants
HmtVar: in silico prediction for tRNA variants'''

def in_mitomap(other_databases_path): #disease-associated variant list, and MITOMAP database allele frequencies
	mitomap_counts = {}
	matches_disease = {}

	with open(other_databases_path + 'MITOMAP_polymorphisms_02022021.cgi') as csv_file:
		var_list = csv.DictReader(csv_file,delimiter='\t')

		for row in var_list: #for each variant in mitomap
			#need to use allele counts to convert to allele frequency, as if count = 1, gbfreq is written as 0.0 (rounded) in MITOMAP, total used for denominator n=51836
			mitomap_counts[str(row["pos"]),str(row["ref"]),str(row["alt"])] = float(decimal.Decimal(row["gbcnt"])/51836)
	
	with open(other_databases_path + 'MITOMAP_disease_02012021.cgi') as csv_file:
		var_list = csv.DictReader(csv_file,delimiter='\t')

		for row in var_list: 
			if (str(row["pos"]),str(row["ref"]),str(row["alt"])) not in mitomap_counts:
				mitomap_counts[str(row["pos"]),str(row["ref"]),str(row["alt"])] = float(decimal.Decimal(row["gbcnt"])/51836)
			if (row["status"] == "Cfrm") or (row["status"] == "Reported"): #gather disease-associated variants
				matches_disease[str(row["pos"]),str(row["ref"]),str(row["alt"])] = (row["status"],row["heteroplasmy"],row["homoplasmy"],row["disease"])

	return(mitomap_counts,matches_disease)

def apogee(insilicos_path): #in silico prediction for non-synonymous variants
	with open(insilicos_path + 'MitImpact_db_3.0.6.txt') as mitimpact:
		mitimpact = csv.DictReader(mitimpact, delimiter='\t')

		matches_apogee = {}

		for row in mitimpact:
			matches_apogee[(row["Start"],row["Ref"],row["Alt"],row["Gene_symbol"])] = (row["APOGEE"])

	return(matches_apogee)

def hmtvar(insilicos_path): #in silico prediction for tRNA variants, retrieved via api using get_hmtvar.py
	with open(insilicos_path + 'hmtvar_annotations.txt') as hmtvar:
		hmtvar = csv.DictReader(hmtvar, delimiter='\t')

		matches_hmtvar = {}

		for row in hmtvar:
			annotation = row["HmtVar"]
			insilico = ""
			#extract the in silico prediction from the annotation
			if len(annotation)>3:
				for character in annotation.split("pathogenicity")[1].split(",")[0]:
					if character.isalnum():
						insilico += character

			matches_hmtvar[(row["POS"],row["REF"],row["ALT"])] = (insilico)

	return(matches_hmtvar)

'''extract relevant data from gnomAD sample sheet and vcf:
haplogroup of each sample, heteroplasmy of each variant call, haplogroup associated with each variant call, allele frequencies and maximum heteroplasmy, and vep annotations'''

def samples(gnomAD_path): #create dictionary of sample ID and haplogroup
	with open(gnomAD_path + 't21/sample_annotations_gnomad.txt') as csv_file:
		samples = csv.DictReader(csv_file,delimiter='\t')

		matches_samples = {}

		for row in samples:
			matches_samples[row["s"]] = (row["hap"])

	return(matches_samples)

def parse_vcf(matches_samples,gnomAD_path): #vcf with genotype and heteroplasmy level data for each sample, not publicly available
	with open(gnomAD_path + 't21/sample_annotations_gnomad.vcf') as tsv_file:
		vcf = csv.reader(tsv_file,delimiter='\t')
		#vep annotations are per header
		header = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|ALLELE_NUM|DISTANCE|STRAND|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info"
		matches_annotations = {} #extracting vep, allele frequencies, maximum heteroplasmy and other annotations from vcf
		matches_base = {} #this is the maximum heteroplasmy of a SNV at each base
		matches_codon = {} #this is the maximum heteroplasmy of a non-synonymous SNV at each codon in protein-coding genes
		matches_insilico = {} #extracting in silico annotations in vcf
		sample_ids = {} #this is used to annotate the haplogroup of the individuals the variant is called in

		for row_index, row in enumerate(itertools.islice(vcf, 68, None)): #row with header and sample IDs is 69 for t21
			all_het = [] #resets every row/variant
			haplos = []
			#this step extracts the heteroplasmy level of variant calls, and the haplogroup of the samples with the variant call
			if (row[6] == "PASS") or (row_index == 0): #pass only sites and header
				for col_index, cell in enumerate(itertools.islice(row, 9, None)): #start from column 10
					if row_index == 0: #ie if the first row (the header)
						sample_ids[col_index] = cell #building a dictionary of the sample IDs
					else:
						genotype = cell.split(':')[0] 
						if (genotype == "0/1") or (genotype == "1/1"): 
							col_header = sample_ids[col_index] #sample
							haplogroup = matches_samples[col_header] #haplogroup of sample
							haplos.append(haplogroup)
							heteroplasmy = cell.split(':')[2] 
							all_het.append(heteroplasmy)

			if (row[6] == "PASS"):
				POS = row[1]
				REF = row[3]
				ALT = row[4]
				info = row[7] #column 8
				max_hl = info.split(';max_hl=')[1].split(';')[0]
				AF_hom = info.split(';AF_hom=')[1].split(';')[0] 
				AF_het = info.split(';AF_het=')[1].split(';')[0] 
				faf_hapmax_hom = info.split(';faf_hapmax_hom=')[1].split(';')[0] 
				AN = info.split(';AN=')[1].split(';')[0] 
				AC_hom = info.split(';AC_hom=')[1].split(';')[0] 
				AC_het = info.split(';AC_het=')[1].split(';')[0] 
				vep = info.split('vep=')[1].split(';')[0]
				VARIANT_CLASS = vep.split('|')[20]
				SYMBOL = vep.split('|')[3]
				BIOTYPE = vep.split('|')[7]
				Consequence = vep.split('|')[1]
				Protein_position = vep.split('|')[14] #residue/codon
				HGVSc = vep.split('|')[10] 
				HGVSp = vep.split('|')[11] 

				max_hap_AF_hom = max_pop_AF_hom = 0.0 #reset each variant
				#maximum AF_hom per haplogroup, and per population
				for value in info.split(';hap_AF_hom=')[1].split(';')[0].split('|'):
					if (value!='.') and (float(value) >= float(max_hap_AF_hom)):
						max_hap_AF_hom = value
				for value in info.split(';pop_AF_hom=')[1].split(';')[0].split('|'):
					if (value!='.') and (float(value) >= float(max_pop_AF_hom)):
						max_pop_AF_hom = value

				matches_annotations[POS,REF,ALT,SYMBOL] = (max_hl,AN,AC_hom,AC_het,AF_hom,AF_het,faf_hapmax_hom,VARIANT_CLASS,SYMBOL,BIOTYPE,Consequence,Protein_position,HGVSc,HGVSp,all_het,haplos,max_hap_AF_hom,max_pop_AF_hom)

				#for max heteroplasmy at base #for SNVs only
				if VARIANT_CLASS == "SNV":
					if (not POS in matches_base) or (matches_base[POS] < max_hl): 
						matches_base[POS] = (max_hl)
				#for max heteroplasmy at codon #for SNVs only
				if (BIOTYPE == "protein_coding") and (VARIANT_CLASS == "SNV") and (Consequence != "synonymous_variant"):
					if (not (SYMBOL,Protein_position) in matches_codon) or (matches_codon[(SYMBOL,Protein_position)] < max_hl): 
						matches_codon[(SYMBOL,Protein_position)] = (max_hl)
				#tRNA in silicos in vcf
				if (BIOTYPE == "Mt_tRNA") and (VARIANT_CLASS == "SNV"): 
					if not (SYMBOL,POS,ALT) in matches_insilico:
						mitotip = info.split(';mitotip_trna_prediction=')[1].split(';')[0]
						pon_mt = info.split(';pon_mt_trna_prediction=')[1].split(';')[0]
						matches_insilico[(SYMBOL,POS,ALT)] = (mitotip,pon_mt)

				if vep.count('|') > header.count('|'): #then is two annotations, variant lies in two genes
					SYMBOL = vep.split('|')[48] #plus 45, header.count('|') = 44
					BIOTYPE = vep.split('|')[52]
					Consequence = vep.split('|')[46]
					Protein_position = vep.split('|')[59] 
					HGVSc = vep.split('|')[55] 
					HGVSp = vep.split('|')[56] 

					#create a second entry for variants in two genes, with two consequences
					matches_annotations[POS,REF,ALT,SYMBOL] = (max_hl,AN,AC_hom,AC_het,AF_hom,AF_het,faf_hapmax_hom,VARIANT_CLASS,SYMBOL,BIOTYPE,Consequence,Protein_position,HGVSc,HGVSp,all_het,haplos,max_hap_AF_hom,max_pop_AF_hom)

					#for max heteroplasmy at codon #for SNVs only
					#to catch variants where one consequence but not other is protein-changing
					if (BIOTYPE == "protein_coding") and (VARIANT_CLASS == "SNV") and (Consequence != "synonymous_variant"): 
						if (not (SYMBOL,Protein_position) in matches_codon) or (matches_codon[(SYMBOL,Protein_position)] < max_hl): 
							matches_codon[(SYMBOL,Protein_position)] = (max_hl)

					#tRNA in silicos in vcf
					#to catch variants in two different tRNA genes
					if (BIOTYPE == "Mt_tRNA") and (VARIANT_CLASS == "SNV"): 
						if not (SYMBOL,POS,ALT) in matches_insilico:
							mitotip = info.split(';mitotip_trna_prediction=')[1].split(';')[0]
							pon_mt = info.split(';pon_mt_trna_prediction=')[1].split(';')[0]
							matches_insilico[(SYMBOL,POS,ALT)] = (mitotip,pon_mt)

	return(matches_annotations,matches_base,matches_codon,matches_insilico)

'''now generate 'reformated.vcf', which is used to produce fig5, fig6, figS4d, figS6, figS7 and table S2
this includes the above gathered annotations from gnomAD vcf and other sources'''

def write_file_for_figures(matches_annotations,matches_disease,matches_base,matches_codon,matches_insilico,matches_apogee,matches_hmtvar,gnomAD_path): 
	with open(gnomAD_path + 't21/sample_annotations_gnomad.vcf') as tsv_file:
		vcf = csv.reader(tsv_file,delimiter='\t')
		
		file=open('reformated.vcf', "w")
		header = "POS	REF	ALT	max_hl	AN	AC_hom	AC_het	AF_hom	AF_het	faf_hapmax_hom	VARIANT_CLASS	SYMBOL	BIOTYPE	Consequence	Protein_position	HGVSc	HGVSp	all_hl	all_haplogroups	max_hap_AF_hom	max_pop_AF_hom	max_hl_SNV_base	max_hl_nonsyn_SNV_codon	Mitomap_dz_status	Mitomap_dz_heteroplasmy	Mitomap_dz_homoplasmy	Mitomap_disease	APOGEE	Mitotip	Pon_mt_trna	Hmtvar"
		file.write(header + '\n')

		for row in vcf:
			if not row[0].startswith('#'):
				if row[6] == "PASS": #pass only sites
					POS = row[1]
					REF = row[3]
					ALT = row[4]
					vep = row[7].split('vep=')[1].split(';')[0]
					SYMBOL = vep.split('|')[3]
					Protein_position = vep.split('|')[14] 
				
					if POS in matches_base: #this is the maximum heteroplasmy of a SNV at each base
						max_base = matches_base[POS]
					else:
						max_base = 0

					if (SYMBOL,Protein_position) in matches_codon: #this is the maximum heteroplasmy of a non-synonymous SNV at each codon in protein-coding genes
						max_codon = matches_codon[(SYMBOL,Protein_position)]
					else:
						max_codon = 0

					if (POS,REF,ALT) in matches_disease: #this is MITOMAP disease-associated variant
						mitomap_status = matches_disease[(POS,REF,ALT)][0]
						mitomap_het = matches_disease[(POS,REF,ALT)][1]
						mitomap_hom = matches_disease[(POS,REF,ALT)][2]
						mitomap_dz = matches_disease[(POS,REF,ALT)][3]
					else:
						mitomap_status = mitomap_het = mitomap_hom = mitomap_dz = ""

					if (SYMBOL,POS,ALT) in matches_insilico: #this is tRNA in silicos from the gnomAD vcf
						mitotip = matches_insilico[(SYMBOL,POS,ALT)][0]
						pon_mt = matches_insilico[(SYMBOL,POS,ALT)][1]
					else:
						mitotip = pon_mt = ""

					if (POS,REF,ALT,SYMBOL) in matches_apogee: #this is additional in silico for non-synonymous variants
						apogee = matches_apogee[(POS,REF,ALT,SYMBOL)]
					else:
						apogee = ""

					if ((POS,REF,ALT) in matches_hmtvar) and (SYMBOL.startswith('MT-T')): #this is additional in silico for tRNA variants
						hmtvar = matches_hmtvar[(POS,REF,ALT)]
					else:
						hmtvar = ""

					file.write(str(POS) + '\t' + str(REF) + '\t' + str(ALT) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][0]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][1]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][2]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][3]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][4]) + 
							'\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][5]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][6]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][7]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][8]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][9]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][10]) +
							'\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][11]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][12]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][13]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][14]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][15]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][16]) + 
							'\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][17]) + '\t' + str(max_base) + '\t' + str(max_codon) + '\t' + str(mitomap_status) + '\t' + str(mitomap_het) + '\t' + str(mitomap_hom) + '\t' + str(mitomap_dz) + 
							'\t' + str(apogee) + '\t' + str(mitotip) + '\t' + str(pon_mt) + '\t' + str(hmtvar) + '\n')

					if vep.count('|') > 44: #if two annotations, variant lies in two genes
						SYMBOL = vep.split('|')[48]
						Protein_position = vep.split('|')[59] 

						if (SYMBOL,Protein_position) in matches_codon: #this is the maximum heteroplasmy of a non-synonymous SNV at each codon in protein-coding genes
							max_codon = matches_codon[(SYMBOL,Protein_position)]
						else:
							max_codon = 0

						if (SYMBOL,POS,ALT) in matches_insilico: #this is tRNA in silicos from the gnomAD vcf
							mitotip = matches_insilico[(SYMBOL,POS,ALT)][0]
							pon_mt = matches_insilico[(SYMBOL,POS,ALT)][1]
						else:
							mitotip = pon_mt = ""

						if (POS,REF,ALT,SYMBOL) in matches_apogee: #this is additional in silico for non-synonymous variants
							apogee = matches_apogee[(POS,REF,ALT,SYMBOL)]
						else:
							apogee = ""

						if ((POS,REF,ALT) in matches_hmtvar) and (SYMBOL.startswith('MT-T')): #this is additional in silico for tRNA variants
							hmtvar = matches_hmtvar[(POS,REF,ALT)]
						else:
							hmtvar = ""

						file.write(str(POS) + '\t' + str(REF) + '\t' + str(ALT) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][0]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][1]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][2]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][3]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][4]) + 
							'\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][5]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][6]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][7]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][8]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][9]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][10]) +
							'\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][11]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][12]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][13]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][14]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][15]) + '\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][16]) + 
							'\t' + str(matches_annotations[POS,REF,ALT,SYMBOL][17]) + '\t' + str(max_base) + '\t' + str(max_codon) + '\t' + str(mitomap_status) + '\t' + str(mitomap_het) + '\t' + str(mitomap_hom) + '\t' + str(mitomap_dz) + 
							'\t' + str(apogee) + '\t' + str(mitotip) + '\t' + str(pon_mt) + '\t' + str(hmtvar) + '\n')

'''now generate a second output file 'annotated_synthetic.vcf', which is used to produce table S3
first use a VEP annotated synthetic vcf for all possible SNVs in the mtDNA to create dictionaries of the consequences and locus of each possible SNV
gather allele frequency and maximum heteroplasmy data from HelixMTdb for annotating, also annoate with MITOMAP allele frequency and disease-assocated variants gathered above'''

def consequences(synthetic_vcf_path): #create dictionaries of the consequences and locus of each possible SNV in mtDNA
	with open(synthetic_vcf_path + 'NC_012920.1_synthetic_vep_splitvarstwogenes.vcf') as vcf: #this file has been edited so that for variants in two genes, the consequence on each gene is listed on a separate line using split_vars_two_genes.py
		vcf = csv.DictReader(vcf,delimiter='\t')

		matches_locus = {} #gene name
		matches_conseq = {} #ie synonymous, missense, etc
		matches_HGVSp = {} #HGVSp
		matches_HGVSc = {} #HGVSc

		for row in vcf:
			base = str(row["POS"])
			alt = row["ALT"]
			gene = row["SYMBOL"]

			if not base in matches_locus:
				matches_locus[base] = []
				if gene: #ie if base lies in a gene
					matches_locus[base] = [gene] 
				else:
					matches_locus[base] = ["non-coding"]		
			else: #if base lies in two genes
				if not gene:
					gene = "non-coding"
				if not gene in matches_locus[base]: #so only adding different gene names
					matches_locus[base].append(gene)

			if not (base,alt) in matches_conseq: #every base position will have a consequence
				matches_conseq[(base,alt)] = [row["Consequence"]]
			else:
				matches_conseq[(base,alt)].append(row["Consequence"])

			if not (base,alt) in matches_HGVSp: 
				matches_HGVSp[(base,alt)] = [row["HGVSp"]]
			else:
				matches_HGVSp[(base,alt)].append(row["HGVSp"])

			if not (base,alt) in matches_HGVSc: 
				matches_HGVSc[(base,alt)] = [row["HGVSc"]]
			else:
				matches_HGVSc[(base,alt)].append(row["HGVSc"])

	return (matches_locus,matches_conseq,matches_HGVSc,matches_HGVSp)

def in_helix(other_databases_path): #allele frequency and maximum heteroplasmy data from HelixMTdb database
	with open(other_databases_path + 'HelixMTdb_20200327.tsv') as csv_file:
		var_list = csv.DictReader(csv_file,delimiter='\t')

		helix_counts = {}

		for row in var_list:
			pos = str(row["locus"].split("chrM:")[1])
			ref = str(row["alleles"].split("\"")[1])
			alt = str(row["alleles"].split("\"")[3])

			if float(row["AF_hom"]) > 0:
				max_het = 1
			elif float(row["AF_hom"])==0:
				max_het = float(row["max_ARF"])

			helix_counts[pos,ref,alt] = (max_het,row["AF_hom"],row["AF_het"])
	return (helix_counts)

def annotate_syn_vcf(matches_annotations,helix_counts,mitomap_counts,matches_locus,matches_conseq,matches_HGVSc,matches_HGVSp,synthetic_vcf_path):
	with open(synthetic_vcf_path + 'NC_012920.1_synthetic_vep.vcf') as csv_file:
		synt_vcf = csv.reader(csv_file,delimiter='\t')

		file=open('annotated_synthetic.vcf', "w")
		header = "POS	REF	ALT	SYMBOL	Consequence	HGVSc	HGVSp	gnomAD_max_hl	gnomAD_AF_hom	gnomAD_AF_het	Helix_max_hl	Helix_af_hom	Helix_af_het	Mitomap_af	Mitomap_dz_status	Mitomap_dz_heteroplasmy	Mitomap_dz_homoplasmy"
		file.write(header + '\n')

		for row in synt_vcf:
			if not row[0].startswith('#'):
				POS = row[1]
				REF = row[3]
				ALT = row[4]
				SYMBOL = row[7].split('|')[3]
				if (POS,REF,ALT,SYMBOL) in matches_annotations: #annotations from gnomAD vcf
					in_gnomad_max = matches_annotations[(POS,REF,ALT,SYMBOL)][0]
					in_gnomad_afhom = matches_annotations[(POS,REF,ALT,SYMBOL)][4]
					in_gnomad_afhet = matches_annotations[(POS,REF,ALT,SYMBOL)][5]
				else:
					in_gnomad_max = in_gnomad_afhom = in_gnomad_afhet = 0

				if (POS,REF,ALT) in helix_counts: #annotations from HelixMTdb
					in_helix_max = helix_counts[(POS,REF,ALT)][0]
					in_helix_afhom = helix_counts[(POS,REF,ALT)][1]
					in_helix_afhet = helix_counts[(POS,REF,ALT)][2]
				else:
					in_helix_max = in_helix_afhom = in_helix_afhet = 0

				if (POS,REF,ALT) in mitomap_counts: #annotations from MITOMAP
					in_mitomap = mitomap_counts[(POS,REF,ALT)]
				else:
					in_mitomap = 0

				if (POS,REF,ALT) in matches_disease: #disease-associated variants from MITOMAP
					mitomap_status = matches_disease[(POS,REF,ALT)][0]
					mitomap_het = matches_disease[(POS,REF,ALT)][1]
					mitomap_hom = matches_disease[(POS,REF,ALT)][2]
				else:
					mitomap_status = mitomap_het = mitomap_hom = ""

				file.write(str(POS) + '\t' + str(REF) + '\t' + str(ALT) + '\t' + str(matches_locus[POS]).strip('[]').replace("'","") + '\t' + str(matches_conseq[(POS,ALT)]).strip('[]').replace("'","") + '\t' + str(matches_HGVSc[(POS,ALT)]).strip('[]').replace("'","") + '\t' + str(matches_HGVSp[(POS,ALT)]).strip('[]').replace("'","") + '\t' + 
					str(in_gnomad_max) + '\t' + str(in_gnomad_afhom) + '\t' + str(in_gnomad_afhet) + '\t' + str(in_helix_max) + '\t' + str(in_helix_afhom) + '\t' + str(in_helix_afhet) + '\t' + str(in_mitomap) + '\t' + str(mitomap_status) + '\t' + str(mitomap_het) + '\t' + str(mitomap_hom) + '\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-other_databases_path', help='path to directory with MITOMAP and HelixMTdb files', required=False, default='final_data_files/other_databases/')
    parser.add_argument('-insilicos_path', help='path to directory with APOGEE and HmtVar files', required=False, default='final_data_files/insilicos/')
    parser.add_argument('-synthetic_vcf_path', help='path to directory with synthetic vcf files', required=False, default='final_data_files/synthetic_vcf/')
    parser.add_argument('-gnomAD_path', help='path to directory with gnomAD files (not publicly available', required=False, default='final_data_files/gnomAD/')
    args = parser.parse_args()

    print("Starting!\nThis script will produce 2 output files: reformated.vcf and annotated_synthetic.vcf, which are the input files for the R scripts used to produce the figures and tables")

    #gather annotations for output "reformated.vcf", these are additional insilico and population frequency annotations
    (mitomap_counts,matches_disease) = in_mitomap(args.other_databases_path)
    matches_apogee = apogee(args.insilicos_path)
    matches_hmtvar = hmtvar(args.insilicos_path)
    print("gathered annotations!\nNow will parse the gnomAD data, this might take ~15 minutes")
    #these are annotations extracted from the gnomAD data
    matches_samples = samples(args.gnomAD_path)
    (matches_annotations,matches_base,matches_codon,matches_insilico) = parse_vcf(matches_samples,args.gnomAD_path)
    #now write the output file "reformated.vcf", which is used to produce fig5, fig6, figS4d, figS6, figS7 and table S2
    write_file_for_figures(matches_annotations,matches_disease,matches_base,matches_codon,matches_insilico,matches_apogee,matches_hmtvar,args.gnomAD_path)
    print("written file reformated.vcf!\nNow will make the second output file annotated_synthetic.vcf")

    #now generate a second output file "annotated_synthetic.vcf", which is used to produce table S3, this output needs additional annotations vs "reformated.vcf", including from HelixMTdb
    (matches_locus,matches_conseq,matches_HGVSc,matches_HGVSp) = consequences(args.synthetic_vcf_path)
    helix_counts = in_helix(args.other_databases_path)
    annotate_syn_vcf(matches_annotations,helix_counts,mitomap_counts,matches_locus,matches_conseq,matches_HGVSc,matches_HGVSp,args.synthetic_vcf_path)
    print("written file annotated_synthetic.vcf!\nScript is now complete")
