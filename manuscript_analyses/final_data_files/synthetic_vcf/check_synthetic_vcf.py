#!/usr/bin/python
# Run as
#python <check_synthetic_vcf.py> 

import sys 

def make_synthetic_vcf():
	with open('NC_012920.1.fasta') as fasta:
		fasta = fasta.read().replace('>NC_012920.1 Homo sapiens mitochondrion, complete genome','').replace('\n','') #get rid of header

		file=open('alternate_synthetic.vcf', "w")

		POS = 0

		for base in fasta:
			POS += 1
			ref = fasta[int(POS)-1] #ie if base is not N, skip pos 3107
			if ref=="A":
				file.write("MT" + '\t' + str(POS) + '\t' + str(POS) + '\t' + (ref + '/' + "T") + '\t' + "+" + '\n')
				file.write("MT" + '\t' + str(POS) + '\t' + str(POS) + '\t' + (ref + '/' + "C") + '\t' + "+" + '\n')
				file.write("MT" + '\t' + str(POS) + '\t' + str(POS) + '\t' + (ref + '/' + "G") + '\t' + "+" + '\n')
			elif ref=="C":
				file.write("MT" + '\t' + str(POS) + '\t' + str(POS) + '\t' + (ref + '/' + "T") + '\t' + "+" + '\n')
				file.write("MT" + '\t' + str(POS) + '\t' + str(POS) + '\t' + (ref + '/' + "A") + '\t' + "+" + '\n')
				file.write("MT" + '\t' + str(POS) + '\t' + str(POS) + '\t' + (ref + '/' + "G") + '\t' + "+" + '\n')
			elif (ref=="G" or ref=="N"):
				file.write("MT" + '\t' + str(POS) + '\t' + str(POS) + '\t' + (ref + '/' + "T") + '\t' + "+" + '\n')
				file.write("MT" + '\t' + str(POS) + '\t' + str(POS) + '\t' + (ref + '/' + "A") + '\t' + "+" + '\n')
				file.write("MT" + '\t' + str(POS) + '\t' + str(POS) + '\t' + (ref + '/' + "C") + '\t' + "+" + '\n')
			elif ref=="T":
				file.write("MT" + '\t' + str(POS) + '\t' + str(POS) + '\t' + (ref + '/' + "C") + '\t' + "+" + '\n')
				file.write("MT" + '\t' + str(POS) + '\t' + str(POS) + '\t' + (ref + '/' + "A") + '\t' + "+" + '\n')
				file.write("MT" + '\t' + str(POS) + '\t' + str(POS) + '\t' + (ref + '/' + "G") + '\t' + "+" + '\n')

make_synthetic_vcf()