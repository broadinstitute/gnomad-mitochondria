#!/usr/bin/python
# using this get Hmtvar annotations

import sys
import csv
import decimal
import urllib
import json


def hmtvar():
    with open(
        "required_files/synthetic_vcf/NC_012920.1_synthetic_vep_noheader.vcf"
    ) as vcf:
        vcf = csv.DictReader(vcf, delimiter="\t")

        f = open("Hmtvar/hmtvar_annotations.txt", "w")
        header = "POS	REF	ALT	HmtVar"
        f.write(header + "\n")

        for row in vcf:
            pos = row["POS"]
            ref = row["REF"]
            alt = row["ALT"]
            variant = ref + pos + alt

            url = "https://www.hmtvar.uniba.it/api/main/mutation/" + variant

            output = urllib.urlopen(url)
            annotation = output.read()
            print(annotation.decode("utf-8"))

            f.write(
                str(pos)
                + "\t"
                + str(ref)
                + "\t"
                + str(alt)
                + "\t"
                + str(annotation)
                + "\n"
            )


hmtvar()
