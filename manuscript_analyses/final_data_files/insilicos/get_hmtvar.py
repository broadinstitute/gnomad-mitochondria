#!/usr/bin/python
# using this to get Hmtvar annotations

import csv
import decimal
import json
import sys
import urllib


def hmtvar():
    """Retrieve annotations for every possible mtDNA SNV from the HmtVar database via API."""
    with open("../synthetic_vcf/NC_012920.1_synthetic_vep.vcf") as vcf:
        vcf = csv.reader(vcf, delimiter="\t")

        f = open("hmtvar_annotations.txt", "w")
        header = "POS	REF	ALT	HmtVar"
        f.write(header + "\n")

        for row in vcf:
            if not row[0].startswith("#"):
                pos = row[1]
                ref = row[3]
                alt = row[4]
                variant = ref + pos + alt

                url = "https://www.hmtvar.uniba.it/api/main/mutation/" + variant

                output = urllib.urlopen(url)
                annotation = output.read()

                f.write(
                    str(pos)
                    + "\t"
                    + str(ref)
                    + "\t"
                    + str(alt)
                    + "\t"
                    + str(annotation)
                )

        f.close()


if __name__ == "__main__":
    print("getting HmtVar annotations!")
    hmtvar()
