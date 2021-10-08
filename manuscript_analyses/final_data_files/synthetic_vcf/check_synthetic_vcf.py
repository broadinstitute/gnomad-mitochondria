#!/usr/bin/python
# Run as
# python <check_synthetic_vcf.py>

import sys


def make_synthetic_vcf():
    """Generate a vcf with every possible single nucleotide variant in the mtDNA."""
    with open("NC_012920.1.fasta") as fasta:
        fasta = (
            fasta.read()
            .replace(">NC_012920.1 Homo sapiens mitochondrion, complete genome", "")
            .replace("\n", "")
        )  # get rid of header

        file = open("alternate_synthetic.vcf", "w")

        pos = 0

        for base in fasta:
            pos += 1
            ref = fasta[int(pos) - 1]
            if ref == "A":
                file.write(
                    "MT"
                    + "\t"
                    + str(pos)
                    + "\t"
                    + str(pos)
                    + "\t"
                    + (ref + "/" + "T")
                    + "\t"
                    + "+"
                    + "\n"
                )
                file.write(
                    "MT"
                    + "\t"
                    + str(pos)
                    + "\t"
                    + str(pos)
                    + "\t"
                    + (ref + "/" + "C")
                    + "\t"
                    + "+"
                    + "\n"
                )
                file.write(
                    "MT"
                    + "\t"
                    + str(pos)
                    + "\t"
                    + str(pos)
                    + "\t"
                    + (ref + "/" + "G")
                    + "\t"
                    + "+"
                    + "\n"
                )
            elif ref == "C":
                file.write(
                    "MT"
                    + "\t"
                    + str(pos)
                    + "\t"
                    + str(pos)
                    + "\t"
                    + (ref + "/" + "T")
                    + "\t"
                    + "+"
                    + "\n"
                )
                file.write(
                    "MT"
                    + "\t"
                    + str(pos)
                    + "\t"
                    + str(pos)
                    + "\t"
                    + (ref + "/" + "A")
                    + "\t"
                    + "+"
                    + "\n"
                )
                file.write(
                    "MT"
                    + "\t"
                    + str(pos)
                    + "\t"
                    + str(pos)
                    + "\t"
                    + (ref + "/" + "G")
                    + "\t"
                    + "+"
                    + "\n"
                )
            elif (
                ref == "G" or ref == "N"
            ):  # arbitrary grouping of N with G, to keep the spacer position m.3107 which is N in the reference sequence
                file.write(
                    "MT"
                    + "\t"
                    + str(pos)
                    + "\t"
                    + str(pos)
                    + "\t"
                    + (ref + "/" + "T")
                    + "\t"
                    + "+"
                    + "\n"
                )
                file.write(
                    "MT"
                    + "\t"
                    + str(pos)
                    + "\t"
                    + str(pos)
                    + "\t"
                    + (ref + "/" + "A")
                    + "\t"
                    + "+"
                    + "\n"
                )
                file.write(
                    "MT"
                    + "\t"
                    + str(pos)
                    + "\t"
                    + str(pos)
                    + "\t"
                    + (ref + "/" + "C")
                    + "\t"
                    + "+"
                    + "\n"
                )
            elif ref == "T":
                file.write(
                    "MT"
                    + "\t"
                    + str(pos)
                    + "\t"
                    + str(pos)
                    + "\t"
                    + (ref + "/" + "C")
                    + "\t"
                    + "+"
                    + "\n"
                )
                file.write(
                    "MT"
                    + "\t"
                    + str(pos)
                    + "\t"
                    + str(pos)
                    + "\t"
                    + (ref + "/" + "A")
                    + "\t"
                    + "+"
                    + "\n"
                )
                file.write(
                    "MT"
                    + "\t"
                    + str(pos)
                    + "\t"
                    + str(pos)
                    + "\t"
                    + (ref + "/" + "G")
                    + "\t"
                    + "+"
                    + "\n"
                )

        file.close()


if __name__ == "__main__":
    print("generating a synthetic vcf for the mtDNA!")
    make_synthetic_vcf()
