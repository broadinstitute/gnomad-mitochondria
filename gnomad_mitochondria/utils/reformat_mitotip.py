#!/usr/bin/env python

import argparse
from subprocess import check_output


def main(args):
    mitotip_scores = args.mitotip_scores
    mt_reference = args.mt_reference
    output_file = args.output_file

    out_file = open(output_file, "w")

    with open(mitotip_scores, "r") as f:
        header = next(f)
        out_file.write(header)
        for line in f:
            items = line.rstrip().split("\t")
            position, ref, alt, score, quartile, count, percent, status = items[0:8]

            # ":" indicates a deletion for the given base in the ref column at that position, need reformatting
            if alt == ":":
                previous_position = int(position) - 1
                ref_info = check_output(
                    [
                        "samtools",
                        "faidx",
                        mt_reference,
                        f"MT:{previous_position}-{previous_position}",
                    ]
                )

                # Second line will contain the reference base from the query
                previous_base = ref_info.decode().split("\n")[1]

                new_ref = f"{previous_base}{ref}"
                new_alt = previous_base

                new_line = "\t".join(
                    [
                        str(previous_position),
                        new_ref,
                        new_alt,
                        score,
                        quartile,
                        count,
                        percent,
                        status,
                    ]
                )
                out_file.write(new_line + "\n")
            else:
                out_file.write(line)

    out_file.close()


if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="This script reformats mitotip data so that single base deletions are compatible with mutect pipeline"
    )
    p.add_argument(
        "-m", "--mitotip-scores", help="Downloaded txt file from mitotip", required=True
    )
    p.add_argument(
        "-r",
        "--mt-reference",
        help="Path to mitochondrial reference fasta",
        required=True,
    )
    p.add_argument(
        "-o", "--output-file", help="Name of output file for results", required=True
    )

    args = p.parse_args()

    main(args)
