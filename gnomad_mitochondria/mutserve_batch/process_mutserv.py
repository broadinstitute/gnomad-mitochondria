#!/usr/bin/env python
import argparse
import logging
import os
import pysam
import re
import statistics
import sys

from io import TextIOWrapper
from os.path import basename, splitext
from subprocess import Popen, PIPE, check_output


logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("run mutserv")
logger.setLevel(logging.INFO)


def evaluate_variant(
    alt: str,
    pos: int,
    variant_level: float,
    dp: int,
    ref: str,
    filters: str,
    results_out: TextIOWrapper,
    line_count: int,
    num_lines: int,
    input_vcf: TextIOWrapper,
    mt_reference: str,
) -> None:
    """
    Check if variant is a SNP or insertion (output as is) or deletion (initiate deletion function).

    Note: Initiating the deletion function will collapse multiple deletion calls into one deletion when certain conditions are met. This functionality is needed as deletions are output at individuals bases, but if it seems likely that several deletions are actually one large deletion, we want to represent the multiple deletions as just one deletion.

    :param alt: Alternate allele
    :param pos: Position of the variant
    :param variant_level: Heteroplasmy level of the variant
    :param dp: Depth of coverage at the given position
    :param ref: Reference allele
    :param filters: Filter fields
    :param results_out: Open filehandle where results should be written
    :param line_count: Line number for current iteration (count for row in the VCF)
    :param num_lines: Total number of lines in file
    :param input_vcf: Open filehandle for input vcf
    :param mt_reference: Path to mitochondria reference fasta
    :return: None
    """
    # If the alternate allele is "*", start the process to collapse the deletion, otherwise output the variant
    if alt == "*":
        initiate_deletion(
            pos,
            variant_level,
            dp,
            ref,
            results_out,
            line_count,
            num_lines,
            input_vcf,
            mt_reference,
        )
    else:
        results_out.write(
            f"MT\t{pos}\t.\t{ref}\t{alt}\t.\t{filters}\tAF=.\tGT:AD:VL:DP\t.:.:{variant_level}:{dp}\n"
        )


def format_variant_info(line_items: list) -> list:
    """
    Reformat VCF to pull out needed info to pass on to downstream functions and to correct types.

    :param line_items: Variant row of the VCF, with fields split into a list
    :return: List of position, reference allele, alternate allele, depth of coverage, variant level, and filters
    """
    pos = int(line_items[1])
    ref = line_items[3]
    alt = line_items[4]
    filters = line_items[6]
    sample_info = line_items[9]
    gt_fields = sample_info.split(":")
    gt = gt_fields[0]
    if gt == "1":
        dp = gt_fields[1]
        vl = 1.0
    else:
        dp = gt_fields[2]
        vl = gt_fields[1]

    dp = int(dp)
    vl = float(vl)

    return (pos, ref, alt, dp, vl, filters)


def close_deletion(
    first_deletion_position: int,
    num_deleted_bases: int,
    deletion_coverage_depths: list,
    deletion_variant_levels: list,
    mt_reference: str,
) -> list:
    """
    Close the deletion.

    Take input of single or multiple deletions and output one deletion variant to the VCF. If mulitple deletions are collapsed into one deletion, the variant level and DP are the average of all positions that constitute the deletion.

    :param first_deletion_position: First position of the deletion
    :param num_deleted_bases: Number of bases that are deleted
    :param deletion_coverage_depths: List of depth of coverage for all the positions of the deletion
    :param deletion_variant_levels: List of heteroplasmy levels for all the positions of the deletion
    :param mt_reference: Path to mitochondria reference fasta
    :return: first deletion position, reference bases, deletion_alt (first base in reference bases), average depth across the deletion positions, average variant level across the deletion positions, last position of the deletion
    """
    end_deletion = first_deletion_position + num_deleted_bases
    fasta_open = pysam.Fastafile(mt_reference)
    ref_bases = fasta_open.fetch("rCRS", first_deletion_position - 1, end_deletion)
    deletion_alt = ref_bases[
        0
    ]  # Grab just the first position of the ref bases to report in the ALT column of the VCF

    # Calculate the average depth across all positions of the deletion
    deletion_dp = int(statistics.mean(deletion_coverage_depths))
    # Calculate the average variant level across all positions of the deletion
    deletion_variant_level = statistics.mean(deletion_variant_levels)

    return (
        first_deletion_position,
        ref_bases,
        deletion_alt,
        deletion_dp,
        deletion_variant_level,
        end_deletion,
    )


def initiate_deletion(
    pos: int,
    variant_level: float,
    dp: int,
    ref: str,
    results_out: TextIOWrapper,
    line_count: int,
    num_lines: int,
    input_vcf: TextIOWrapper,
    mt_reference: str,
) -> None:
    """
    Initiate a deletion.

    Check next variant and if it's a deletion at the consecutive position with a varaint level +/- 10% of the current deletion,
    append that deletion to the current one (count as one larger deletion).
    For example, given three deletions at positions 2, 3, 4, such as (in format of POS REF ALT):
    2 A *
    3 A *
    4 C *
    If these deletions had similar variants levels, the three deletions would be reformatted to:
    1 TAAC T

    :param pos: Position
    :param variant_level: Heteroplasmy level
    :param dp: Depth of coverage at the given position
    :param ref: Reference allele
    :param results_out: Open filehandle where results should be written
    :param line_count: Line number for current iteration
    :param num_lines: Total number of lines in file
    :param input_vcf: Open filehandle for input vcf
    :param mt_reference: Path to mitochondria reference fasta
    :return: None
    """
    first_deletion_position = pos - 1
    num_deleted_bases = 1
    deletion_variant_levels = [variant_level]
    deletion_coverage_depths = [dp]
    deletion_pos = pos
    deletion_variant_level = variant_level

    # Close deletion if you've hit the last line of the VCF
    if line_count == num_lines:
        (
            first_deletion_position,
            ref_bases,
            deletion_alt,
            deletion_dp,
            deletion_variant_level,
            end_deletion,
        ) = close_deletion(
            first_deletion_position,
            num_deleted_bases,
            deletion_coverage_depths,
            deletion_variant_levels,
            mt_reference,
        )
        results_out.write(
            f"MT\t{first_deletion_position}\t.\t{ref_bases}\t{deletion_alt}\t.\tPASS\tAF=.\tGT:AD:VL:DP\t.:.:{deletion_variant_level}:{deletion_dp}\n"
        )

    else:
        for next_variant in input_vcf:
            line_count += 1
            next_variant = next_variant.rstrip()
            next_items = next_variant.split("\t")

            pos, ref, alt, dp, variant_level, filters = format_variant_info(next_items)

            # Extend the deletion if the bases are consecutive and the variant levels differ by no more than 10% heteroplasmy
            if (
                alt == "*"
                and pos == (deletion_pos + 1)
                and variant_level < (deletion_variant_level + 0.10)
                and variant_level > (deletion_variant_level - 0.10)
            ):
                deletion_pos = pos
                deletion_variant_level = variant_level
                num_deleted_bases += 1
                deletion_variant_levels.append(variant_level)
                deletion_coverage_depths.append(dp)

                # Close deletion if you've hit the last line of the VCF
                if line_count == num_lines:
                    (
                        first_deletion_position,
                        ref_bases,
                        deletion_alt,
                        deletion_dp,
                        deletion_variant_level,
                        end_deletion,
                    ) = close_deletion(
                        first_deletion_position,
                        num_deleted_bases,
                        deletion_coverage_depths,
                        deletion_variant_levels,
                        mt_reference,
                    )
                    results_out.write(
                        f"MT\t{first_deletion_position}\t.\t{ref_bases}\t{deletion_alt}\t.\tPASS\tAF=.\tGT:AD:VL:DP\t.:.:{deletion_variant_level}:{deletion_dp}\n"
                    )
            # When the bases of the deletion are not consecutive and/or the variant levels differ by more than 10% heteroplasmy, close the deletion
            else:
                (
                    first_deletion_position,
                    ref_bases,
                    deletion_alt,
                    deletion_dp,
                    deletion_variant_level,
                    end_deletion,
                ) = close_deletion(
                    first_deletion_position,
                    num_deleted_bases,
                    deletion_coverage_depths,
                    deletion_variant_levels,
                    mt_reference,
                )
                # Handle exception where deletion is at first reference position
                # From VCF spec: "must include the base before the event (which must be reflected in the POS field), unless the event occurs at position 1 on the contig in which case it must include the base after the event"
                # For example, a deletion of A at the first position would be represented as (in the order of CHROM POS REF ALT): chromosome_name 1  AG G
                if first_deletion_position == 0:
                    first_deletion_position = 1
                    one_over = end_deletion + 1

                    fasta_open = pysam.Fastafile(mt_reference)
                    ref_bases = fasta_open.fetch(
                        "rCRS", first_deletion_position - 1, one_over
                    )
                    deletion_alt = ref_bases[
                        0
                    ]  # Grab just the first position of the ref bases to report in the ALT column of the VCF

                results_out.write(
                    f"MT\t{first_deletion_position}\t.\t{ref_bases}\t{deletion_alt}\t.\tPASS\tAF=.\tGT:AD:VL:DP\t.:.:{deletion_variant_level}:{deletion_dp}\n"
                )

                evaluate_variant(
                    alt,
                    pos,
                    variant_level,
                    dp,
                    ref,
                    filters,
                    results_out,
                    line_count,
                    num_lines,
                    input_vcf,
                    mt_reference,
                )
                break


def main(args):
    input_file = args.input_file
    output_file = args.output_file
    mt_reference = args.mt_reference

    # Get the count of the number of lines in the file
    num_lines = check_output(["wc", "-l", f"{input_file}"])
    num_lines = num_lines.decode().split(" ")[0]
    num_lines = int(num_lines)

    # Open the file where output will be written no
    results_out = open(output_file, "w")

    logger.info("Reformatting VCF...")
    # Reformat VCF to be compatible with the combine_vcfs.py script
    line_count = 1
    with open(input_file, "r") as input_vcf:
        for line in input_vcf:
            line_count += 1
            if line.startswith("#"):
                results_out.write(line)
            else:
                line = line.rstrip()
                items = line.split("\t")
                # Format the VCF content
                pos, ref, alt, dp, vl, filters = format_variant_info(items)
                # Evaluate the variant (check if it's a SNP, insertion, or deletion) and collapse multiple deletions into one deletion if necessary
                evaluate_variant(
                    alt,
                    pos,
                    vl,
                    dp,
                    ref,
                    filters,
                    results_out,
                    line_count,
                    num_lines,
                    input_vcf,
                    mt_reference,
                )
    results_out.close()


if __name__ == "__main__":
    p = argparse.ArgumentParser("This script runs mutserv and reformats the output")
    p.add_argument(
        "-i",
        "--input-file",
        required=True,
        help="The VCF output by mtDNA-Server with multiallelic sites split",
    )
    p.add_argument(
        "-o", "--output-file", required=True, help="Name to use for the output file"
    )
    p.add_argument(
        "-r", "--mt-reference", required=True, help="Mitochondria rCRS reference fasta"
    )

    args = p.parse_args()

    main(args)
