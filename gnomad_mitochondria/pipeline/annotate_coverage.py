import argparse
import logging
import math
import os
import re
import sys

import hail as hl

from os.path import dirname
from gnomad.utils.slack import slack_notifications
from hail.utils.java import info

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("Annotate coverage")
logger.setLevel(logging.INFO)


def multi_way_union_mts(mts: list, temp_dir: str, chunk_size: int) -> hl.MatrixTable:
    """
    Hierarchically join together MatrixTables in the provided list.

    :param mts: List of MatrixTables to join together
    :param temp_dir: Path to temporary directory for intermediate results
    :param chunk_size: Number of MatrixTables to join per chunk (the number of individual VCFs that should be combined at a time)
    :return: Joined MatrixTable
    """
    # Convert the MatrixTables to tables where entries are an array of structs
    staging = [mt.localize_entries("__entries", "__cols") for mt in mts]
    stage = 0
    while len(staging) > 1:
        # Calculate the number of jobs to run based on the chunk size
        n_jobs = int(math.ceil(len(staging) / chunk_size))
        info(f"multi_way_union_mts: stage {stage}: {n_jobs} total jobs")
        next_stage = []

        for i in range(n_jobs):
            # Grab just the tables for the given job
            to_merge = staging[chunk_size * i : chunk_size * (i + 1)]
            info(
                f"multi_way_union_mts: stage {stage} / job {i}: merging {len(to_merge)} inputs"
            )

            # Multiway zip join will produce an __entries annotation, which is an array where each element is a struct containing the __entries annotation (array of structs) for that sample
            merged = hl.Table.multi_way_zip_join(to_merge, "__entries", "__cols")
            # Flatten __entries while taking into account different entry lengths at different samples/variants (samples lacking a variant will be NA)
            merged = merged.annotate(
                __entries=hl.flatten(
                    hl.range(hl.len(merged.__entries)).map(
                        # Coalesce will return the first non-missing argument, so if the entry info is not missing, use that info, but if it is missing, create an entries struct with the correct element type for each null entry annotation (such as int32 for DP)
                        lambda i: hl.coalesce(
                            merged.__entries[i].__entries,
                            hl.range(hl.len(merged.__cols[i].__cols)).map(
                                lambda j: hl.null(
                                    merged.__entries.__entries.dtype.element_type.element_type
                                )
                            ),
                        )
                    )
                )
            )

            # Flatten col annotation from array<struct{__cols: array<struct{s: str}>} to array<struct{s: str}>
            merged = merged.annotate_globals(
                __cols=hl.flatten(merged.__cols.map(lambda x: x.__cols))
            )

            next_stage.append(
                merged.checkpoint(
                    os.path.join(temp_dir, f"stage_{stage}_job_{i}.ht"), overwrite=True
                )
            )
        info(f"Completed stage {stage}")
        stage += 1
        staging.clear()
        staging.extend(next_stage)

    # Unlocalize the entries, and unfilter the filtered entries and populate fields with missing values
    return (
        staging[0]
        ._unlocalize_entries("__entries", "__cols", list(mts[0].col_key))
        .unfilter_entries()
    )


def main(args):  # noqa: D103
    input_tsv = args.input_tsv
    output_ht = args.output_ht
    temp_dir = args.temp_dir
    chunk_size = args.chunk_size
    overwrite = args.overwrite

    if args.overwrite == False and hl.hadoop_exists(output_ht):
        logger.warn(
            "Overwrite is set to False but file already exists at %s, script will run but output will not be written",
            output_ht,
        )
    # Ensure that user supplied ht extension for output_ht
    if not output_ht.endswith(".ht"):
        sys.exit("Path supplied as output_ht must end with .ht extension")

    mt_list = []
    logger.info(
        "Reading in individual coverage files as matrix tables and adding to a list of matrix tables..."
    )
    with hl.hadoop_open(input_tsv, "r") as f:
        next(f)
        for line in f:
            line = line.rstrip()
            items = line.split("\t")
            participant_id, base_level_coverage_metrics, sample = items[0:3]
            mt = hl.import_matrix_table(
                base_level_coverage_metrics,
                delimiter="\t",
                row_fields={"chrom": hl.tstr, "pos": hl.tint, "target": hl.tstr},
                row_key=["chrom", "pos"],
            ).drop("target")
            mt = mt.rename({"x": "coverage"})
            mt = mt.key_cols_by(s=sample)
            mt_list.append(mt)

    logger.info("Joining individual coverage mts...")
    out_dir = dirname(output_ht)

    cov_mt = multi_way_union_mts(mt_list, temp_dir, chunk_size)
    n_samples = cov_mt.count_cols()

    logger.info("Adding coverage annotations...")
    # Calculate the mean and median coverage as well the fraction of samples above 100x or 1000x coverage at each base
    cov_mt = cov_mt.annotate_rows(
        locus=hl.locus(cov_mt.chrom, cov_mt.pos, reference_genome="GRCh38"),
        mean=hl.float(hl.agg.mean(cov_mt.coverage)),
        median=hl.median(hl.agg.collect(cov_mt.coverage)),
        over_100=hl.float((hl.agg.count_where(cov_mt.coverage > 100) / n_samples)),
        over_1000=hl.float((hl.agg.count_where(cov_mt.coverage > 1000) / n_samples)),
    )
    cov_mt.show()

    cov_mt = cov_mt.key_rows_by("locus").drop("chrom", "pos")

    output_mt = re.sub("\.ht$", ".mt", output_ht)
    output_tsv = re.sub("\.ht$", ".tsv", output_ht)
    output_samples = re.sub("\.ht$", "_sample_level.txt", output_ht)

    logger.info("Writing sample level coverage...")
    sample_mt = cov_mt.key_rows_by(pos=cov_mt.locus.position)
    sample_mt.coverage.export(output_samples)

    logger.info("Writing coverage mt and ht...")
    cov_mt.write(output_mt, overwrite=overwrite)
    cov_ht = cov_mt.rows()
    cov_ht = cov_ht.checkpoint(output_ht, overwrite=overwrite)
    cov_ht.export(output_tsv)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script combines individual mitochondria coverage files and outputs a hail table with coverage annotations"
    )
    parser.add_argument(
        "-i",
        "--input-tsv",
        help="Input file with coverage files to combine in tab-delimited format of participant_id, base_level_coverage_metrics, sample",
        required=True,
    )
    parser.add_argument(
        "-o", "--output-ht", help="Name of ht to write output", required=True
    )
    parser.add_argument(
        "-t",
        "--temp-dir",
        help="Temporary directory to use for intermediate outputs",
        required=True,
    )
    parser.add_argument(
        "--slack-token", help="Slack token that allows integration with slack",
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to",
    )
    parser.add_argument(
        "--chunk-size",
        help="Chunk size to use for combining VCFs (the number of individual VCFs that should be combined at a time)",
        type=int,
        default=100,
    )
    parser.add_argument(
        "--overwrite", help="Overwrites existing files", action="store_true"
    )

    args = parser.parse_args()

    # Both a slack token and slack channel must be supplied to receive notifications on slack
    if args.slack_channel and args.slack_token:
        with slack_notifications(args.slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
