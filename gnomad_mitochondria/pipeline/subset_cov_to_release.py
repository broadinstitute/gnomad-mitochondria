#!/usr/bin/env python
import argparse
import logging
import re

import hail as hl

from gnomad.utils.slack import slack_notifications


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("subset cov")
logger.setLevel(logging.INFO)


def main(args):  # noqa: D103
    input_mt_path = args.input_mt_path
    cov_mt_path = args.cov_mt_path
    out_tsv_path = args.out_tsv_path

    logger.info("Subsetting coverage mt to samples in the input mt...")
    input_mt = hl.read_matrix_table(input_mt_path)
    cov_mt = hl.read_matrix_table(cov_mt_path)
    samples_to_subset = hl.literal(input_mt.s.collect())
    cov_mt = cov_mt.filter_cols(samples_to_subset.contains(cov_mt["s"]))

    n_samples_input = input_mt.count_cols()
    n_samples_cov = cov_mt.count_cols()
    logger.info(
        "Input mt has %d samples and subsetted coverage table has %d samples",
        n_samples_input,
        n_samples_cov,
    )

    # Calculate the mean and median coverage as well the fraction of samples above 100x or 1000x coverage at each base
    cov_mt = cov_mt.annotate_rows(
        mean=hl.float(hl.agg.mean(cov_mt.coverage)),
        median=hl.median(hl.agg.collect(cov_mt.coverage)),
        over_100=hl.float((hl.agg.count_where(cov_mt.coverage > 100) / n_samples_cov)),
        over_1000=hl.float(
            (hl.agg.count_where(cov_mt.coverage > 1000) / n_samples_cov)
        ),
    )
    cov_ht = cov_mt.rows()

    output_ht = re.sub(r"\.tsv$", ".ht", out_tsv_path)
    cov_ht = cov_ht.checkpoint(output_ht, overwrite=True)

    logger.info("Writing results to tsv...")
    cov_ht.export(out_tsv_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script subsets the coverage files to the sample in the supplied mt"
    )
    parser.add_argument(
        "-i",
        "--input-mt-path",
        help='Path to MatrixTable containing a sample "s" column (coverage mt will be subset to the samples in this mt)',
        required=True,
    )
    parser.add_argument(
        "-c",
        "--cov-mt-path",
        help="Path to MatrixTable of per base coverages for all samples (per sample and per variant",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--out-tsv-path",
        help="Path to which resulting tsv should be written",
        required=True,
    )
    parser.add_argument(
        "--slack-token", help="Slack token that allows integration with slack",
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to",
    )

    args = parser.parse_args()

    # Both a slack token and slack channel must be supplied to receive notifications on slack
    if args.slack_channel and args.slack_token:
        with slack_notifications(args.slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
