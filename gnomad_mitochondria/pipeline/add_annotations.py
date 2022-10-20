#!/usr/bin/env python
import argparse
import hail as hl
import logging

from gnomad.utils.annotations import age_hists_expr
from gnomad.utils.reference_genome import add_reference_sequence
from gnomad.utils.slack import slack_notifications
from gnomad_qc.v3.resources.meta import meta  # pylint: disable=import-error
from gnomad_mitochondria.pipeline.annotation_descriptions import (
    add_descriptions,
    adjust_descriptions,
)
from gnomad_mitochondria.utils.annotations import (
    add_age_and_pop,
    add_annotations_by_hap_and_pop,
    add_filter_annotations,
    add_genotype,
    add_gnomad_metadata,
    add_hap_defining,
    add_quality_histograms,
    add_rsids,
    add_sample_annotations,
    add_terra_metadata,
    add_trna_predictions,
    add_variant_context,
    add_vep,
    generate_expressions,
)
from gnomad_mitochondria.utils.exports import (
    export_simplified_variants,
    format_vcf,
    generate_output_paths,
    report_stats,
)
from gnomad_mitochondria.utils.filtering import (
    filter_by_contamination,
    filter_by_copy_number,
    filter_genotypes,
)

# Github repo locations for imports:
# gnomad: https://github.com/broadinstitute/gnomad_methods
# gnomad_qc: https://github.com/broadinstitute/gnomad_qc


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("add annotations")
logger.setLevel(logging.INFO)

logger.info("Setting hail flag to avoid array index out of bounds error...")
# Setting this flag isn't generally recommended, but is needed (since at least Hail version 0.2.75) to avoid an array index out of bounds error until changes are made in future versions of Hail
# TODO: reassess if this flag is still needed for future versions of Hail
hl._set_flags(no_whole_stage_codegen="1")


def main(args):  # noqa: D103
    mt_path = args.mt_path
    output_dir = args.output_dir
    participant_data = args.participant_data
    vep_results = args.vep_results
    min_hom_threshold = args.min_hom_threshold
    vaf_filter_threshold = args.vaf_filter_threshold
    min_het_threshold = args.min_het_threshold
    gnomad_subset = args.subset_to_gnomad_release
    keep_all_samples = args.keep_all_samples
    run_vep = args.run_vep

    logger.info("Cutoff for homoplasmic variants is set to %.2f...", min_hom_threshold)

    # Define mt path, output directory, subset name
    subset_name = ""

    logger.info("Adding genotype annotation...")
    mt = add_genotype(mt_path, min_hom_threshold)

    logger.info("Adding annotations from Terra...")
    mt = add_terra_metadata(mt, participant_data)

    logger.info("Annotating haplogroup-defining variants...")
    mt = add_hap_defining(mt)

    logger.info("Annotating tRNA predictions...")
    mt = add_trna_predictions(mt)

    # If 'subset-to-gnomad-release' is set, 'age' and 'pop' are added by the add_gnomad_metadata function.
    # If 'subset-to-gnomad-release' is not set, the user should include an 'age' and 'pop' column in the file supplied to `participant-data`.
    if gnomad_subset:
        logger.info("Adding gnomAD metadata sample annotations...")
        mt = add_gnomad_metadata(mt)
    else:
        logger.info("Adding age and pop annotations...")
        mt = add_age_and_pop(mt, participant_data)

    logger.info("Adding variant context annotations...")
    mt = add_variant_context(mt)

    # If specified, subet to only the gnomAD samples in the current release
    if gnomad_subset:
        logger.warning("Subsetting results to gnomAD release samples...")
        subset_name = "_gnomad"

        # Subset to release samples and filter out rows that no longer have at least one alt call
        mt = mt.filter_cols(mt.release)  # Filter to cols where release is true
        mt = mt.filter_rows(hl.agg.any(mt.HL > 0))

    logger.info("Checking for samples with low/high mitochondrial copy number...")
    mt, n_removed_below_cn, n_removed_above_cn = filter_by_copy_number(
        mt, keep_all_samples
    )

    logger.info("Checking for contaminated samples...")
    mt, n_contaminated = filter_by_contamination(mt, output_dir, keep_all_samples)

    logger.info("Switch build and checkpoint...")
    # Switch build 37 to build 38
    mt = mt.key_rows_by(
        locus=hl.locus("chrM", mt.locus.position, reference_genome="GRCh38"),
        alleles=mt.alleles,
    )
    mt = mt.checkpoint(f"{output_dir}/prior_to_vep.mt", overwrite=args.overwrite)

    logger.info("Adding vep annotations...")
    mt = add_vep(mt, run_vep, vep_results)

    logger.info("Adding dbsnp annotations...")
    mt = add_rsids(mt)

    logger.info("Setting up output paths...")
    annotated_mt_path = generate_output_paths(
        output_dir, "annotated_combined", subset_name, "mt"
    )
    sites_ht_path = generate_output_paths(
        output_dir, "combined_sites_only", subset_name, "ht"
    )
    sites_txt_path = generate_output_paths(
        output_dir, "combined_sites_only", subset_name, "txt"
    )
    sites_vcf_path = generate_output_paths(
        output_dir, "combined_sites_only", subset_name, "vcf.bgz"
    )
    samples_txt_path = generate_output_paths(
        output_dir, "sample_annotations", subset_name, "txt"
    )
    samples_vcf_path = generate_output_paths(
        output_dir, "sample_vcf", subset_name, "vcf.bgz"
    )

    logger.info("Results will be output to the following files:")
    print(
        "\n".join(
            [
                annotated_mt_path,
                sites_ht_path,
                sites_txt_path,
                sites_vcf_path,
                samples_txt_path,
                samples_vcf_path,
            ]
        )
    )

    logger.info("Annotating MT...")
    mt, n_het_below_min_het_threshold = add_filter_annotations(
        mt, vaf_filter_threshold, min_het_threshold
    )

    mt = mt.checkpoint(
        f"{output_dir}/prior_to_filter_genotypes.mt", overwrite=args.overwrite
    )

    mt = filter_genotypes(mt)
    # Add variant annotations such as AC, AF, and AN
    mt = mt.annotate_rows(**dict(generate_expressions(mt, min_hom_threshold)))
    # Checkpoint to help avoid Hail errors from large queries
    mt = mt.checkpoint(f"{output_dir}/temp.mt", overwrite=args.overwrite)
    mt = add_quality_histograms(mt)
    mt = add_annotations_by_hap_and_pop(mt)
    mt = add_descriptions(
        mt, min_hom_threshold, vaf_filter_threshold, min_het_threshold
    )
    mt = mt.checkpoint(
        annotated_mt_path, overwrite=args.overwrite
    )  # Full matrix table for internal use

    logger.info("Generating summary statistics reports...")
    report_stats(
        mt,
        output_dir,
        False,
        n_removed_below_cn,
        n_removed_above_cn,
        n_contaminated,
        n_het_below_min_het_threshold,
    )
    report_stats(
        mt,
        output_dir,
        True,
        n_removed_below_cn,
        n_removed_above_cn,
        n_contaminated,
        n_het_below_min_het_threshold,
    )

    logger.info("Writing ht...")
    variant_ht = mt.rows()
    variant_ht = variant_ht.drop("region", "variant_context")
    variant_ht = adjust_descriptions(variant_ht)
    variant_ht.export(sites_txt_path)  # Sites-only txt file for external use
    variant_ht.write(
        sites_ht_path, overwrite=args.overwrite
    )  # Sites-only ht for external use

    logger.info("Writing sample annotations...")
    mt = add_sample_annotations(mt, min_hom_threshold)
    sample_ht = mt.cols()
    sample_ht.group_by(sample_ht.hap).aggregate(n=hl.agg.count()).export(
        f"{output_dir}/haplogroup_counts.txt"
    )  # Counts of top level haplogroups
    sample_ht.export(samples_txt_path)  # Sample annotations txt file for internal use

    logger.info("Formatting and writing VCF...")
    rows_ht = mt.rows()
    export_simplified_variants(rows_ht, output_dir)
    vcf_mt, vcf_meta, vcf_header_file = format_vcf(mt, output_dir, min_hom_threshold)
    hl.export_vcf(
        vcf_mt,
        samples_vcf_path,
        metadata=vcf_meta,
        append_to_header=vcf_header_file,
        tabix=True,
    )  # Full VCF for internal use
    vcf_variant_ht = vcf_mt.rows()
    rows_mt = hl.MatrixTable.from_rows_table(vcf_variant_ht).key_cols_by(s="foo")
    hl.export_vcf(
        rows_mt,
        sites_vcf_path,
        metadata=vcf_meta,
        append_to_header=vcf_header_file,
        tabix=True,
    )  # Sites-only VCF for external use

    logger.info("All annotation steps are completed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script adds variant annotations to the mitochondria VCF/MT"
    )
    parser.add_argument("-m", "--mt-path", help="Path to combined mt", required=True)
    parser.add_argument(
        "-d",
        "--output-dir",
        help="Path to directory to which output should be written",
        required=True,
    )
    parser.add_argument(
        "-a",
        "--participant-data",
        help="Output file that results from Terra data download",
        required=True,
    )
    parser.add_argument(
        "-v",
        "--vep-results",
        help="MatrixTable path to output vep results (either the existing results or where to ouput new vep results if also setting run_vep)",
        required=True,
    )
    parser.add_argument(
        "--slack-token",
        help="Slack token that allows integration with slack",
    )
    parser.add_argument(
        "--slack-channel",
        help="Slack channel to post results and notifications to",
    )
    parser.add_argument(
        "--min-het-threshold",
        help="Minimum heteroplasmy level to define a variant as a PASS heteroplasmic variant, genotypes below this threshold will count towards the heteroplasmy_below_min_het_threshold filter and be set to missing",
        type=float,
        default=0.10,
    )
    parser.add_argument(
        "--min-hom-threshold",
        help="Minimum heteroplasmy level to define a variant as homoplasmic",
        type=float,
        default=0.95,
    )
    parser.add_argument(
        "--vaf-filter-threshold",
        help="Should match vaf_filter_threshold supplied to Mutect2, variants below this value will be set to homoplasmic reference after calculating the common_low_heteroplasmy filter",
        type=float,
        default=0.01,
    )
    parser.add_argument(
        "--subset-to-gnomad-release",
        help="Set to True to only include released gnomAD samples",
        action="store_true",
    )
    parser.add_argument(
        "--keep-all-samples",
        help="Set to True to keep all samples (will skip steps that filter samples because of contamination and/or mitochondrial copy number)",
        action="store_true",
    )
    parser.add_argument(
        "--run-vep", help="Set to True to run/rerun vep", action="store_true"
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
