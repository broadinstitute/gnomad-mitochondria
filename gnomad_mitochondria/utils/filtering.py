import hail as hl
import logging


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("add annotations")
logger.setLevel(logging.INFO)


def filter_by_copy_number(
    input_mt: hl.MatrixTable, keep_all_samples: bool = False
) -> hl.MatrixTable:
    """
    Calculate the mitochondrial copy number based on mean mitochondrial coverage and median nuclear coverage. Filter out samples with more extreme copy numbers.

    Note that median and mean coverage for mitochondria are very similar. Mean mitochondria coverage was used based on metrics available at the time, but releases will switch to using median mitochondria coverage.

    :param hl.MatrixTable input_mt: MatrixTable
    :param keep_all_samples: If True, keep all samples (calculate mitochondrial copy number, but do not filter any samples based on this metric)
    :return: MatrixTable filtered to samples with a copy number of at least 50 and less than 500, number samples below 50 removed, number samples above 500 removed
    """
    # Calculate mitochondrial copy number, if median autosomal coverage is not present default to a wgs_median_coverage of 30x
    input_mt = input_mt.annotate_cols(
        mito_cn=2
        * input_mt.mt_mean_coverage
        / hl.if_else(
            hl.is_missing(input_mt.wgs_median_coverage),
            30,
            input_mt.wgs_median_coverage,
        )
    )
    n_removed_below_cn = input_mt.aggregate_cols(
        hl.agg.count_where(input_mt.mito_cn < 50)
    )
    n_removed_above_cn = input_mt.aggregate_cols(
        hl.agg.count_where(input_mt.mito_cn > 500)
    )

    if not keep_all_samples:
        # Remove samples with a mitochondrial copy number below 50 or greater than 500
        input_mt = input_mt.filter_cols(
            (input_mt.mito_cn >= 50) & (input_mt.mito_cn <= 500)
        )
    input_mt = input_mt.filter_rows(hl.agg.any(input_mt.HL > 0))

    return input_mt, n_removed_below_cn, n_removed_above_cn


def filter_by_contamination(
    input_mt: hl.MatrixTable, output_dir: str, keep_all_samples: bool = False
) -> hl.MatrixTable:
    """
    Calculate contamination based on internal algorithm and filter out samples with contamination above 2%.

    Contamination takes into account:
    a) mitochondria contamination output by HaploCheck
    b) nuclear contamination (freemix) output by VerifyBamID
    c) an internal algorithm with utilizes the PASS haplogroup-defining variants which should be homoplasmic (100% alternate alleles), but in contaminated samples show multiple alleles with heteroplasmy 85-99.8%

    :param input_mt: MatrixTable
    :param output_dir: Output directory to which results should be written
    :param keep_all_samples: If True, keep all samples (calculate contamination, but do not filter any samples based on this metric)
    :return: MatrixTable filtered to samples without contamination, number of contaminated samples removed
    """
    # Generate expression for genotypes with >= 85% heteroplasmy and no FT filters at haplogroup-defining sites that are not filtered as artifact-prone sites
    over_85_expr = (
        (input_mt.HL >= 0.85)
        & (input_mt.FT == {"PASS"})
        & input_mt.hap_defining_variant
        & ~hl.str(input_mt.filters).contains("artifact_prone_site")
    )

    input_mt = input_mt.annotate_cols(
        over_85_mean=hl.agg.filter(over_85_expr, hl.agg.mean(input_mt.HL)),
        over_85_count=hl.agg.filter(
            over_85_expr, hl.agg.count_where(hl.is_defined(input_mt.HL))
        ),
        bt_85_and_99_mean=hl.agg.filter(
            over_85_expr & (input_mt.HL <= 0.998), hl.agg.mean(input_mt.HL)
        ),
        bt_85_and_99_count=hl.agg.filter(
            over_85_expr & (input_mt.HL <= 0.998),
            hl.agg.count_where(hl.is_defined(input_mt.HL)),
        ),
    )

    input_mt = input_mt.annotate_cols(
        contam_high_het=hl.if_else(
            input_mt.bt_85_and_99_count >= 3,
            1 - input_mt.bt_85_and_99_mean,
            1 - input_mt.over_85_mean,
        )
    )

    # If contam_high_het is nan, set to 0 (to avoid filtering out missing values which would be more common with haplogroups closer to the reference haplogroup)
    input_mt = input_mt.annotate_cols(
        contam_high_het=hl.if_else(
            hl.is_nan(input_mt.contam_high_het), 0, input_mt.contam_high_het
        )
    )

    # Find samples on border of .02 that may flip between < 0.02 and > 0.02 from issues with floating point precision and mark these samples for removal
    epsilon = 0.000001
    border_samples = input_mt.aggregate_cols(
        hl.agg.filter(
            (input_mt.contam_high_het > (0.02 - epsilon))
            & (input_mt.contam_high_het < (0.02 + epsilon)),
            hl.agg.collect((input_mt.s)),
        )
    )

    border_samples = (
        hl.literal(border_samples) if border_samples else hl.empty_array(hl.tstr)
    )

    # Add annotation to keep only samples with a contamination less than 2%
    input_mt = input_mt.annotate_cols(
        keep=(input_mt.contamination < 0.02)
        & (input_mt.freemix_percentage < 2)
        & (input_mt.contam_high_het < 0.02)
        & ~border_samples.contains(input_mt.s)
    )
    # Save sample contamination information to separate file
    n_contaminated = input_mt.aggregate_cols(hl.agg.count_where(~input_mt.keep))

    sample_data = input_mt.select_cols(
        "contamination",
        "freemix_percentage",
        "contam_high_het",
        "over_85_mean",
        "over_85_count",
        "bt_85_and_99_mean",
        "bt_85_and_99_count",
        "keep",
    )
    data_export = sample_data.cols()
    data_export.export(f"{output_dir}/sample_contamination.tsv")

    if not keep_all_samples:
        logger.info(
            "Removing %d samples with contamination above 2 percent", n_contaminated
        )
        input_mt = input_mt.filter_cols(input_mt.keep)
    input_mt = input_mt.drop("keep")

    input_mt = input_mt.filter_rows(hl.agg.any(input_mt.HL > 0))

    return input_mt, n_contaminated


def filter_genotypes_below_min_het_threshold(
    input_mt: hl.MatrixTable, min_het_threshold: float = 0.10
) -> hl.MatrixTable:
    """
    Filter out genotypes with a heteroplasmy below the min_het_threshold.

    This filter is a genotype level filter to remove variants with a heteroplasmy level below the specified min_het_threshold
    NOTE: Should later parameterize this function to allow other heteroplasmy cutoffs?

    :param input_mt: MatrixTable
    :param min_het_threshold: Minimum heteroplasmy level to define a variant as a PASS heteroplasmic variant, genotypes below this threshold will count towards the heteroplasmy_below_min_het_threshold filter and be set to missing
    :return: MatrixTable with the heteroplasmy_below_min_het_threshold in the FT field added where applicable
    """
    input_mt = input_mt.annotate_entries(
        FT=hl.if_else(
            (input_mt.HL < min_het_threshold) & (input_mt.GT.is_het()),
            input_mt.FT.add("heteroplasmy_below_min_het_threshold"),
            input_mt.FT,
        )
    )

    # Remove "PASS" from FT if it's not the only filter
    input_mt = input_mt.annotate_entries(
        FT=hl.if_else(input_mt.FT != {"PASS"}, input_mt.FT.remove("PASS"), input_mt.FT)
    )

    return input_mt


def filter_genotypes(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Set all genotype field values to missing if the variant is not "PASS" for that sample.

    :param input_mt: MatrixTable
    :return: MatrixTable with filtered genotype fields set to missing
    """
    pass_expr = input_mt.FT == {"PASS"}

    input_mt = input_mt.annotate_entries(
        GT=hl.or_missing(pass_expr, input_mt.GT),
        DP=hl.or_missing(pass_expr, input_mt.DP),
        HL=hl.or_missing(pass_expr, input_mt.HL),
        FT=hl.or_missing(pass_expr, input_mt.FT),
        MQ=hl.or_missing(pass_expr, input_mt.MQ),
        TLOD=hl.or_missing(pass_expr, input_mt.TLOD),
    )

    return input_mt


def remove_low_allele_frac_genotypes(
    input_mt: hl.MatrixTable, vaf_filter_threshold: float = 0.01
) -> hl.MatrixTable:
    """
    Remove low_allele_frac genotypes and sets the call to homoplasmic reference.

    NOTE: vaf_filter_threshold should match what was supplied to the vaf_filter_threshold when running Mutect2, variants below this value will be set to homoplasmic reference after calculating the common_low_heteroplasmy filter, Mutect2 will have flagged these variants as "low_allele_frac"

    :param input_mt: MatrixTable
    :param vaf_filter_threshold: Should match vaf_filter_threshold supplied to Mutect2, variants below this value will be set to homoplasmic reference after calculating the common_low_heteroplasmy filter
    :return: MatrixTable with genotypes below the vaf_filter_threshold set to homoplasmic reference
    """
    # Set HL to 0 if < vaf_filter_threshold and remove variants that no longer have at least one alt call
    input_mt = input_mt.annotate_entries(
        HL=hl.if_else(
            (input_mt.HL > 0) & (input_mt.HL < vaf_filter_threshold), 0, input_mt.HL
        )
    )
    # Filter field for all variants with a heteroplasmy of 0 should be set to PASS
    # This step is needed to prevent homref calls that are filtered
    input_mt = input_mt.annotate_entries(
        FT=hl.if_else(input_mt.HL < vaf_filter_threshold, {"PASS"}, input_mt.FT)
    )
    input_mt = input_mt.annotate_entries(
        GT=hl.if_else(
            input_mt.HL < vaf_filter_threshold, hl.parse_call("0/0"), input_mt.GT
        )
    )

    # Check that variants no longer contain the "low_allele_frac" filter (vaf_filter_threshold should be set to appropriate level to remove these variants)
    laf_rows = input_mt.filter_rows(
        hl.agg.any(hl.str(input_mt.FT).contains("low_allele_frac"))
    )
    n_laf_rows = laf_rows.count_rows()
    if n_laf_rows > 0:
        sys.exit(
            "low_allele_frac filter should no longer be present after applying vaf_filter_threshold (vaf_filter_threshold should equal the vaf_filter_threshold supplied to Mutect2)"
        )
    input_mt = input_mt.filter_rows(hl.agg.any(input_mt.HL > 0))

    return input_mt


def apply_indel_stack_filter(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Apply the indel_stack filter to the MatrixTable.

    The indel_stack filter marks alleles where all samples with the variant call had at least 2 different indels called at the position

    :param input_mt: MatrixTable
    :return: MatrixTable with the indel_stack filter added
    """
    # Add variant-level indel_stack at any indel allele where all samples with a variant call had at least 2 different indels called at that position
    # If any sample had a solo indel at that position, do not filter
    indel_expr = get_indel_expr(input_mt)
    input_mt = input_mt.annotate_cols(
        indel_pos_counter=hl.agg.filter(
            indel_expr, hl.agg.counter(input_mt.locus.position)
        )
    )
    indel_expr = get_indel_expr(input_mt)
    input_mt = input_mt.annotate_entries(
        indel_occurences=(
            hl.case()
            .when(
                (
                    indel_expr
                    & (input_mt.indel_pos_counter.get(input_mt.locus.position) >= 2)
                ),
                "stack",
            )
            .when(
                (
                    indel_expr
                    & (input_mt.indel_pos_counter.get(input_mt.locus.position) == 1)
                ),
                "solo",
            )
            .or_missing()
        )
    )

    # If stack is true and solo is false, the indel is stack only and should be filtered out
    input_mt = input_mt.annotate_rows(
        filters=hl.if_else(
            hl.agg.any(input_mt.indel_occurences == "stack")
            & ~hl.agg.any(input_mt.indel_occurences == "solo"),
            input_mt.filters.add("indel_stack"),
            input_mt.filters,
        )
    )

    return input_mt


def apply_common_low_het_flag(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Apply the common_low_heteroplasmy flag to the MatrixTable.

    The common_low_heteroplasmy flag marks variants where the overall frequency is > 0.001 for samples with a heteroplasmy level > 0 and < 0.50 and either "low_allele_frac" or "PASS" for the genotype filter

    NOTE: The "low_allele_frac" is applied by Mutect2 to variants with a heteroplasmy level below the supplied vaf_filter_threshold

    :param input_mt: MatrixTable
    :return: MatrixTable with the common_low_heteroplasmy flag added
    """
    input_mt = input_mt.annotate_rows(
        AC_mid_het=hl.agg.count_where(
            (input_mt.HL < 0.50)
            & (input_mt.HL > 0.0)
            & ((input_mt.FT == {"PASS"}) | (input_mt.FT == {"low_allele_frac"}))
        )
    )
    input_mt = input_mt.annotate_rows(
        AF_mid_het=input_mt.AC_mid_het
        / hl.agg.count_where(
            hl.is_defined(input_mt.HL)
            & ((input_mt.FT == {"PASS"}) | (input_mt.FT == {"low_allele_frac"}))
        )
    )
    input_mt = input_mt.annotate_rows(
        common_low_heteroplasmy=input_mt.AF_mid_het > 0.001
    )

    return input_mt


def apply_npg_filter(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Apply the npg filter to the MatrixTable.

    The npg (no pass genotypes) filter marks sites that don't have at least one pass alt call

    :param input_mt: MatrixTable
    :return: MatrixTable with the npg filter added
    """
    input_mt = input_mt.annotate_rows(
        filters=hl.if_else(
            ~(hl.agg.any((input_mt.HL > 0.0) & (input_mt.FT == {"PASS"}))),
            input_mt.filters.add("npg"),
            input_mt.filters,
        )
    )

    return input_mt
