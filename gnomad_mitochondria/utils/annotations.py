import hail as hl
import logging
import re
import sys

from gnomad.resources.grch38.gnomad import POPS
from gnomad.resources.grch38.reference_data import dbsnp, _import_dbsnp
from gnomad.utils.annotations import age_hists_expr
from gnomad.utils.reference_genome import add_reference_sequence
from gnomad_qc.v3.resources.meta import meta  # pylint: disable=import-error


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("add annotations")
logger.setLevel(logging.INFO)

# Include NA in POPS to account for cases where population annotations are missing
POPS.append("NA")

RESOURCES = {
    "variant_context": "gs://gcp-public-data--gnomad/resources/mitochondria/variant_context/chrM_pos_ref_alt_context_categories.txt",
    "phylotree": "gs://gcp-public-data--gnomad/resources/mitochondria/phylotree/rCRS-centered_phylo_vars_final_update.txt",
    "pon_mt_trna": "gs://gcp-public-data--gnomad/resources/mitochondria/trna_predictions/pon_mt_trna_predictions_08_27_2020.txt",
    "mitotip": "gs://gcp-public-data--gnomad/resources/mitochondria/trna_predictions/mitotip_scores_08_27_2020.txt",
}


def add_genotype(mt_path: str, min_hom_threshold: float = 0.95) -> hl.MatrixTable:
    """
    Add in genotype annotation based on heteroplasmy level.

    If the heteroplasmy level is above the min_hom_threshold, set the genotype to 1/1.
    If the heteroplasmy level is less than the min_hom_threshold, but greater than 0, set the genotype to 0/1.
    Otherwise set the genotype to 0/0.

    :param mt_path: Path to the MatrixTable (this MatrixTable can be generated by running combine_vcfs.py)
    :param min_hom_threshold: Minimum heteroplasmy level to define a variant as homoplasmic
    :return: MatrixTable with GT field added
    """
    logger.info("Reading in MT...")
    mt = hl.read_matrix_table(mt_path)

    # Add in genotype (GT) based on min_hom_threshold
    mt = mt.annotate_entries(
        GT=(
            hl.case()
            .when((mt.HL < min_hom_threshold) & (mt.HL > 0.0), hl.parse_call("0/1"))
            .when(mt.HL >= min_hom_threshold, hl.parse_call("1/1"))
            .when(mt.HL == 0, hl.parse_call("0/0"))
            .default(hl.missing(hl.tcall))
        ),
    )

    return mt


def add_variant_context(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Add variant context annotations to the MatrixTable.

    This fucntion adds in information on regions/strand for SNPs that can be useful for determining mutational signatures.

    :param input_mt: MatrixTable
    :return: MatrixTable with variant context information added
    """
    # Read in variant context data
    vc_ht = hl.import_table(RESOURCES["variant_context"], impute=True)

    # Split columns into separate annotations
    vc_ht = vc_ht.annotate(
        ref=vc_ht["POS.REF.ALT"].split(r"\.")[1],
        alt=vc_ht["POS.REF.ALT"].split(r"\.")[2],
        strand=vc_ht.Context_category.split("_")[-1],
        variant=vc_ht.Context_category.split("_")[0],
    )

    # Rename and select certain columns
    vc_ht = vc_ht.rename({"MT_POS": "pos", "Annotation": "region"})
    vc_ht = vc_ht.select("pos", "ref", "alt", "strand", "region", "variant")

    # Key by locus and allele
    vc_ht = vc_ht.key_by(
        locus=hl.locus("MT", vc_ht.pos, reference_genome="GRCh37"),
        alleles=[vc_ht.ref, vc_ht.alt],
    )

    # Annotate original mt with variant context information
    input_mt = input_mt.annotate_rows(**vc_ht[input_mt.locus, input_mt.alleles])
    input_mt = input_mt.annotate_rows(
        variant_context=hl.str(input_mt.variant) + "_" + hl.str(input_mt.strand)
    )
    input_mt = input_mt.drop("pos", "ref", "alt", "strand", "variant")

    return input_mt


def add_gnomad_metadata(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Add select gnomAD metadata to the MatrixTable.

    :param input_mt: MatrixTable
    :return: MatrixTable with select gnomAD metadata added
    """
    # TODO: Add option here to accomodate non-gnomAD metadata
    genome_meta_ht = meta.versions["3.1"].ht()

    genome_meta_struct = genome_meta_ht[input_mt.s]

    input_mt = input_mt.annotate_cols(
        release=genome_meta_struct.release,
        hard_filters=genome_meta_struct.sample_filters.hard_filters,
        research_project=genome_meta_struct.project_meta.research_project,
        project_id=genome_meta_struct.project_meta.project_id,
        product=genome_meta_struct.project_meta.product,
        sample_pi=genome_meta_struct.project_meta.sample_pi,
        sex_karyotype=genome_meta_struct.sex_imputation.sex_karyotype,
        age=hl.if_else(
            hl.is_defined(genome_meta_struct.project_meta.age),
            genome_meta_struct.project_meta.age,
            genome_meta_struct.project_meta.age_alt,
        ),
        broad_external=genome_meta_struct.project_meta.broad_external,
        pop=genome_meta_struct.population_inference.pop,
    )

    return input_mt


def add_age_and_pop(input_mt: hl.MatrixTable, participant_data: str) -> hl.MatrixTable:
    """
    Add sample-level metadata for age and pop to `input_mt`.

    :param input_mt: MatrixTable
    :param participant_data: Path to metadata file downloaded from Terra that contains sample age and pop information
    :return: MatrixTable with select age and pop annotations added
    """
    ht = hl.import_table(
        participant_data,
        types={"age": hl.tint32, "pop": hl.tstr},
    ).key_by("s")

    ht = ht.select("age", "pop")

    input_mt = input_mt.annotate_cols(**ht[input_mt.col_key])

    # If a sample doesn't have an annotated population, set it to the string "NA"
    input_mt = input_mt.annotate_cols(
        pop=hl.if_else(hl.is_missing(input_mt.pop), "NA", input_mt.pop)
    )

    return input_mt


def add_terra_metadata(
    input_mt: hl.MatrixTable, participant_data: str
) -> hl.MatrixTable:
    """
    Add Terra metadata to the MatrixTable.

    The participant_data file can be obtained by downloading the participant data after running Mutect2 in Terra. This file should contain the following columns:
        - entity:participant_id: Participant ID uploaded to Terra by user
        - s: Sample ID uploaded to Terra by user
        - contamination: Output by Mutect2, gives the estimate of mitochondrial contamination
        - freemix_percentage: Uploaded to Terra by user, can be calculated with VerifyBamID
        - major_haplogroup: Output by Mutect2 which utilizes Haplogrep
        - wgs_median_coverage: Uploaded to Terra by user, can be calculated with Picard's CollectWgsMetrics
        - mt_mean_coverage: Output by Mutect2, gives the mean mitochondrial coverage

    :param input_mt: MatrixTable
    :param participant_data: Path to metadata file downloaded from Terra
    :return: MatrixTable with Terra metadata annotations added
    """
    # Add haplogroup and Mutect2/Terra output annotations
    ht = hl.import_table(
        participant_data,
        types={
            "contamination": hl.tfloat64,
            "freemix_percentage": hl.tfloat64,
            "mt_mean_coverage": hl.tfloat64,
            "wgs_median_coverage": hl.tfloat64,
        },
        missing="",
    ).key_by("s")
    ht = ht.rename({"entity:participant_id": "participant_id"})

    ht = ht.select(
        "participant_id",
        "contamination",
        "freemix_percentage",
        "major_haplogroup",
        "wgs_median_coverage",
        "mt_mean_coverage",
    )

    input_mt = input_mt.annotate_cols(**ht[input_mt.s])

    # Annotate the high level haplogroup by taking the first letter, with the exception of H and L haplogroups which are more commonly referred to using the first two letters
    input_mt = input_mt.annotate_cols(
        hap=hl.if_else(
            input_mt.major_haplogroup.startswith("HV")
            | input_mt.major_haplogroup.startswith("L"),
            input_mt.major_haplogroup[0:2],
            input_mt.major_haplogroup[0],
        )
    )

    return input_mt


def add_hap_defining(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Add bool on whether or not a variant is a haplogroup-defining variant to the MatrixTable.

    Haplogroup-defining annotations were obtained from PhyloTree Build 17.

    :param input_mt: MatrixTable
    :return: MatrixTable with annotation on whether or not the variant is haplogroup-defining added
    """
    # TODO: move dataset location
    hap_defining_variants = hl.import_table(RESOURCES["phylotree"])

    hap_defining = hl.literal(set(hap_defining_variants.variant.collect()))
    input_mt = input_mt.annotate_rows(
        variant_collapsed=input_mt.alleles[0]
        + hl.str(input_mt.locus.position)
        + input_mt.alleles[1]
    )
    input_mt = input_mt.annotate_rows(
        hap_defining_variant=hap_defining.contains(input_mt.variant_collapsed)
    )  # set hap_defining_variant to True or False

    return input_mt


def add_trna_predictions(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Add tRNA predictions on pathogenicity from PON-mt-tRNA and MitoTIP to the MatrixTable.

    :param input_mt: MatrixTable
    :return: MatrixTable with tRNA predictions of pathogenicity added
    """
    # Add PON-mt-tRNA predictions
    pon_predictions = hl.import_table(RESOURCES["pon_mt_trna"])

    # If reference allele from fasta doesn't match Reference_nucleotide, PON-mt-tRNA is reporting the allele of opposite strand and need to get reverse complement for ref and alt
    add_reference_sequence(hl.get_reference("GRCh37"))
    pon_predictions = pon_predictions.annotate(
        ref=hl.get_sequence(
            "MT", hl.int(pon_predictions.mtDNA_position), reference_genome="GRCh37"
        )
    )
    pon_predictions = pon_predictions.annotate(
        alt=hl.if_else(
            pon_predictions.Reference_nucleotide == pon_predictions.ref,
            pon_predictions.New_nucleotide,
            hl.reverse_complement(pon_predictions.New_nucleotide),
        )
    )
    pon_predictions = pon_predictions.key_by(
        variant_id=pon_predictions.ref
        + hl.str(pon_predictions.mtDNA_position)
        + pon_predictions.alt
    )
    input_mt = input_mt.annotate_rows(
        pon_mt_trna_prediction=pon_predictions[input_mt.variant_collapsed]
        .Classification.lower()
        .replace(" ", "_"),
        pon_ml_probability_of_pathogenicity=hl.float(
            pon_predictions[input_mt.variant_collapsed].ML_probability_of_pathogenicity
        ),
    )

    # Add MitoTIP predictions
    mitotip_predictions = hl.import_table(RESOURCES["mitotip"])
    mitotip_predictions = mitotip_predictions.key_by(
        variant_id=mitotip_predictions.rCRS
        + hl.str(mitotip_predictions.Position)
        + mitotip_predictions.Alt
    )
    input_mt = input_mt.annotate_rows(
        mitotip_score=hl.float(
            mitotip_predictions[input_mt.variant_collapsed].MitoTIP_Score
        )
    )
    # Set pathogenicity based on MitoTIP scores, classifications obtained from MitoTIP's website
    input_mt = input_mt.annotate_rows(
        mitotip_trna_prediction=(
            hl.case()
            .when(input_mt.mitotip_score > 16.25, "likely_pathogenic")
            .when(
                (input_mt.mitotip_score <= 16.25) & (input_mt.mitotip_score > 12.66),
                "possibly_pathogenic",
            )
            .when(
                (input_mt.mitotip_score <= 12.66) & (input_mt.mitotip_score >= 8.44),
                "possibly_benign",
            )
            .when((input_mt.mitotip_score < 8.44), "likely_benign")
            .or_missing()
        )
    )

    return input_mt


def generate_expressions(
    input_mt: hl.MatrixTable, min_hom_threshold: float = 0.95
) -> hl.MatrixTable:
    """
    Create expressions to use for annotating the MatrixTable.

    The expressions include AC, AN, AF, filtering allele frequency (FAF) split by homplasmic/heteroplasmic, haplgroup, and population.
    Also includes calcuations of mean DP, MQ, and TLOD.

    :param input_mt: MatrixTable
    :param min_hom_threshold: Minimum heteroplasmy level to define a variant as homoplasmic
    :return: Tuple of hail expressions
    """
    # Calculate AC and AN
    AC = hl.agg.count_where((input_mt.HL > 0.0))
    AN = hl.agg.count_where(hl.is_defined(input_mt.HL))
    # Note: if AN is zero, AFs will evaluate to NaN, which may need to be converted to zero for downstream tools
    AF = AC / AN

    # Calculate AC for het and hom variants, and histogram for HL
    AC_hom = hl.agg.count_where(input_mt.HL >= min_hom_threshold)
    AC_het = hl.agg.count_where((input_mt.HL < min_hom_threshold) & (input_mt.HL > 0.0))
    HL_hist = hl.agg.filter(input_mt.HL > 0, hl.agg.hist(input_mt.HL, 0, 1, 10))
    DP_hist_alt = hl.agg.filter(
        input_mt.GT.is_non_ref(), hl.agg.hist(input_mt.DP, 0, 2000, 10)
    )
    DP_hist_all = hl.agg.hist(input_mt.DP, 0, 2000, 10)
    DP_mean = hl.agg.mean(input_mt.DP)
    MQ_mean = hl.agg.mean(input_mt.MQ)
    TLOD_mean = hl.agg.mean(input_mt.TLOD)

    # Calculate AF
    # Note: if AN is zero, AFs will evaluate to NaN, which may need to be converted to zero for downstream tools
    AF_hom = AC_hom / AN
    AF_het = AC_het / AN

    # Calculate max individual heteroplasmy
    max_HL = hl.agg.max(input_mt.HL)

    # Haplogroup annotations
    pre_hap_AC = hl.agg.group_by(input_mt.hap, AC)
    pre_hap_AN = hl.agg.group_by(input_mt.hap, AN)
    pre_hap_AF = hl.agg.group_by(input_mt.hap, AF)
    pre_hap_AC_het = hl.agg.group_by(input_mt.hap, AC_het)
    pre_hap_AC_hom = hl.agg.group_by(input_mt.hap, AC_hom)
    pre_hap_AF_hom = hl.agg.group_by(input_mt.hap, AF_hom)
    pre_hap_AF_het = hl.agg.group_by(input_mt.hap, AF_het)
    pre_hap_HL_hist = hl.agg.group_by(input_mt.hap, HL_hist.bin_freq)
    pre_hap_FAF = hl.agg.group_by(
        input_mt.hap,
        hl.experimental.filtering_allele_frequency(hl.int32(AC), hl.int32(AN), 0.95),
    )
    pre_hap_FAF_hom = hl.agg.group_by(
        input_mt.hap,
        hl.experimental.filtering_allele_frequency(
            hl.int32(AC_hom), hl.int32(AN), 0.95
        ),
    )

    # population annotations
    pre_pop_AC = hl.agg.group_by(input_mt.pop, AC)
    pre_pop_AN = hl.agg.group_by(input_mt.pop, AN)
    pre_pop_AF = hl.agg.group_by(input_mt.pop, AF)
    pre_pop_AC_het = hl.agg.group_by(input_mt.pop, AC_het)
    pre_pop_AC_hom = hl.agg.group_by(input_mt.pop, AC_hom)
    pre_pop_AF_hom = hl.agg.group_by(input_mt.pop, AF_hom)
    pre_pop_AF_het = hl.agg.group_by(input_mt.pop, AF_het)
    pre_pop_HL_hist = hl.agg.group_by(input_mt.pop, HL_hist.bin_freq)

    return hl.struct(
        AC=AC,
        AN=AN,
        AF=AF,
        AC_hom=AC_hom,
        AC_het=AC_het,
        hl_hist=HL_hist,
        dp_hist_all=DP_hist_all,
        dp_hist_alt=DP_hist_alt,
        dp_mean=DP_mean,
        mq_mean=MQ_mean,
        tlod_mean=TLOD_mean,
        AF_hom=AF_hom,
        AF_het=AF_het,
        max_hl=max_HL,
        pre_hap_AC=pre_hap_AC,
        pre_hap_AN=pre_hap_AN,
        pre_hap_AF=pre_hap_AF,
        pre_hap_AC_het=pre_hap_AC_het,
        pre_hap_AF_het=pre_hap_AF_het,
        pre_hap_AC_hom=pre_hap_AC_hom,
        pre_hap_AF_hom=pre_hap_AF_hom,
        pre_hap_hl_hist=pre_hap_HL_hist,
        pre_hap_faf=pre_hap_FAF,
        pre_hap_faf_hom=pre_hap_FAF_hom,
        pre_pop_AN=pre_pop_AN,
        pre_pop_AC_het=pre_pop_AC_het,
        pre_pop_AF_het=pre_pop_AF_het,
        pre_pop_AC_hom=pre_pop_AC_hom,
        pre_pop_AF_hom=pre_pop_AF_hom,
        pre_pop_hl_hist=pre_pop_HL_hist,
    )


def standardize_haps(
    input_mt: hl.MatrixTable, annotation: str, haplogroup_order: list
) -> list:
    """
    Convert the dictionary of haplogroup annotations into an array of values in a predefined haplogroup order.

    :param input_mt: MatrixTable
    :param annotation: Annotation to convert and sort
    :param haplogroup_order: Order in which to sort the haplogroups
    :return: Sorted list of haplogroup annotations (the values of the dictionary)
    """
    # Converts haplogroup dictionary to sorted array
    value = [input_mt[annotation][x] for x in haplogroup_order]

    return value


def standardize_pops(
    input_mt: hl.MatrixTable, annotation: str, population_order: list
) -> list:
    """
    Convert the dictionary of population annotations into an array of values in a predefined population order.

    :param input_mt: MatrixTable
    :param annotation: Annotation to convert and sort
    :param population_order: Order in which to sort the populations
    :return: Sorted list of population annotations (the values of the dictionary)
    """
    # Converts haplogroup dictionary to sorted array
    value = [input_mt[annotation][x] for x in population_order]

    return value


def add_quality_histograms(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Add histogram annotations for quality metrics to the MatrixTable.

    :param input_mt: MatrixTable
    :return: MatrixTable annotated with quality metric histograms
    """
    # Generate histogram for site quality metrics across all variants
    # TODO: decide on bin edges
    dp_hist_all_variants = input_mt.aggregate_rows(
        hl.agg.hist(input_mt.dp_mean, 0, 4000, 40)
    )
    input_mt = input_mt.annotate_globals(
        dp_hist_all_variants_bin_freq=dp_hist_all_variants.bin_freq,
        dp_hist_all_variants_n_larger=dp_hist_all_variants.n_larger,
        dp_hist_all_variants_bin_edges=dp_hist_all_variants.bin_edges,
    )

    mq_hist_all_variants = input_mt.aggregate_rows(
        hl.agg.hist(input_mt.mq_mean, 0, 80, 40)
    )  # is 80 the actual max value here?
    input_mt = input_mt.annotate_globals(
        mq_hist_all_variants_bin_freq=mq_hist_all_variants.bin_freq,
        mq_hist_all_variants_n_larger=mq_hist_all_variants.n_larger,
        mq_hist_all_variants_bin_edges=mq_hist_all_variants.bin_edges,
    )

    tlod_hist_all_variants = input_mt.aggregate_rows(
        hl.agg.hist(input_mt.tlod_mean, 0, 40000, 40)
    )
    input_mt = input_mt.annotate_globals(
        tlod_hist_all_variants_bin_freq=tlod_hist_all_variants.bin_freq,
        tlod_hist_all_variants_n_larger=tlod_hist_all_variants.n_larger,
        tlod_hist_all_variants_bin_edges=tlod_hist_all_variants.bin_edges,
    )

    # Generate histogram for overall age distribution
    age_hist_all_samples = input_mt.aggregate_cols(
        hl.agg.hist(input_mt.age, 30, 80, 10)
    )
    input_mt = input_mt.annotate_globals(
        age_hist_all_samples_bin_freq=age_hist_all_samples.bin_freq,
        age_hist_all_samples_n_larger=age_hist_all_samples.n_larger,
        age_hist_all_samples_n_smaller=age_hist_all_samples.n_smaller,
        age_hist_all_samples_bin_edges=age_hist_all_samples.bin_edges,
    )

    # Add age histograms per variant type (heteroplasmic or homoplasmic)
    age_data = age_hists_expr(True, input_mt.GT, input_mt.age)
    input_mt = input_mt.annotate_rows(
        age_hist_hom=age_data.age_hist_hom, age_hist_het=age_data.age_hist_het
    )

    return input_mt


def add_annotations_by_hap_and_pop(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Add variant annotations (such as AC, AN, AF, heteroplasmy histogram, and filtering allele frequency) split by haplogroup and population.

    :param input_mt: MatrixTable
    :return: MatrixTable with variant annotations
    """
    # Order the haplogroup-specific annotations
    list_hap_order = list(set(input_mt.hap.collect()))
    input_mt = input_mt.annotate_globals(hap_order=sorted(list_hap_order))

    # Sanity check for haplogroups (make sure that they at least start with a letter)
    for i in list_hap_order:
        if not re.match("^[A-Z]", i):
            sys.exit(f"Invalid haplogroup {i}, does not start with a letter")

    pre_hap_annotation_labels = [
        "pre_hap_AC",
        "pre_hap_AN",
        "pre_hap_AF",
        "pre_hap_AC_het",
        "pre_hap_AC_hom",
        "pre_hap_AF_hom",
        "pre_hap_AF_het",
        "pre_hap_hl_hist",
        "pre_hap_faf",
        "pre_hap_faf_hom",
    ]

    for i in pre_hap_annotation_labels:
        final_annotation = re.sub(
            "pre_", "", i
        )  # remove "pre" prefix for final annotations
        input_mt = input_mt.annotate_rows(
            **{final_annotation: standardize_haps(input_mt, i, sorted(list_hap_order))}
        )

    # Get a list of indexes where AC of the haplogroup is greater than 0, then get the list of haplogroups with that index
    input_mt = input_mt.annotate_rows(
        alt_haps=hl.enumerate(input_mt.hap_AC)
        .filter(lambda x: x[1] > 0)
        .map(lambda x: input_mt.hap_order[x[0]])
    )
    # Count number of haplogroups containing an alt allele
    input_mt = input_mt.annotate_rows(n_alt_haps=hl.len(input_mt.alt_haps))

    # Calculate hapmax
    input_mt = input_mt.annotate_rows(
        hapmax_AF_hom=input_mt.hap_order[(hl.argmax(input_mt.hap_AF_hom, unique=True))],
        hapmax_AF_het=input_mt.hap_order[(hl.argmax(input_mt.hap_AF_het, unique=True))],
    )

    # Calculate faf hapmax
    input_mt = input_mt.annotate_rows(
        faf_hapmax=hl.max(input_mt.hap_faf), faf_hapmax_hom=hl.max(input_mt.hap_faf_hom)
    )

    # Add populatation annotations
    found_pops = set(input_mt.pop.collect())
    # Order according to POPS
    final_pops = [x for x in POPS if x in found_pops]

    if len(found_pops - set(POPS)) > 0:
        sys.exit("Invalid population found")
    input_mt = input_mt.annotate_globals(pop_order=final_pops)

    pre_pop_annotation_labels = [
        "pre_pop_AN",
        "pre_pop_AC_het",
        "pre_pop_AC_hom",
        "pre_pop_AF_hom",
        "pre_pop_AF_het",
        "pre_pop_hl_hist",
    ]

    for i in pre_pop_annotation_labels:
        # Remove "pre" prefix for final annotations
        final_annotation = re.sub("pre_", "", i)
        input_mt = input_mt.annotate_rows(
            **{final_annotation: standardize_pops(input_mt, i, final_pops)}
        )

    # Drop intermediate annotations
    annotations_to_drop = [
        "pre_hap_AC",
        "pre_hap_AN",
        "pre_hap_AF",
        "pre_hap_AC_het",
        "pre_hap_AC_hom",
        "pre_hap_AF_hom",
        "pre_hap_AF_het",
        "pre_hap_hl_hist",
        "pre_hap_faf",
        "pre_hap_faf_hom",
        "AC_mid_het",
        "AF_mid_het",
        "pre_pop_AN",
        "pre_pop_AC_het",
        "pre_pop_AC_hom",
        "pre_pop_AF_hom",
        "pre_pop_AF_het",
        "pre_pop_hl_hist",
    ]

    input_mt = input_mt.drop(*annotations_to_drop)
    # Last-minute drops (ever add back in?)
    input_mt = input_mt.drop(
        "AC",
        "AF",
        "hap_AC",
        "hap_AF",
        "hap_faf",
        "faf_hapmax",
        "alt_haps",
        "n_alt_haps",
    )

    input_mt = input_mt.annotate_rows(
        filters=hl.if_else(
            input_mt.filters == {"PASS"}, hl.empty_set(hl.tstr), input_mt.filters
        )
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


def get_indel_expr(input_mt: hl.MatrixTable) -> hl.expr.BooleanExpression:
    """
    Generate expression for filtering to indels that should be used to evaluate indel stacks.

    To be considered a variant to be used to evaluate indel stacks, the variant should:
    a) be an indel
    b) have a heteroplasmy level >= 0.01 and <= 0.95
    c) have a PASS genotype

    :param input_mt: MatrixTable
    :return: Expression to be used for determining if a variant is an indel that should to be used to evaluate indel stacks
    """
    indel_expr = (
        hl.is_indel(input_mt.alleles[0], input_mt.alleles[1])
        & (input_mt.HL <= 0.95)
        & (input_mt.HL >= 0.01)
        & (input_mt.FT == {"PASS"})
    )

    return indel_expr


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


def generate_filter_histogram(
    input_mt: hl.MatrixTable, filter_name: str
) -> hl.ArrayExpression:
    """
    Generate histogram for number of indiviudals with the specified sample-level filter at different heteroplasmy levels.

    :param input_mt: MatrixTable
    :param filter_name: Name of sample-filter for which to generate a histogram
    :return: Histogram containing the counts of individuals with a variant filtered by the specified filter name across binned heteroplasmy levels
    """
    filter_histogram = hl.agg.filter(
        hl.str(input_mt.FT).contains(filter_name), hl.agg.hist(input_mt.HL, 0, 1, 10)
    ).bin_freq

    return filter_histogram


def add_filter_annotations(
    input_mt: hl.MatrixTable,
    vaf_filter_threshold: float = 0.01,
    min_het_threshold: float = 0.10,
) -> hl.MatrixTable:
    """
    Generate histogram for number of individuals with the specified sample-level filter at different heteroplasmy levels.

    :param input_mt: MatrixTable
    :param vaf_filter_threshold: Should match vaf_filter_threshold supplied to Mutect2, variants below this value will be set to homoplasmic reference after calculating the common_low_heteroplasmy filter
    :param min_het_threshold: Minimum heteroplasmy level to define a variant as a PASS heteroplasmic variant, genotypes below this threshold will count towards the heteroplasmy_below_min_het_threshold filter and be set to missing
    :return: MatrixTable with added annotations for sample and variant level filters and number of genotypes with heteroplasmy_below_min_het_threshold
    """
    # TODO: pull these from header instead?
    filters = [
        "base_qual",
        "position",
        "strand_bias",
        "weak_evidence",
        "contamination",
        "heteroplasmy_below_min_het_threshold",
    ]

    logger.info("Applying common low heteroplasmy flag...")
    input_mt = apply_common_low_het_flag(input_mt)

    logger.info("Removing low_allele_frac genotypes...")
    input_mt = remove_low_allele_frac_genotypes(input_mt, vaf_filter_threshold)

    logger.info("Applying indel_stack filter...")
    input_mt = apply_indel_stack_filter(input_mt)

    logger.info(
        "Filtering genotypes below with heteroplasmy below the min_het_threshold..."
    )
    input_mt = filter_genotypes_below_min_het_threshold(input_mt, min_het_threshold)
    n_het_below_min_het_threshold = input_mt.aggregate_entries(
        hl.agg.count_where(
            hl.str(input_mt.FT).contains("heteroplasmy_below_min_het_threshold")
        )
    )

    logger.info("Applying npg filter...")
    input_mt = apply_npg_filter(input_mt)

    logger.info("Generating filter histograms and calculating excluded_AC...")
    for i in filters:
        annotation_name = i + "_hist"
        input_mt = input_mt.annotate_rows(
            **{annotation_name: generate_filter_histogram(input_mt, i)}
        )
    input_mt = input_mt.annotate_rows(
        excluded_AC=hl.agg.count_where(input_mt.FT != {"PASS"})
    )

    # Remove "PASS" from filters column if it's not the only filter
    input_mt = input_mt.annotate_rows(
        filters=hl.if_else(
            input_mt.filters != {"PASS"},
            input_mt.filters.remove("PASS"),
            input_mt.filters,
        )
    )

    return input_mt, n_het_below_min_het_threshold


def add_sample_annotations(
    input_mt: hl.MatrixTable, min_hom_threshold: float = 0.95
) -> hl.MatrixTable:
    """
    Add sample annotations to the MatrixTable.

    These sample annotations include the callrate, number of heteroplasmic/homoplasmic SNPs/indels, and the number of singletons.
    Filtered variants (such as artifact_prone_site, npg, and indel_stack) are excluded from these calculations.

    :param input_mt: MatrixTable
    :param min_hom_threshold: Minimum heteroplasmy level to define a variant as homoplasmic
    :return: MatrixTable with sample annotations added
    """
    # Count number of variants
    num_rows = input_mt.count_rows()

    # Add sample qc annotations
    filter_expr = hl.len(input_mt.filters) == 0
    input_mt = input_mt.annotate_cols(
        callrate=hl.agg.filter(
            filter_expr, (hl.agg.count_where(hl.is_defined(input_mt.HL))) / num_rows
        ),
        n_singletons_het=hl.agg.filter(
            filter_expr,
            hl.agg.count_where(
                (input_mt.AC_het == 1)
                & ((input_mt.HL < min_hom_threshold) & (input_mt.HL > 0.0))
            ),
        ),
        n_singletons_hom=hl.agg.filter(
            filter_expr,
            hl.agg.count_where(
                (input_mt.AC_hom == 1) & (input_mt.HL >= min_hom_threshold)
            ),
        ),
        n_snp_het=hl.agg.filter(
            filter_expr,
            hl.agg.count_where(
                (input_mt.HL < min_hom_threshold)
                & (input_mt.HL > 0.0)
                & hl.is_snp(input_mt.alleles[0], input_mt.alleles[1])
            ),
        ),
        n_snp_hom=hl.agg.filter(
            filter_expr,
            hl.agg.count_where(
                (input_mt.HL >= min_hom_threshold)
                & hl.is_snp(input_mt.alleles[0], input_mt.alleles[1])
            ),
        ),
        n_indel_het=hl.agg.filter(
            filter_expr,
            hl.agg.count_where(
                (input_mt.HL < min_hom_threshold)
                & (input_mt.HL > 0.0)
                & (~hl.is_snp(input_mt.alleles[0], input_mt.alleles[1]))
            ),
        ),
        n_indel_hom=hl.agg.filter(
            filter_expr,
            hl.agg.count_where(
                (input_mt.HL >= min_hom_threshold)
                & (~hl.is_snp(input_mt.alleles[0], input_mt.alleles[1]))
            ),
        ),
    )

    return input_mt


def add_vep(input_mt: hl.MatrixTable, run_vep: bool, vep_output: str) -> hl.MatrixTable:
    """
    Add vep annotations to the MatrixTable.

    :param input_mt: MatrixTable
    :param run_vep: Whether or not to run vep
    :param vep_output: Path to the MatrixTable output vep results (either the existing results or where to ouput new vep results)
    :return: MatrixTable with vep annotations
    """
    if run_vep:
        vep_mt = hl.vep(input_mt)
        vep_mt = vep_mt.checkpoint(vep_output, overwrite=True)
    else:
        vep_mt = hl.read_matrix_table(vep_output)

    input_mt = input_mt.annotate_rows(
        vep=vep_mt.index_rows(input_mt.locus, input_mt.alleles).vep
    )
    # TODO: get vep version directly from config file
    input_mt = input_mt.annotate_globals(vep_version="v101")

    # If only filter is END_TRUNC, change lof for LC to HC and remove the END_TRUNC filter
    # Remove SINGLE_EXON flags because all exons are single exon in the mitochondria
    input_mt = input_mt.annotate_rows(
        vep=input_mt.vep.annotate(
            transcript_consequences=input_mt.vep.transcript_consequences.map(
                lambda x: x.annotate(
                    lof=hl.if_else(x.lof_filter == "END_TRUNC", "HC", x.lof),
                    lof_filter=hl.if_else(
                        x.lof_filter == "END_TRUNC", hl.missing(hl.tstr), x.lof_filter
                    ),
                    lof_flags=hl.if_else(
                        x.lof_flags == "SINGLE_EXON", hl.missing(hl.tstr), x.lof_flags
                    ),
                )
            )
        )
    )

    end_trunc_count = input_mt.filter_rows(
        hl.str(input_mt.vep.transcript_consequences[0].lof_filter).contains("END_TRUNC")
    ).count_rows()
    if end_trunc_count > 0:
        sys.exit(
            f"END_TRUNC filter should no longer be present but was found for {end_trunc_count} variants"
        )

    single_exon_count = input_mt.filter_rows(
        hl.str(input_mt.vep.transcript_consequences[0].lof_flags).contains(
            "SINGLE_EXON"
        )
    ).count_rows()
    if single_exon_count > 0:
        sys.exit(
            f"SINGLE_EXON flag should no longer be present but was found for {single_exon_count} variants"
        )

    return input_mt


def add_rsids(input_mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Add rsid annotations to the MatrixTable.

    :param input_mt: MatrixTable
    :return: MatrixTable with rsid annotations added
    """
    dbsnp_import_args = dbsnp.versions["b154"].import_args
    # Replace the contig recoding with just the chrM mapping
    dbsnp_import_args.update({"contig_recoding": {"NC_012920.1": "chrM"}})
    dbsnp_ht = _import_dbsnp(**dbsnp_import_args)

    input_mt = input_mt.annotate_rows(
        rsid=dbsnp_ht[input_mt.locus, input_mt.alleles].rsid
    )
    input_mt = input_mt.annotate_globals(dbsnp_version="b154")

    return input_mt
