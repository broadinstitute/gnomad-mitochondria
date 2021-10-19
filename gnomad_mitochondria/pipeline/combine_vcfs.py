import argparse
import logging
import math
import os

import hail as hl
from hail.utils.java import info
from typing import Dict

META_DICT = {
    "filter": {
        "artifact_prone_site": {
            "Description": "Variant overlaps an artifact-prone site"
        }
    },
    "format": {
        "DP": {"Description": "Depth of coverage", "Number": "1", "Type": "Integer"},
        "FT": {
            "Description": "Sample-level genotype filters",
            "Number": ".",
            "Type": "String",
        },
        "HL": {"Description": "Heteroplasmy level", "Number": "1", "Type": "Float"},
        "MQ": {"Description": "Mapping quality", "Number": "1", "Type": "Float"},
        "TLOD": {
            "Description": "Log 10 likelihood ratio score of variant existing versus not existing",
            "Number": "1",
            "Type": "Float",
        },
    },
}


logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("combine_mitochondria_vcfs_into_mt")
logger.setLevel(logging.INFO)

logger.info("Setting hail flag to avoid array index out of bounds error error...")
# Setting this flag isn't generally recommended, but is needed (since at least Hail version 0.2.75) to avoid an array index out of bounds error until changes are made in future versions of Hail
# TODO: reassess if this flag is still needed for future versions of Hail
hl._set_flags(no_whole_stage_codegen='1')


def collect_vcf_paths(
    participant_data: str, vcf_col_name: str, participants_to_subset: str = None,
) -> Dict[str, str]:
    """
    Create dictionary of VCF paths for only the samples specified in participants_to_subset.

    .. note::
        Participant data should be a tab-delimited file with at minimum columns for:
        - 'entity:participant_id': sample name with prohibited characters replaced with underscores
        - 's': sample name
        - path to the Mutect2 VCF output, where name of this column is supplied to the `vcf_col_name` parameter

    :param participant_data: Participant data (the downloaded data tab from Terra)
    :param vcf_col_name: Name of column that contains VCF output
    :param participants_to_subset: Path to file of participant_ids to which the data should be subset
    :return: Dictionary with sample name as key and path to VCF as value
    """
    vcf_paths = {}
    # Load in data from Terra
    participant_ht = hl.import_table(participant_data)

    # Remove participants that don't have VCF output
    participant_ht.filter(participant_ht[vcf_col_name] != "")

    # Subset participants if specified
    if participants_to_subset:
        participants_of_interest = hl.import_table(
            participants_to_subset
        ).participant.collect()
        participant_ht = participant_ht.filter(
            hl.literal(participants_of_interest).contains(
                participant_ht["entity:participant_id"]
            )
        )

    # Add the vcf path to a dictionary with sample name as key
    df = participant_ht.to_pandas()

    for _, row in df.iterrows():
        vcf_paths[row["s"]] = row[vcf_col_name]

    return vcf_paths


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


def join_mitochondria_vcfs_into_mt(
    vcf_paths: Dict[str, str], output_bucket: str, temp_dir: str, chunk_size: int = 100
) -> hl.MatrixTable:
    """
    Reformat and join individual mitochondrial VCFs into one MatrixTable.

    :param vcf_paths: Dictionary of samples to combine (sample as key, path to VCF as value)
    :param output_bucket: Path to bucket to which results should be written
    :param temp_dir: Path to temporary directory for intermediate results
    :param chunk_size: Number of MatrixTables to join per chunk (the number of individual VCFs that should be combined at a time)
    :return: Joined MatrixTable of samples given in vcf_paths dictionary
    """
    mt_list = []
    for sample, vcf_path in vcf_paths.items():
        try:
            mt = hl.import_vcf(vcf_path, reference_genome="GRCh38")
        except:
            raise ValueError(f"vcf path {vcf_path} does not exist for sample {sample}")

        # Because the vcfs are split, there is only one AF value, although misinterpreted as an array because Number=A in VCF header
        # Second value of MMQ is the value of the mapping quality for the alternate allele
        # Add FT annotation for sample genotype filters (pull these from filters annotations of the single-sample VCFs)
        mt = mt.select_entries("DP", HL=mt.AF[0])
        mt = mt.annotate_entries(
            MQ=hl.float(mt.info["MMQ"][1]),
            TLOD=mt.info["TLOD"][0],
            FT=hl.if_else(hl.len(mt.filters) == 0, {"PASS"}, mt.filters),
        )
        # Use GRCh37 as reference as most external resources added in downstream scripts use GRCh37 contig names (although note that the actual sequences of the mitochondira in both GRCh37 and GRCh38 are the same)
        mt = mt.key_rows_by(
            locus=hl.locus("MT", mt.locus.position, reference_genome="GRCh37"),
            alleles=mt.alleles,
        )
        mt = mt.key_cols_by(s=sample)
        mt = mt.select_rows()
        mt_list.append(mt)

    combined_mt = multi_way_union_mts(mt_list, temp_dir, chunk_size)

    return combined_mt


def remove_genotype_filters(
    mt: hl.MatrixTable,
    filters_to_remove: set = {
        "possible_numt",
        "mt_many_low_hets",
        "FAIL",
        "blacklisted_site",
    },
) -> hl.MatrixTable:
    """
    Remove unneeded sample-level genotype filters (in FT field of the VCF) specified by the filters_to_remove parameter.

    By default, remove the 'possible_numt', 'mt_many_low_hets', and 'FAIL' filters because these filters were found to have low performance.
    Also remove the 'blacklisted_site' filter because this filter did not always behave as expected in early GATK versions. This filter can be reimplemented with the apply_mito_artifact_filter function.

    :param mt: MatrixTable containing genotype filters in the FT field of the VCF that should be removed
    :param filters_to_remove: List of genptype filters (in FT field of VCF) that should be removed from the entries
    :return: MatrixTable with specific genotype filters (in FT field of VCF) removed
    """
    mt = mt.annotate_entries(FT=mt.FT.difference(filters_to_remove))

    # If no filters exist after removing those specified above, set the FT field to PASS
    mt = mt.annotate_entries(FT=hl.if_else(hl.len(mt.FT) == 0, {"PASS"}, mt.FT))

    return mt


def determine_hom_refs(
    mt: hl.MatrixTable, coverage_mt_path: str, minimum_homref_coverage: int = 100
) -> hl.MatrixTable:
    """
    Use coverage to distinguish between homref and missing sites.

    :param mt: MatrixTable from initial multi-sample merging, without homref sites determined
    :param coverage_mt_path: MatrixTable of sample level coverage at each position (per-sample and per-base; can be generated by running annotate_coverage.py)
    :param minimum_homref_coverage: Minimum depth of coverage required to call a genotype homoplasmic reference rather than missing
    :return: MatrixTable with missing genotypes converted to homref depending on coverage
    """
    # Convert coverage to build GRCh37 to match contig names
    # Note: the mitochondrial reference genome is the same for GRCh38 and GRCh37
    coverages = hl.read_matrix_table(coverage_mt_path)
    coverages = coverages.key_rows_by(
        locus=hl.locus("MT", coverages.locus.position, reference_genome="GRCh37")
    )

    mt = mt.annotate_entries(
        DP=hl.if_else(hl.is_missing(mt.HL), coverages[mt.locus, mt.s].coverage, mt.DP)
    )

    hom_ref_expr = hl.is_missing(mt.HL) & (mt.DP > minimum_homref_coverage)

    mt = mt.annotate_entries(
        HL=hl.if_else(hom_ref_expr, 0.0, mt.HL),
        FT=hl.if_else(hom_ref_expr, {"PASS"}, mt.FT),
        DP=hl.if_else(
            hl.is_missing(mt.HL) & (mt.DP <= minimum_homref_coverage),
            hl.null(hl.tint32),
            mt.DP,
        ),
    )

    return mt


def apply_mito_artifact_filter(
    mt: hl.MatrixTable, artifact_prone_sites_path: str,
) -> hl.MatrixTable:
    """
    Add in artifact_prone_site filter.

    :param mt: MatrixTable to be annotated with artifact_prone_sites filter
    :param artifact_prone_sites_path: Path to BED file of artifact_prone_sites to flag in the filters column
    :return: MatrixTable with artifact_prone_sites filter
    """
    # Apply "artifact_prone_site" filter to any SNP or deletion that spans a known problematic site
    bed = hl.import_bed(artifact_prone_sites_path)
    bed = bed.annotate(target="artifact")

    # Create a region annotation containing the interval that the variant overlaps (for SNP will be one position, but will be longer for deletions based on the length of the deletion)
    mt = mt.annotate_rows(
        region=hl.interval(
            hl.locus("MT", mt.locus.position, reference_genome="GRCh37"),
            hl.locus(
                "MT",
                mt.locus.position + hl.len(mt.alleles[0]) - 1,
                reference_genome="GRCh37",
            ), includes_end = True
        )
    )

    # Annotate if the start of the variant overlaps an interval in the bed file
    mt = mt.annotate_rows(start_overlaps=bed.index(mt.region.start, all_matches=True))

    # Annotate if the end of the variant overlaps an interval in the bed file
    mt = mt.annotate_rows(end_overlaps=bed.index(mt.region.end, all_matches=True))

    # Create struct containing locus and allele (need to the check if any position of the allele overlaps an artifact-prone site, not just the locus)
    mt_temp = mt.annotate_rows(variant=hl.struct(locus=mt.locus, alleles=mt.alleles))
    mt_temp = mt_temp.key_rows_by(mt_temp.region)

    # Need to account for cases where the start and end of the variant interval don't fall within a bed interval, but start before and after the interval (the bed interval falls completely within the variant interval)
    bed_temp = bed.annotate(
        contained_mt_alleles=mt_temp.index_rows(
            bed.interval.start, all_matches=True
        ).variant
    )

    # Explode so that each allele is on its own row and create locus and allele annotations
    bed_temp = bed_temp.explode(bed_temp.contained_mt_alleles).rename(
        {"contained_mt_alleles": "contained_mt_allele"}
    )
    bed_temp = bed_temp.annotate(
        locus=bed_temp.contained_mt_allele.locus,
        alleles=bed_temp.contained_mt_allele.alleles,
    )
    bed_temp = bed_temp.key_by(bed_temp.locus, bed_temp.alleles)

    # Annotate back onto the original mt cases where the bed interval falls completely within the variant interval
    mt = mt.annotate_rows(start_and_end_span=bed_temp[mt.locus, mt.alleles].target)

    # Add artifact-prone site filter to any SNP/deletion that starts within, ends within, or completely overlaps an artifact-prone site
    mt = mt.annotate_rows(
        filters=hl.if_else(
            (hl.len(mt.start_overlaps) > 0)
            | (hl.len(mt.end_overlaps) > 0)
            | (hl.is_defined(mt.start_and_end_span)),
            {"artifact_prone_site"},
            {"PASS"},
        )
    )

    mt = mt.drop("region", "start_overlaps", "end_overlaps", "start_and_end_span")

    return mt


def main(args):
    participant_data = args.participant_data
    coverage_mt_path = args.coverage_mt_path
    vcf_col_name = args.vcf_col_name
    artifact_prone_sites_path = args.artifact_prone_sites_path
    output_bucket = args.output_bucket
    file_name = args.file_name
    minimum_homref_coverage = args.minimum_homref_coverage
    chunk_size = args.chunk_size

    output_path_mt = f"{output_bucket}/raw_combined.mt"

    if args.overwrite == False and hl.hadoop_exists(output_path_mt):
        logger.warn(
            "Overwrite is set to False but file already exists at %s, script will run but output will not be written",
            output_path_mt,
        )

    logger.info("Collecting VCF paths for samples to subset...")
    vcf_paths = collect_vcf_paths(
        participant_data, vcf_col_name, args.participants_to_subset
    )

    logger.info("Combining VCFs...")
    combined_mt = join_mitochondria_vcfs_into_mt(
        vcf_paths, output_bucket, args.temp_dir, chunk_size
    )
    combined_mt = combined_mt.checkpoint(output_path_mt, overwrite=args.overwrite)

    logger.info("Removing select sample-level filters...")
    combined_mt = remove_genotype_filters(combined_mt)

    logger.info("Determining homoplasmic reference sites...")
    combined_mt = determine_hom_refs(
        combined_mt, coverage_mt_path, minimum_homref_coverage
    )

    logger.info("Applying artifact_prone_site fiter...")
    combined_mt = apply_mito_artifact_filter(combined_mt, artifact_prone_sites_path)

    logger.info("Writing combined MT and VCF...")
    # Set the file names for output files
    out_vcf = f"{output_bucket}/{file_name}.vcf.bgz"
    out_mt = f"{output_bucket}/{file_name}.mt"

    combined_mt = combined_mt.checkpoint(out_mt, overwrite=args.overwrite)
    # For the VCF output, join FT values by semicolon
    combined_mt = combined_mt.annotate_entries(
        FT=hl.str(";").join(hl.array(combined_mt.FT))
    )
    hl.export_vcf(combined_mt, out_vcf, metadata=META_DICT)


if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="This script combines individual mitochondria VCF files into one MatrixTable, determines homoplasmic reference sites, and applies an artifact_prone_site filter"
    )
    p.add_argument(
        "-p",
        "--participant-data",
        help="Participant data (the downloaded data tab from Terra), should be a tab-delimited file with at minimum columns for 'entity:participant_id' (sample name with prohibited characters replaced with underscores), 's' (sample name), and VCF output (path to the Mutect2 VCF output, where the name of this column is supplied to the `vcf_col_name` parameter)",
        required=True,
    )
    p.add_argument(
        "-c",
        "--coverage-mt-path",
        help="Path to MatrixTable of sample-level coverage (per-sample and per-base, can be generated by running annotate_coverage.py)",
        required=True,
    )
    p.add_argument(
        "-v",
        "--vcf-col-name",
        help="Name of column in participant data file that contains the path to the VCF output by Mutect2",
        required=True,
    )
    p.add_argument(
        "-a",
        "--artifact-prone-sites-path",
        help="Path to BED file of artifact-prone sites to flag in the FILTER column",
        required=True,
    )
    p.add_argument(
        "-o",
        "--output-bucket",
        help="Path to bucket to which results should be written",
        required=True,
    )
    p.add_argument(
        "-t",
        "--temp-dir",
        help="Temporary directory to use for intermediate outputs",
        required=True,
    )
    p.add_argument(
        "-f",
        "--file-name",
        help="File name to use for output files (will be used the the .vcf.bgz and .mt outputs)",
        required=True,
    )
    p.add_argument(
        "-s",
        "--participants-to-subset",
        help="Path to txt file of participant_ids to which the data should be subset (file should contain header (named 'participant') and one line for each participant_id matching the 'entity:participant_id's supplied in Terra",
    )
    p.add_argument(
        "--minimum-homref-coverage",
        help="Minimum depth of coverage required to call a genotype homoplasmic reference rather than missing",
        type=int,
        default=100,
    )
    p.add_argument(
        "--chunk-size",
        help="Chunk size to use for combining VCFs (the number of individual VCFs that should be combined at a time)",
        type=int,
        default=100,
    )
    p.add_argument("--overwrite", help="Overwrites existing files", action="store_true")

    args = p.parse_args()

    main(args)
