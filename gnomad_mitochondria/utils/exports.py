import hail as hl

from collections import Counter
from textwrap import dedent

from gnomad.utils.vep import vep_struct_to_csq


def generate_output_paths(
    output_dir: str, file_name: str, subset_name: str, extension: str
) -> list:
    """
    Generate output paths for results files based on the given output directory, file name, subset name, and extension.

    :param output_dir: Output directory to which results should be output
    :param file_name: Name of the file, preceeds subset_name
    :param subset_name: Name that should be appended to output file names
    :param extension: Extension for the output file
    :return: Path for the file
    """
    # set up output paths for callset
    file_path = f"{output_dir}/{file_name}{subset_name}.{extension}"

    return file_path


def report_stats(
    input_mt: hl.MatrixTable,
    output_dir: str,
    pass_only: bool,
    n_samples_below_cn: int,
    n_samples_above_cn: int,
    n_samples_contam: int,
    n_het_below_min_het_threshold: int,
    min_het_threshold: float = 0.10,
    min_hom_threshold: float = 0.95,
) -> None:
    """
    Generate output report with basic stats.

    :param input_mt: MatrixTable
    :param output_dir: Output directory to which results should be output
    :param pass_only: Whether or not directory should be filtered to pass_only variants
    :param n_samples_below_cn: Number of samples removed because mitochondrial number is less than 50
    :param n_samples_above_cn: Number of samples removed because mitochondrial number is above 500
    :param n_samples_contam: Number of samples removed because of contamination
    :param n_het_below_min_het_threshold: Number of genotypes filtered because the heteroplasmy levels was below the min_het_threshold
    :param min_het_threshold: Minimum heteroplasmy level to define a variant as a PASS heteroplasmic variant, genotypes below this threshold will count towards the heteroplasmy_below_min_het_threshold filter and be set to missing
    :param min_hom_threshold: Minimum heteroplasmy level to define a variant as homoplasmic
    :return: None
    """
    if pass_only:
        suffix = "_pass"
        input_mt = input_mt.filter_rows(hl.len(input_mt.filters) == 0)
    else:
        suffix = ""
    out_stats = hl.hadoop_open(f"{output_dir}/stats{suffix}.txt", "w")

    if pass_only:
        out_stats.write("Below metrics are for PASS-only variants\n\n")

    # Report numbers of filtered samples/genotypes
    out_stats.write(
        f"Number of samples removed because contamination above 2%: {n_samples_contam}\n"
    )
    out_stats.write(
        f"Number of samples removed because mitochondrial copy number below 50: {n_samples_below_cn}\n"
    )
    out_stats.write(
        f"Number of samples removed because mitochondrial copy number above 500: {n_samples_above_cn}\n"
    )
    out_stats.write(
        f'Number of genotypes filtered because "heteroplasmy_below_min_het_threshold": {n_het_below_min_het_threshold}\n\n'
    )

    # Count variant, samples, bases
    unique_variants, samples = input_mt.count()
    out_stats.write(f"Number of samples: {samples}\n")
    out_stats.write(f"Number of unique variants: {unique_variants}\n")
    bases_w_variant = len(set(input_mt.locus.position.collect()))
    out_stats.write(f"Number of bases with variation: {bases_w_variant}\n\n")

    # Count number of filters
    for filter_name, filter_count in Counter(
        [i for sublist in input_mt.filters.collect() for i in sublist]
    ).items():
        out_stats.write(
            f'Number of variants with "{filter_name}" filter: {filter_count} variants\n'
        )

    # Calculate row stats
    row_stats = input_mt.aggregate_rows(
        hl.struct(
            common_low_het_count=hl.agg.count_where(input_mt.common_low_heteroplasmy),
            het_only_sites=hl.agg.count_where(
                (input_mt.AC_het > 0) & (input_mt.AC_hom == 0)
            ),
            hom_only_sites=hl.agg.count_where(
                (input_mt.AC_hom > 0) & (input_mt.AC_het == 0)
            ),
            het_and_hom_sites=hl.agg.count_where(
                (input_mt.AC_hom > 0) & (input_mt.AC_het > 0)
            ),
            hap_defining_sites=hl.agg.count_where(input_mt.hap_defining_variant),
            snps=hl.agg.count_where(
                hl.is_snp(input_mt.alleles[0], input_mt.alleles[1])
            ),
            indels=hl.agg.count_where(
                hl.is_indel(input_mt.alleles[0], input_mt.alleles[1])
            ),
            transitions=hl.agg.count_where(
                hl.is_transition(input_mt.alleles[0], input_mt.alleles[1])
            ),
            transversions=hl.agg.count_where(
                hl.is_transversion(input_mt.alleles[0], input_mt.alleles[1])
            ),
        )
    )

    # Calculate col stats
    col_stats = input_mt.aggregate_cols(
        hl.struct(
            unique_haplogroups=hl.len(hl.agg.collect_as_set(input_mt.major_haplogroup)),
            unique_top_level_haplogroups=hl.len(hl.agg.collect_as_set(input_mt.hap)),
        )
    )

    # Calculate entry stats
    entry_stats = input_mt.aggregate_entries(
        hl.struct(
            total_variants=hl.agg.count_where(input_mt.HL > 0),
            total_hom_variants=hl.agg.count_where(input_mt.HL >= min_hom_threshold),
            total_het_variants=hl.agg.count_where(
                (input_mt.HL < min_hom_threshold) & (input_mt.HL >= min_het_threshold)
            ),
            min_hl=hl.agg.filter(input_mt.HL > 0, hl.agg.min(input_mt.HL)),
            max_hl=hl.agg.filter(input_mt.HL > 0, hl.agg.max(input_mt.HL)),
        )
    )

    # Count number of flags
    out_stats.write(
        f'Number of variants with "common_low_heteroplasmy" flag: {row_stats["common_low_het_count"]} variants\n\n'
    )

    # Count variants
    out_stats.write(f'Total number of variants: {entry_stats["total_variants"]}\n')

    # Count of homoplasmic/heteroplasmic variants
    out_stats.write(
        f'Number of homoplasmic-only sites: {row_stats["hom_only_sites"]}\n'
    )
    out_stats.write(
        f'Number of heteroplasmic-only sites: {row_stats["het_only_sites"]}\n'
    )
    out_stats.write(f'Number of het and hom sites: {row_stats["het_and_hom_sites"]}\n')

    percent_hom = round(
        entry_stats["total_hom_variants"] / entry_stats["total_variants"], 2
    )
    percent_het = round(
        entry_stats["total_het_variants"] / entry_stats["total_variants"], 2
    )
    out_stats.write(
        f'Total number of homoplasmic variants: {entry_stats["total_hom_variants"]}\n'
    )
    out_stats.write(f"Percent homoplasmic variants: {percent_hom}\n")
    out_stats.write(
        f'Total number of heteroplasmic variants: {entry_stats["total_het_variants"]}\n'
    )
    out_stats.write(f"Percent heteroplasmic variants: {percent_het}\n\n")

    out_stats.write(f'Minimum heteroplasmy detected: {entry_stats["min_hl"]}\n')
    out_stats.write(f'Maximum heteroplasmy detected: {entry_stats["max_hl"]}\n\n')

    # Count number of snps and indels
    out_stats.write(f'Number of SNPs: {row_stats["snps"]}\n')
    out_stats.write(f'Number of indels: {row_stats["indels"]}\n')

    # Count number of transitions and transversions
    out_stats.write(f'Number of transitions: {row_stats["transitions"]}\n')
    out_stats.write(f'Number of transversions: {row_stats["transversions"]}\n')
    out_stats.write(
        f'Number of haplogroup defining variants: {row_stats["hap_defining_sites"]}\n\n'
    )

    # Count number of haplogroups
    out_stats.write(
        f'Number of unique haplogroups: {col_stats["unique_haplogroups"]}\n'
    )
    out_stats.write(
        f'Number of top-level haplogroups: {col_stats["unique_top_level_haplogroups"]}\n'
    )

    out_stats.close()


def export_simplified_variants(input_ht: hl.Table, output_dir: str) -> None:
    """
    Export a text file containing only several high-level variant annotations.

    :param input_ht: Hail Table of variants
    :param output_dir: Output directory to which results should be output
    :return: None
    """
    reduced_ht = (
        input_ht.key_by(
            chromosome=input_ht.locus.contig,
            position=input_ht.locus.position,
            ref=input_ht.alleles[0],
            alt=input_ht.alleles[1],
        )
        .select("filters", "AC_hom", "AC_het", "AF_hom", "AF_het", "AN", "max_hl")
        .rename({"max_hl": "max_observed_heteroplasmy"})
    )
    reduced_ht = reduced_ht.annotate(
        filters=hl.if_else(
            hl.len(reduced_ht.filters) == 0,
            "PASS",
            hl.str(",").join(hl.array(reduced_ht.filters)),
        )
    )

    reduced_ht.export(f"{output_dir}/reduced_annotations.txt")


def change_to_grch38_chrm(input_mt: hl.MatrixTable) -> None:
    """
    Change build to GRCh38 and filters reference genome to chrM.

    :param input_mt: MatrixTable
    :return: MatrixTable with GRCh38 reference genome subsetted to just chrM (so that extraneous contigs will be excluded from VCF output)
    """
    ref = hl.get_reference("GRCh38")
    my_ref = hl.ReferenceGenome(
        "GRCh38_chrM", contigs=["chrM"], lengths={"chrM": ref.lengths["chrM"]}
    )
    assert "chrM" in ref.contigs
    input_mt = input_mt.key_rows_by(
        locus=hl.locus("chrM", input_mt.locus.position, reference_genome="GRCh38_chrM"),
        alleles=input_mt.alleles,
    )

    return input_mt


def format_vcf(
    input_mt: hl.MatrixTable,
    output_dir: str,
    min_hom_threshold: float = 0.95,
    vaf_filter_threshold: float = 0.01,
    min_het_threshold: float = 0.10,
) -> dict:
    """
    Generate dictionary for VCF header annotations.

    :param input_mt: MatrixTable
    :param min_hom_threshold: Minimum heteroplasmy level to define a variant as homoplasmic
    :param vaf_filter_threshold: Should match vaf_filter_threshold supplied to Mutect2, variants below this value will be set to homoplasmic reference after calculating the common_low_heteroplasmy filter
    :param min_het_threshold: Minimum heteroplasmy level to define a variant as a PASS heteroplasmic variant, genotypes below this threshold will count towards the heteroplasmy_below_min_het_threshold filter and be set to missing
    :param output_dir: Output directory to which appended header info should be written
    :return: MatrixTable with VCF annotations in the info field and dictionary of filter, info, and format fields to be output in the VCF header; path of VCF headers to append
    """
    input_mt = change_to_grch38_chrm(input_mt)

    haplogroup_order = hl.eval(input_mt.hap_order)
    population_order = hl.eval(input_mt.pop_order)

    age_hist_hom_bin_edges = input_mt.age_hist_hom.bin_edges.take(1)[0]
    age_hist_het_bin_edges = input_mt.age_hist_het.bin_edges.take(1)[0]
    hl_hist_bin_edges = input_mt.hl_hist.bin_edges.take(1)[0]
    dp_hist_all_bin_edges = input_mt.dp_hist_all.bin_edges.take(1)[0]
    dp_hist_alt_bin_edges = input_mt.dp_hist_alt.bin_edges.take(1)[0]

    input_mt = input_mt.annotate_rows(
        hl_hist=input_mt.hl_hist.bin_freq,
        age_hist_hom_bin_freq=input_mt.age_hist_hom.bin_freq,
        age_hist_hom_n_smaller=input_mt.age_hist_hom.n_smaller,
        age_hist_hom_n_larger=input_mt.age_hist_hom.n_larger,
        age_hist_het_bin_freq=input_mt.age_hist_het.bin_freq,
        age_hist_het_n_smaller=input_mt.age_hist_het.n_smaller,
        age_hist_het_n_larger=input_mt.age_hist_het.n_larger,
        dp_hist_all_n_larger=input_mt.dp_hist_all.n_larger,
        dp_hist_alt_n_larger=input_mt.dp_hist_alt.n_larger,
        dp_hist_all_bin_freq=input_mt.dp_hist_all.bin_freq,
        dp_hist_alt_bin_freq=input_mt.dp_hist_alt.bin_freq,
    )

    input_mt = input_mt.annotate_rows(vep=vep_struct_to_csq(input_mt.vep))

    # Get length of annotations to use in Number fields in the VCF where necessary
    len_hap_hl_hist = len(input_mt.hap_hl_hist.take(1)[0])
    len_pop_hl_hist = len(input_mt.pop_hl_hist.take(1)[0])

    # Output appended header info to file
    vcf_header_file = output_dir + "/extra_fields_for_header.tsv"
    appended_vcf_header = dedent(
        f"""
    ##VEP version={hl.eval(input_mt.vep_version)}
    ##dbSNP version={hl.eval(input_mt.dbsnp_version)}
    ##age distributions=bin_edges:{hl.eval(input_mt.age_hist_all_samples_bin_edges)}, bin_freq:{hl.eval(input_mt.age_hist_all_samples_bin_freq)}, n_smaller:{hl.eval(input_mt.age_hist_all_samples_n_smaller)}, n_larger:{hl.eval(input_mt.age_hist_all_samples_n_larger)}

    """
    )
    with hl.hadoop_open(vcf_header_file, "w") as out:
        out.write(appended_vcf_header)

    # Drop intermediate annotations
    input_mt = input_mt.drop(
        "region",
        "variant_context",
        "age_hist_het",
        "age_hist_hom",
        "dp_hist_all",
        "dp_hist_alt",
    )

    ht = input_mt.rows()
    # Move row annotations into info struct
    input_mt = input_mt.annotate_rows(info=hl.struct())
    input_mt = input_mt.annotate_rows(
        info=input_mt.info.annotate(**ht[input_mt.row_key]).drop("rsid")
    ).select_rows(
        "rsid", "filters", "info"
    )  # create info annotation

    # Convert "rsid" array to str for VCF output
    input_mt = input_mt.annotate_rows(rsid=hl.str(";").join(input_mt.rsid))

    # Convert "," to "|" for array annotations
    for key, value in input_mt.row_value.info.items():
        if str(value).startswith("<Array"):
            if str(value.dtype).startswith("array<array"):
                # If value is an array of arrays, only replace the commas within each individual array
                input_mt = input_mt.annotate_rows(
                    info=input_mt.info.annotate(
                        **{
                            key: hl.map(
                                lambda x: hl.delimit(x, delimiter="|"),
                                input_mt.info[key],
                            )
                        }
                    )
                )
            else:
                input_mt = input_mt.annotate_rows(
                    info=input_mt.info.annotate(
                        **{key: hl.delimit(input_mt.info[key], delimiter="|")}
                    )
                )

    meta_dict = {
        "filter": {
            "artifact_prone_site": {
                "Description": "Variant overlaps site that is commonly reported in literature to be artifact prone"
            },
            "npg": {
                "Description": "No-pass-genotypes site (no individuals were PASS for the variant)"
            },
            "indel_stack": {
                "Description": "Allele where all samples with the variant call had at least 2 different heteroplasmic indels called at the position"
            },
        },
        "info": {
            "variant_collapsed": {
                "Description": "Variant in format of RefPosAlt",
                "Number": "1",
                "Type": "String",
            },
            "hap_defining_variant": {
                "Description": "Present if variant is present as a haplogroup defining variant in PhyloTree build 17",
                "Number": "0",
                "Type": "Flag",
            },
            "common_low_heteroplasmy": {
                "Description": f"Present if variant is found at an overall frequency of .001 across all samples with a heteroplasmy level > 0 and < 0.50 (includes variants <{vaf_filter_threshold} heteroplasmy which are subsequently filtered)",
                "Number": "0",
                "Type": "Flag",
            },
            "AN": {
                "Description": "Overall allele number (number of samples with non-missing genotype)",
                "Number": "1",
                "Type": "Integer",
            },
            "AC_hom": {
                "Description": f"Allele count restricted to variants with a heteroplasmy level >= {min_hom_threshold}",
                "Number": "1",
                "Type": "Integer",
            },  # should put in threshold variable
            "AC_het": {
                "Description": f"Allele count restricted to variants with a heteroplasmy level >= 0.10 and < {min_hom_threshold}",
                "Number": "1",
                "Type": "Integer",
            },
            "hl_hist": {
                "Description": f"Histogram of heteroplasmy levels; bin edges are: {hl_hist_bin_edges}",
                "Number": "1",
                "Type": "String",
            },
            "hap_hl_hist": {
                "Description": f"Histogram of heteroplasmy levels for each haplogroup; bin edges are: {hl_hist_bin_edges}, haplogroup order: {haplogroup_order}",
                "Number": f"{len_hap_hl_hist}",
                "Type": "String",
            },
            "AF_hom": {
                "Description": f"Allele frequency restricted to variants with a heteroplasmy level >= {min_hom_threshold}",
                "Number": "1",
                "Type": "Float",
            },
            "AF_het": {
                "Description": f"Allele frequency restricted to variants with a heteroplasmy level >= 0.10 and < {min_hom_threshold}",
                "Number": "1",
                "Type": "Float",
            },
            "max_hl": {
                "Description": "Maximum heteroplasmy level observed among all samples for that variant",
                "Number": "1",
                "Type": "Float",
            },
            "hap_AN": {
                "Description": f"List of overall allele number for each haplogroup, haplogroup order: {haplogroup_order}",
                "Number": "1",
                "Type": "String",
            },
            "hap_AC_het": {
                "Description": f"List of AC_het for each haplogroup, haplogroup order: {haplogroup_order}",
                "Number": "1",
                "Type": "String",
            },
            "hap_AC_hom": {
                "Description": f"List of AC_hom for each haplogroup, haplogroup order: {haplogroup_order}",
                "Number": "1",
                "Type": "String",
            },
            "hap_AF_hom": {
                "Description": f"List of AF_hom for each haplogroup, haplogroup order: {haplogroup_order}",
                "Number": "1",
                "Type": "String",
            },
            "hap_AF_het": {
                "Description": f"List of AF_het for each haplogroup, haplogroup order: {haplogroup_order}",
                "Number": "1",
                "Type": "String",
            },
            "hap_faf_hom": {
                "Description": f"List of filtering allele frequency for each haplogroup restricted to homoplasmic variants, haplogroup order: {haplogroup_order}",
                "Number": "1",
                "Type": "String",
            },
            "hapmax_AF_hom": {
                "Description": "Haplogroup with maximum AF_hom",
                "Number": "1",
                "Type": "String",
            },
            "hapmax_AF_het": {
                "Description": "Haplogroup with maximum AF_het",
                "Number": "1",
                "Type": "String",
            },
            "faf_hapmax_hom": {
                "Description": "Maximum filtering allele frequency across haplogroups restricted to homoplasmic variants",
                "Number": "1",
                "Type": "Float",
            },
            "bin_edges_hl_hist": {
                "Description": "Bin edges for histogram of heteroplasmy levels",
                "Number": "1",
                "Type": "String",
            },
            "pop_AN": {
                "Description": f"List of overall allele number for each population, population order: {population_order}",
                "Number": "1",
                "Type": "String",
            },
            "pop_AC_het": {
                "Description": f"List of AC_het for each population, population order: {population_order}",
                "Number": "1",
                "Type": "String",
            },
            "pop_AC_hom": {
                "Description": f"List of AC_hom for each population, population order: {population_order}",
                "Number": "1",
                "Type": "String",
            },
            "pop_AF_hom": {
                "Description": f"List of AF_hom for each population, population order: {population_order}",
                "Number": "1",
                "Type": "String",
            },
            "pop_AF_het": {
                "Description": f"List of AF_het for each population, population order: {population_order}",
                "Number": "1",
                "Type": "String",
            },
            "pop_hl_hist": {
                "Description": f"Histogram of heteroplasmy levels for each population; bin edges are: {hl_hist_bin_edges}, population order: {population_order}",
                "Number": f"{len_pop_hl_hist}",
                "Type": "String",
            },
            "age_hist_hom_bin_freq": {
                "Description": f"Histogram of ages of individuals with a homoplasmic variant; bin edges are: {age_hist_hom_bin_edges}",
                "Number": "1",
                "Type": "String",
            },
            "age_hist_hom_n_smaller": {
                "Description": "Count of age values falling below lowest histogram bin edge for individuals with a homoplasmic variant",
                "Number": "1",
                "Type": "String",
            },
            "age_hist_hom_n_larger": {
                "Description": "Count of age values falling above highest histogram bin edge for individuals with a homoplasmic variant",
                "Number": "1",
                "Type": "String",
            },
            "age_hist_het_bin_freq": {
                "Description": f"Histogram of ages of individuals with a heteroplasmic variant; bin edges are: {age_hist_het_bin_edges}",
                "Number": "1",
                "Type": "String",
            },
            "age_hist_het_n_smaller": {
                "Description": "Count of age values falling below lowest histogram bin edge for individuals with a heteroplasmic variant",
                "Number": "1",
                "Type": "String",
            },
            "age_hist_het_n_larger": {
                "Description": "Count of age values falling above highest histogram bin edge for individuals with a heteroplasmic variant",
                "Number": "1",
                "Type": "String",
            },
            "dp_hist_all_n_larger": {
                "Description": "Count of dp values falling above highest histogram bin edge for all individuals",
                "Number": "1",
                "Type": "String",
            },
            "dp_hist_alt_n_larger": {
                "Description": "Count of dp values falling above highest histogram bin edge for individuals with the alternative allele",
                "Number": "1",
                "Type": "String",
            },
            "dp_hist_all_bin_freq": {
                "Description": f"Histogram of dp values for all individuals; bin edges are: {dp_hist_all_bin_edges}",
                "Number": "1",
                "Type": "String",
            },
            "dp_hist_alt_bin_freq": {
                "Description": f"Histogram of dp values for individuals with the alternative allele; bin edges are: {dp_hist_alt_bin_edges}",
                "Number": "1",
                "Type": "String",
            },
            "dp_mean": {
                "Description": "Mean depth across all individuals for the site",
                "Number": "1",
                "Type": "Float",
            },
            "mq_mean": {
                "Description": "Mean MMQ (median mapping quality) across individuals with a variant for the site",
                "Number": "1",
                "Type": "Float",
            },
            "tlod_mean": {
                "Description": "Mean TLOD (Log 10 likelihood ratio score of variant existing versus not existing) across individuals with a variant for the site",
                "Number": "1",
                "Type": "Float",
            },
            "pon_mt_trna_prediction": {
                "Description": "tRNA pathogenicity classification from PON-mt-tRNA",
                "Number": "1",
                "Type": "String",
            },
            "pon_ml_probability_of_pathogenicity": {
                "Description": "tRNA ML_probability_of_pathogenicity from PON-mt-tRNA",
                "Number": "1",
                "Type": "Float",
            },
            "mitotip_score": {
                "Description": "MitoTip raw score",
                "Number": "1",
                "Type": "Float",
            },
            "mitotip_trna_prediction": {
                "Description": "MitoTip score interpretation",
                "Number": "1",
                "Type": "String",
            },
            "vep": {
                "Description": "Consequence annotations from Ensembl VEP; note that the SINGLE_EXON flag and END_TRUNC filters have been removed from the LOFTEE annotations to avoid misinterpretation in context of the mitochondrial genome. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|ALLELE_NUM|DISTANCE|STRAND|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info",
                "Number": ".",
                "Type": "String",
            },
            "filters": {
                "Description": "Site-level filters",
                "Number": ".",
                "Type": "String",
            },
            "base_qual_hist": {
                "Description": f"Histogram of number of individuals failing the base_qual filter (alternate allele median base quality) across heteroplasmy levels, bin edges are: {hl_hist_bin_edges}",
                "Number": "1",
                "Type": "String",
            },
            "heteroplasmy_below_min_het_threshold_hist": {
                "Description": f"Histogram of number of individuals with a heteroplasmy level below {min_het_threshold}, bin edges are: {hl_hist_bin_edges}",
                "Number": "1",
                "Type": "String",
            },
            "position_hist": {
                "Description": f"Histogram of number of individuals failing the position filter (median distance of alternate variants from end of reads) across heteroplasmy levels, bin edges are: {hl_hist_bin_edges}",
                "Number": "1",
                "Type": "String",
            },
            "strand_bias_hist": {
                "Description": f"Histogram of number of individuals failing the strand_bias filter (evidence for alternate allele comes from one read direction only) across heteroplasmy levels, bin edges are: {hl_hist_bin_edges}",
                "Number": "1",
                "Type": "String",
            },
            "weak_evidence_hist": {
                "Description": f"Histogram of number of individuals failing the weak_evidence filter (mutation does not meet likelihood threshold) across heteroplasmy levels, bin edges are: {hl_hist_bin_edges}",
                "Number": "1",
                "Type": "String",
            },
            "contamination_hist": {
                "Description": f"Histogram of number of individuals failing the contamination filter across heteroplasmy levels, bin edges are: {hl_hist_bin_edges}",
                "Number": "1",
                "Type": "String",
            },
            "excluded_AC": {
                "Description": "Excluded allele count (number of individuals in which the variant was filtered out)",
                "Number": "1",
                "Type": "String",
            },
        },
        "format": {
            "GT": {
                "Description": f"Genotype, 1/1 if heteroplasmy level >= {min_hom_threshold}, and 0/1 if heteroplasmy level < {min_hom_threshold}",
                "Number": "1",
                "Type": "String",
            },
            "DP": {
                "Description": "Depth of coverage",
                "Number": "1",
                "Type": "Integer",
            },
            "FT": {
                "Description": "Sample-level filters",
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

    return input_mt, meta_dict, vcf_header_file
