import hail as hl


def add_descriptions(
    input_mt: hl.MatrixTable,
    min_hom_threshold: float = 0.95,
    vaf_filter_threshold: float = 0.01,
    min_het_threshold: float = 0.10,
) -> hl.MatrixTable:
    """
    Add descriptions of annotations to globals.

    :param input_mt: MatrixTable
    :param min_hom_threshold: minimum cutoff to define a variant as homoplasmic
    :param min_het_threshold: minimum cutoff to define a variant as a PASS heteroplasmic variant, genotypes below this threshold will count towards the below_min_het_threshold
    :param vaf_filter_threshold: should match vaf_filter_threshold supplied to Mutect2, variants below this value will be set to homoplasmic reference after calculating the common_low_heteroplasmy filter    :return: MatrixTable with descriptions of annotations stored in  globals
    :rtype: MatrixTable

    """
    # pull out variables that are needed in description dictionaries
    hap_order = hl.eval(input_mt.hap_order)
    population_order = hl.eval(input_mt.pop_order)
    hl_hist_bin_edges = input_mt.hl_hist.bin_edges.take(1)[0]

    global_annotation_dict = hl.struct(
        vep_version=hl.struct(Description="VEP version"),
        dbsnp_version=hl.struct(Description="Version of dbSNP used"),
        hap_order=hl.struct(
            Description="The order in which haplogroups are reported for haplogroup-related annotations"
        ),
        dp_hist_all_variants_bin_freq=hl.struct(
            Description="Histogram values for depth (DP) across all variants"
        ),
        dp_hist_all_variants_n_larger=hl.struct(
            Description="Number of values greater than largest histogram bin for depth (DP) across all variants"
        ),
        dp_hist_all_variants_bin_edges=hl.struct(
            Description="Values of bin edges for histogram of depth (DP) across all variants"
        ),
        mq_hist_all_variants_bin_freq=hl.struct(
            Description="Histogram values for mapping quality (MQ) across all variants"
        ),
        mq_hist_all_variants_n_larger=hl.struct(
            Description="Number of values greater than largest histogram bin for mapping quality (MQ) across all variants"
        ),
        mq_hist_all_variants_bin_edges=hl.struct(
            Description="Values of bin edges for histogram of mapping quality (MQ) across all variants"
        ),
        tlod_hist_all_variants_bin_freq=hl.struct(
            Description="Histogram values for tumor log odds (TLOD) across all variants"
        ),
        tlod_hist_all_variants_n_larger=hl.struct(
            Description="Number of values greater than largest histogram bin for tumor log odds (TLOD) across all variants"
        ),
        tlod_hist_all_variants_bin_edges=hl.struct(
            Description="Values of bin edges for histogram of tumor log odds (TLOD) across all variants"
        ),
        age_hist_all_samples_bin_freq=hl.struct(
            Description="Histogram values for ages across all samples"
        ),
        age_hist_all_samples_n_larger=hl.struct(
            Description="Number of values greater than largest histogram bin for ages across all samples"
        ),
        age_hist_all_samples_n_smaller=hl.struct(
            Description="Number of values less than smallest histogram bin for ages across all samples"
        ),
        age_hist_all_samples_bin_edges=hl.struct(
            Description="Values of bin edges for histogram of ages across all samples"
        ),
        pop_order=hl.struct(
            Description="The order in which populations are reported for population-related annotations"
        ),
    )

    col_annotation_dict = hl.struct(
        s=hl.struct(Description="Sample ID"),
        participant_id=hl.struct(
            Description="Participant ID which is used on the terra platform, not always the same as sample ID"
        ),
        contamination=hl.struct(
            Description="Estimate of mitochondrial contamination that is output by Haplocheck, reported as a proportion"
        ),
        freemix_percentage=hl.struct(
            Description="Estimate of nuclear contamination that is output by verifyBamID, reported as a percentage"
        ),
        major_haplogroup=hl.struct(
            Description="The major haplogroup that is output by Haplogrep"
        ),
        wgs_median_coverage=hl.struct(
            Description="The median depth of coverage on the autosomes for whole genome sequencing data which is output by Picard’s CollectWgsMetrics tool"
        ),
        mt_mean_coverage=hl.struct(
            Description="Mean mitochondrial depth of coverage that is output by MuTect2"
        ),
        hap=hl.struct(
            Description="The top level haplogroup assigned to the samples, typically the first character or two of the major haplogroup"
        ),
        release=hl.struct(Description="True if the sample is in the gnomAD release"),
        hard_filters=hl.struct(
            Description="A set containing any of the hard filters that were applied to the sample for QC on the nuclear genome"
        ),
        research_project=hl.struct(
            Description="Description of the research project to which the sample belongs"
        ),
        project_id=hl.struct(
            Description="The Project ID for the sample, typically a RP or G project for internal Broad samples or a short description for external samples"
        ),
        product=hl.struct(Description="Sequencing platform"),
        sample_pi=hl.struct(
            Description="Principal investigator of the research project of the sample"
        ),
        sex_karyotype=hl.struct(
            Description="Sample’s sex karyotype (combined X and Y karyotype)"
        ),
        age=hl.struct(Description="Age of the sample"),
        broad_external=hl.struct(
            Description='Whether the sample was sequenced internally ("broad") or externally ("external")'
        ),
        pop=hl.struct(
            Description="The inferred population (nuclear ancestry) of the sample"
        ),
        mito_cn=hl.struct(
            Description="The estimated mitochondrial copy number of the sample, calculated as 2*mt_mean_coverage/wgs_median_coverage"
        ),
        over_85_mean=hl.struct(
            Description="Mean heteroplasmy level (restricted to heteroplasmy levels >= 0.85) for the sample for PASS variants that are haplogroup-defining variants and are not at an artifact-prone site"
        ),
        over_85_count=hl.struct(
            Description="Number of variants with heteroplasmy levels >= 0.85 for the sample for PASS variants that are haplogroup-defining variants and are not at an artifact-prone site"
        ),
        bt_85_and_99_mean=hl.struct(
            Description="Mean heteroplasmy level (restricted to heteroplasmy levels >= 0.85 and <= 0.998) for the sample for PASS variants that are haplogroup-defining variants and are not at an artifact-prone site"
        ),
        bt_85_and_99_count=hl.struct(
            Description="Number of variants with heteroplasmy levels >= 0.85 and <= 0.998 for the sample for PASS variants that are haplogroup-defining variants and are not at an artifact-prone site"
        ),
        contam_high_het=hl.struct(
            Description="Internal estimate of contamination. It is defined for each sample as one minus the mean heteroplasmy of any haplogroup-defining variants observed with heteroplasmy 85-99.8 percent alternate alleles if at least 3 such variants are present; otherwise estimated contamination is defined as one minus the mean heteroplasmy of haplogroup-defining variants with heteroplasmy 85-100 percent"
        ),
        indel_pos_counter=hl.struct(
            Description="Dictionary containing the count of the number of indels at each position that contains at least one indel for the sample"
        ),
    )

    row_annotation_dict = hl.struct(
        locus=hl.struct(
            Description="Hail LocusExpression containing contig and position"
        ),
        alleles=hl.struct(
            Description="Alternate allele (multiallelic sites are split)"
        ),
        rsid=hl.struct(Description="The RSID obtained from dbSNP"),
        filters=hl.struct(
            Description="Site or allele-specific filters applied to the variant"
        ),
        variant_collapsed=hl.struct(Description="Variant in format of RefPosAlt"),
        hap_defining_variant=hl.struct(
            Description="Present if variant is present as a haplogroup defining variant in PhyloTree build 17"
        ),
        pon_mt_trna_prediction=hl.struct(
            Description="tRNA pathogenicity classification from PON-mt-tRNA"
        ),
        pon_ml_probability_of_pathogenicity=hl.struct(
            Description="tRNA ML_probability_of_pathogenicity from PON-mt-tRNA"
        ),
        mitotip_score=hl.struct(Description="MitoTip raw score"),
        mitotip_trna_prediction=hl.struct(Description="MitoTip score interpretation"),
        region=hl.struct(
            Description="Region (control, non-control, or NA) as obtained from the file supplied in add_variant_context"
        ),
        variant_context=hl.struct(
            Description="Variant in format of REF>ALT with strand (heavy or light) appended"
        ),
        vep=hl.struct(Description="Annotations from VEP"),
        common_low_heteroplasmy=hl.struct(
            Description=f"Flag, true if the overall allele frequency is > 0.001 for samples with a heteroplasmy > 0 and < 0.50 (includes variants < {vaf_filter_threshold} heteroplasmy which are subsequently filtered), is evaluated before other sample-level filters are applied"
        ),
        base_qual_hist=hl.struct(
            Description=f"Histogram of number of individuals failing the base_qual filter (alternate allele median base quality) across heteroplasmy levels, bin edges are: {hl_hist_bin_edges}"
        ),
        position_hist=hl.struct(
            Description=f"Histogram of number of individuals failing the position filter (median distance of alternate variants from end of reads) across heteroplasmy levels, bin edges are: {hl_hist_bin_edges}"
        ),
        strand_bias_hist=hl.struct(
            Description=f"Histogram of number of individuals failing the strand_bias filter (evidence for alternate allele comes from one read direction only) across heteroplasmy levels, bin edges are: {hl_hist_bin_edges}"
        ),
        weak_evidence_hist=hl.struct(
            Description=f"Histogram of number of individuals failing the weak_evidence filter (mutation does not meet likelihood threshold) across heteroplasmy levels, bin edges are: {hl_hist_bin_edges}"
        ),
        contamination_hist=hl.struct(
            Description=f"Histogram of number of individuals failing the contamination filter across heteroplasmy levels, bin edges are: {hl_hist_bin_edges}"
        ),
        heteroplasmy_below_min_het_threshold_hist=hl.struct(
            Description=f"Histogram of number of individuals with a heteroplasmy level below {min_het_threshold}, bin edges are: {hl_hist_bin_edges}"
        ),
        excluded_AC=hl.struct(
            Description="Excluded allele count (number of individuals in which the variant was filtered out)"
        ),
        AN=hl.struct(
            Description="Overall allele number (number of samples with non-missing genotype)"
        ),
        AC_hom=hl.struct(
            Description=f"Allele count restricted to variants with a heteroplasmy level >= {min_hom_threshold}"
        ),
        AC_het=hl.struct(
            Description=f"Allele count restricted to variants with a heteroplasmy level >= 0.10 and < {min_hom_threshold}"
        ),
        hl_hist=hl.struct(
            Description="Histogram of heteroplasmy levels",
            sub_annotations=hl.struct(
                bin_edges=hl.struct(Description="Bin edges of the histogram"),
                bin_freq=hl.struct(
                    Description="Count for number of values in each bin"
                ),
                n_smaller=hl.struct(
                    Description="Number of values smaller than the start of the first bin"
                ),
                n_larger=hl.struct(
                    Description="Number of values larger than the end of the last bin"
                ),
            ),
        ),
        dp_hist_all=hl.struct(
            Description="Histogram of dp values for all individuals",
            sub_annotations=hl.struct(
                bin_edges=hl.struct(Description="Bin edges of the histogram"),
                bin_freq=hl.struct(
                    Description="Count for number of values in each bin"
                ),
                n_smaller=hl.struct(
                    Description="Number of values smaller than the start of the first bin"
                ),
                n_larger=hl.struct(
                    Description="Number of values larger than the end of the last bin"
                ),
            ),
        ),
        dp_hist_alt=hl.struct(
            Description="Histogram of dp values for individuals with the alternative allele",
            sub_annotations=hl.struct(
                bin_edges=hl.struct(Description="Bin edges of the histogram"),
                bin_freq=hl.struct(
                    Description="Count for number of values in each bin"
                ),
                n_smaller=hl.struct(
                    Description="Number of values smaller than the start of the first bin"
                ),
                n_larger=hl.struct(
                    Description="Number of values larger than the end of the last bin"
                ),
            ),
        ),
        dp_mean=hl.struct(Description="Mean depth across all individuals for the site"),
        mq_mean=hl.struct(
            Description="Mean MMQ (median mapping quality) across individuals with a variant for the site"
        ),
        tlod_mean=hl.struct(
            Description="Mean TLOD (Log 10 likelihood ratio score of variant existing versus not existing) across individuals with a variant for the site"
        ),
        AF_hom=hl.struct(
            Description=f"Allele frequency restricted to variants with a heteroplasmy level >= {min_hom_threshold}"
        ),
        AF_het=hl.struct(
            Description=f"Allele frequency restricted to variants with a heteroplasmy level >= 0.10 and < {min_hom_threshold}"
        ),
        max_hl=hl.struct(
            Description="Maximum heteroplasmy level observed among all samples with the variant"
        ),
        hap_AN=hl.struct(
            Description=f"List of overall allele number for each haplogroup, haplogroup order: {hap_order}"
        ),
        hap_AC_het=hl.struct(
            Description=f"List of AC_het for each haplogroup, haplogroup order: {hap_order}"
        ),
        hap_AC_hom=hl.struct(
            Description=f"List of AC_hom for each haplogroup, haplogroup order: {hap_order}"
        ),
        hap_AF_hom=hl.struct(
            Description=f"List of AF_hom for each haplogroup, haplogroup order: {hap_order}"
        ),
        hap_AF_het=hl.struct(
            Description=f"List of AF_het for each haplogroup, haplogroup order: {hap_order}"
        ),
        hap_hl_hist=hl.struct(
            Description=f"Histogram of heteroplasmy levels for each haplogroup; bin edges are: {hl_hist_bin_edges}, haplogroup order: {hap_order}"
        ),
        hap_faf_hom=hl.struct(
            Description=f"List of filtering allele frequency for each haplogroup restricted to homoplasmic variants, haplogroup order: {hap_order}"
        ),
        hapmax_AF_hom=hl.struct(Description="Haplogroup with maximum AF_hom"),
        hapmax_AF_het=hl.struct(Description="Haplogroup with maximum AF_het"),
        faf_hapmax_hom=hl.struct(
            Description="Maximum filtering allele frequency across haplogroups restricted to homoplasmic variants"
        ),
        pop_AN=hl.struct(
            Description=f"List of overall allele number for each population, population order: {population_order}"
        ),
        pop_AC_het=hl.struct(
            Description=f"List of AC_het for each population, population order: {population_order}",
        ),
        pop_AC_hom=hl.struct(
            Description=f"List of AC_hom for each population, population order: {population_order}"
        ),
        pop_AF_hom=hl.struct(
            Description=f"List of AF_hom for each population, population order: {population_order}"
        ),
        pop_AF_het=hl.struct(
            Description=f"List of AF_het for each population, population order: {population_order}"
        ),
        pop_hl_hist=hl.struct(
            Description=f"Histogram of heteroplasmy levels for each population; bin edges are: {hl_hist_bin_edges}, population order: {population_order}"
        ),
        age_hist_hom=hl.struct(
            Description="Histogram of ages of individuals with a homoplasmic variant",
            sub_annotations=hl.struct(
                bin_edges=hl.struct(Description="Bin edges of the histogram"),
                bin_freq=hl.struct(
                    Description="Count for number of values in each bin"
                ),
                n_smaller=hl.struct(
                    Description="Number of values smaller than the start of the first bin"
                ),
                n_larger=hl.struct(
                    Description="Number of values larger than the end of the last bin"
                ),
            ),
        ),
        age_hist_het=hl.struct(
            Description="Histogram of ages of individuals with a heteroplasmic variant",
            sub_annotations=hl.struct(
                bin_edges=hl.struct(Description="Bin edges of the histogram"),
                bin_freq=hl.struct(
                    Description="Count for number of values in each bin"
                ),
                n_smaller=hl.struct(
                    Description="Number of values smaller than the start of the first bin"
                ),
                n_larger=hl.struct(
                    Description="Number of values larger than the end of the last bin"
                ),
            ),
        ),
    )

    # add descriptions to the matrix table
    input_mt = input_mt.annotate_globals(
        global_annotation_descriptions=hl.literal(global_annotation_dict),
        col_annotation_descriptions=hl.literal(col_annotation_dict),
        row_annotation_descriptions=hl.literal(row_annotation_dict),
    )

    return input_mt


def adjust_descriptions(input_ht: hl.Table) -> hl.Table:
    """
    Remove descriptions for globals that are no longer present in the ht.

    :param input_ht: Hail Table of variants
    :return: Hail Table with globals descriptions edited to remove descriptions for annotations no longer present in the table
    :rtype: hl.Table
    """
    row_descriptions = [x for x in input_ht.row_annotation_descriptions.keys()]
    row_fields = list(input_ht.row)
    rows_descriptions_to_drop = set(row_descriptions) - set(row_fields)

    global_descriptions = [x for x in input_ht.global_annotation_descriptions.keys()]
    gobal_fields = list(input_ht.globals)
    global_descriptions_to_drop = set(global_descriptions) - set(gobal_fields)

    input_ht = input_ht.annotate_globals(
        row_annotation_descriptions=input_ht.row_annotation_descriptions.drop(
            *rows_descriptions_to_drop
        ),
        global_annotation_descriptions=input_ht.global_annotation_descriptions.drop(
            *global_descriptions_to_drop
        ),
    )
    input_ht = input_ht.drop("col_annotation_descriptions")

    return input_ht
