import os
import subprocess

from collections import defaultdict, Counter
from pathlib import Path
from pysam import VariantFile


# Get path of parent directory
ROOT_DIR = Path(__file__).absolute().parents[1]

FIELDS_TO_SPLIT = [
    'SYMBOL', 'Consequence', 'gnomADe_AF', 'gnomADg_AF', 'TWE_AF',
    'ClinVar_CLNSIG', 'ClinVar_CLNSIGCONF', 'SpliceAI_pred_DS_AG',
    'SpliceAI_pred_DS_AL', 'SpliceAI_pred_DS_DG', 'SpliceAI_pred_DS_DL'
]

FIELDS_TO_COLLAPSE = ",".join([
    f"INFO/CSQ_{field}" for field in FIELDS_TO_SPLIT
])


def bgzip(file) -> None:
    """
    Call bgzip on a given file

    Parameters
    ----------
    file : file to compress

    Outputs
    -------
    compressed file

    Raises
    ------
    AssertionError
        Raised when non-zero exit code returned by bgzip
    """
    print(f"Calling bgzip on {file}")

    output = subprocess.run(
        f"bgzip --force {file} -c > {file}.gz",
        shell=True,
        capture_output=True
    )

    assert output.returncode == 0, (
        f"\n\tError in compressing file with bgzip. File: {file}"
        f"\n\tExitcode:{output.returncode}"
        f"\n\t{output.stderr.decode()}"
    )


def bcftools_pre_process(input_vcf) -> str:
    """
    Decompose multiple transcript annotations to individual records, and split VEP CSQ string fields for SYMBOL, Consequence and gnomADe_AF to individual INFO keys. Adds a 'CSQ_' prefix to these fields extracted from the CSQ string to stop potential conflicts with existing INFO fields.

    Parameters
    ----------
    input_vcf : file
        path to VCF file to be split

    Returns
    -------
    output_vcf : str
        name of VCF split by bcftools
    """

    print(
        f"Splitting necessary fields from {input_vcf} "
        "with bcftools +split-vep"
    )
    output_vcf = f"{Path(input_vcf).stem}.split.vcf"

    # Check total rows before splitting out columns
    pre_split_total = subprocess.run(
        f"zgrep -v '^#' {input_vcf} | wc -l",
        shell=True,
        capture_output=True
    )
    # Split out the SYMBOL, Consequence, gnomADe_AF fields from the CSQ
    # string and name them with 'CSQ_{field}' as separate INFO fields
    cmd = (
        f"bcftools +split-vep -c {','.join(FIELDS_TO_SPLIT)} -Ou -p 'CSQ_'"
        f" -d {input_vcf} -o {output_vcf}"
    )

    output = subprocess.run(cmd, shell=True, capture_output=True)

    assert output.returncode == 0, (
        f"\n\tError in splitting VCF with bcftools +split-vep. VCF: {input_vcf}"
        f"\n\tExitcode:{output.returncode}"
        f"\n\t{output.stderr.decode()}"
    )

    # Check total rows after splitting
    post_split_total = subprocess.run(
        f"zgrep -v '^#' {output_vcf} | wc -l",
        shell=True,
        capture_output=True
    )

    print(
        f"Total lines before splitting: {pre_split_total.stdout.decode()}"
        f"Total lines after splitting: {post_split_total.stdout.decode()}"
    )

    return output_vcf


def bcftools_filter(split_vcf, filter_command, filter_vcf):
    """
    Filter given VCF using bcftools command provided in config

    Parameters
    ----------
    split_vcf : pathlib.PosixPath
        path to vcf file to filter
    filter_command : str
        the full bcftools filter command
    filter_vcf : str
        name for output filtered vcf

    Outputs
    -------
    filter_vcf : file
        vcf file with PASS/EXCLUDE added to FILTER columns
    """
    command = f"{filter_command} {split_vcf} -o {filter_vcf}"

    print(
        f"\nFiltering {split_vcf} with the command: \n\t{command}\n"
    )

    output = subprocess.run(command, shell=True, capture_output=True)

    assert output.returncode == 0, (
        f"\n\tError in filtering VCF with bcftools\n"
        f"\n\tVCF: {split_vcf}\n"
        f"\n\tExitcode:{output.returncode}\n"
        f"\n\tbcftools filter command used: {filter_command}\n"
        f"\n\t{output.stderr.decode()}"
    )


def read_in_vcf(filter_vcf, flag_name):
    """
    Read in the VCF file with the pysam package and add a new header line
    for the optimised filtering flag

    Parameters
    ----------
    filename : string
        name of the annotated VCF file to be check against filtering parameters

    Returns
    -------
    vcf_contents: pysam VariantFile object
        contents of the VCF file as a pysam object
    sample_name : str
        the full name of the sample
    """
    print(f"Reading in the split VCF {filter_vcf} with pysam")

    # Read in and create pysam object of the VCF
    vcf_contents = VariantFile(filter_vcf, 'r')

    # Get the name of the sample from the VCF
    sample_name = list(vcf_contents.header.samples)[0]

    # Add in our new flag to VCF header
    vcf_contents.header.info.add(
        flag_name, ".", "String",
        "Flag for optimised filtering based on AF thresholds for a gene's "
        "MOI and zygosity counts"
    )

    # Add the reason the variant was not prioritised as INFO field
    vcf_contents.header.info.add(
        "Filter_reason", ".", "String",
        "Flag explaining why variant has not been prioritised"
    )

    return vcf_contents, sample_name


def add_filtering_flag(
        sample_name, vcf_contents, panel_dict, rules, flag_name
) -> dict:
    """
    Add the flags to each variant which will be used for filtering

    Parameters
    ----------
    sample_name : str
        name of the sample being tested (so we can obtain the genotypes)
    vcf_contents : pysam.VariantFile object
        pysam object containing all the VCF's info
    panel_dict : dict
        default dict with gene symbol as key and gene info as val
    rules : dict
        dict of the filtering rules for each of the inheritance types
    flag_name : str
        name of the info field to add for the flag

    Returns
    -------
    gene_variant_dict : dict
        dictionary of each gene with all variants in that gene as val
        and all filtering flags added
    """
    # Add each variant in a gene/entity to a dict, with gene as key and list
    # of variants as value
    gene_variant_dict = defaultdict(list)
    for record in vcf_contents:
        gene = record.info['CSQ_SYMBOL'][0]
        gene_variant_dict[gene].append(record)

    # For each gene, check whether certain gene info is available
    for gene, variant_list in gene_variant_dict.items():
        gene_present = gene_moi = af_threshold = False
        if gene in panel_dict:
            gene_present = True
            gene_moi = panel_dict[gene].get('mode_of_inheritance')
            if gene_moi:
                if gene_moi in rules:
                    af_threshold = rules[gene_moi].get('af')

        # Iterate over all of the variants called in that gene
        # If variant previously passed bcftools filtering and gene present
        # in our panel_dict and gene_moi and af threshold for that moi
        # are present, check var passes af_threshold for that MOI
        variants_passing_af_filter = []
        for variant in variant_list:
            if variant.filter.keys()[0] == 'PASS':
                if (gene_present) and (gene_moi) and (af_threshold):

                    exome_af = variant.info['CSQ_gnomADe_AF'][0]
                    genome_af = variant.info['CSQ_gnomADg_AF'][0]

                    if not exome_af:
                        exome_af = 0.0
                    if not genome_af:
                        genome_af = 0.0

                    if (exome_af < af_threshold and genome_af < af_threshold):
                        variants_passing_af_filter.append(variant)
                    else:
                        variant.info[flag_name] = 'NOT_PRIORITISED'
                        variant.info['Filter_reason'] = (
                            'gnomAD_AF_exceeds_MOI_threshold'
                        )
                else:
                    variant.info[flag_name] = 'NOT_ASSESSED'
                    variant.info['Filter_reason'] = 'Gene_info_not_available'
            else:
                variant.info[flag_name] = 'NOT_PRIORITISED'
                variant.info['Filter_reason'] = 'BCFtools_filtered'

        # If any variants pass filters for the gene's MOI, get genotypes of all
        if variants_passing_af_filter:
            gtypes = [
                [genotype for genotype in variant.samples[sample_name]['GT']]
                for variant in variants_passing_af_filter
            ]
            # Create list of a count of 1 for each variant
            gt_counts = [genotype.count(1) for genotype in gtypes]
            count_het_homs = Counter(gt_counts)

            # Get number of hets and homs needed for the gene's MOI
            hets_needed = rules[gene_moi]['HET']
            homs_needed = rules[gene_moi]['HOM']

            # Count how many het (1) and hom (2) variants there are which pass
            # filters for that gene
            het_count = count_het_homs.get(1, 0)
            hom_count = count_het_homs.get(2, 0)

            # If either enough hets or enough homs then add our PRIORITISED flag
            if ((het_count >= hets_needed) or (hom_count >= homs_needed)):
                for variant in variants_passing_af_filter:
                    variant.info[flag_name] = 'PRIORITISED'
            # If not enough hets or homs add NOT_PRIORITISED flag
            else:
                for variant in variants_passing_af_filter:
                    variant.info[flag_name] = 'NOT_PRIORITISED'
                    variant.info['Filter_reason'] = (
                        'Zygosity_count_does_not_fit_MOI'
                    )

    return gene_variant_dict


def write_out_flagged_vcf(flagged_vcf, gene_variant_dict, vcf_contents):
    """
    Write out each record to VCF using pysam

    Parameters
    ----------
    flagged_vcf : str
        Name of the VCF to be written out with flags added
    gene_variant_dict : dict
        dictionary with each gene and value as all of the variants in that gene
        as list
    vcf_contents : pysam.VariantFile object
        the contents of the VCF as a pysam object
    """
    print(f"Writing out flagged variants to VCF: {flagged_vcf}")

    with VariantFile(flagged_vcf, 'w', header=vcf_contents.header) as out_vcf:
        # For each gene, write out each of the variants to the VCF with the flag
        for gene, variant_list in gene_variant_dict.items():
            for variant in variant_list:
                out_vcf.write(variant)


def bcftools_remove_csq_annotation(input_vcf):
    """
    Remove expanded CSQ strings which were used to check rules, as having them
    expanded would break eggd_generate_variant workbook and they are still
    present in the VEP CSQ string

    Parameters
    ----------
    input_vcf : str
        Name of VCF file to have CSQ annotations removed

    Returns
    -------
    output_vcf : str
        Name of output VCF file
    """
    output_vcf = f"{Path(input_vcf).stem}.G2P.vcf"

    # check total rows before splitting
    pre_split_total = subprocess.run(
        f"zgrep -v '^#' {input_vcf} | wc -l", shell=True,
        capture_output=True
    )

    cmd = (
        f"bcftools annotate -x {FIELDS_TO_COLLAPSE} {input_vcf} -o {output_vcf}"
    )

    output = subprocess.run(cmd, shell=True, capture_output=True)

    assert output.returncode == 0, (
        f"\n\tError in annotating VCF with bcftools. VCF: {input_vcf}"
        f"\n\tExitcode:{output.returncode}"
        f"\n\t{output.stderr.decode()}"
    )

    # check total rows after splitting
    post_split_total = subprocess.run(
        f"zgrep -v '^#' {output_vcf} | wc -l", shell=True,
        capture_output=True
    )

    print(
        f"Total lines before unsplitting: {pre_split_total.stdout.decode()}"
        f"Total lines after unsplitting: {post_split_total.stdout.decode()}"
    )

    return output_vcf


def add_annotation(
        flag_name, rules, input_vcf, panel_dict, filter_command
    ):
    """
    Main function to take a VCF and add the flags required for filtering

    Parameters
    ----------
    flag_name : str
        Name of the flag to add in
    rules : dict
        dict of the filtering rules for each of the inheritance types
    input_vcf : str
        name of the input VCF
    panel_dict : dict
        default dict with each gene on panel as key and the gene info as val
    """
    split_vcf = f"{Path(input_vcf).stem}.split.vcf"
    filter_vcf = f"{Path(input_vcf).stem}.filter.vcf"
    flagged_vcf = f"{Path(input_vcf).stem}.flagged.vcf"

    bcftools_pre_process(input_vcf)
    bcftools_filter(split_vcf, filter_command, filter_vcf)

    vcf_contents, sample_name = read_in_vcf(filter_vcf, flag_name)
    gene_var_dict = add_filtering_flag(
        sample_name, vcf_contents, panel_dict, rules, flag_name
    )
    write_out_flagged_vcf(flagged_vcf, gene_var_dict, vcf_contents)
    final_vcf = bcftools_remove_csq_annotation(flagged_vcf)
    bgzip(final_vcf)
    os.remove(split_vcf)
    os.remove(filter_vcf)
    os.remove(flagged_vcf)
    os.remove(final_vcf)
