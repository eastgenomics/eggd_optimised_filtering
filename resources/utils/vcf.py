"""
Functions related to processing and reading in the VCF
"""
import os
import subprocess

from collections import defaultdict
from pathlib import Path
from pysam import VariantFile


# Get path of parent directory
ROOT_DIR = Path(__file__).absolute().parents[1]


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
    Decompose multiple transcript annotations to individual records, and split
    VEP CSQ string fields to individual INFO keys. Adds a 'CSQ_' prefix to
    these fields to stop potential conflicts with existing INFO fields.

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
        f"bcftools +split-vep --columns - -a CSQ -Ou -p 'CSQ_'"
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


def read_in_vcf(vcf_file):
    """
    Read in the VCF file with the pysam package and add a new header line
    for the MOI flag

    Parameters
    ----------
    vcf_file : string
        name of the annotated VCF file to be check against filtering parameters

    Returns
    -------
    vcf_contents: pysam VariantFile object
        contents of the VCF file as a pysam object
    sample_name : str
        the full name of the sample
    csq_fields_to_collapse : str
        string containing VEP split INFO fields to remove so we don't break
        eggd_generate_variant_workbooks which requires unsplit VEP fields
    """
    print(f"Reading in the split VCF {vcf_file} with pysam")

    # Read in and create pysam object of the VCF
    vcf_contents = VariantFile(vcf_file, 'r')

    # Get the name of the sample from the VCF
    sample_name = list(vcf_contents.header.samples)[0]

    # Get a list of each of the VEP fields present in the VCF from the header
    csq_string_list = vcf_contents.header.info['CSQ'].description.split(
        'Format: '
    )[1].split('|')
    # Create single comma separated string of all of them with INFO/CSQ_ prefix
    csq_fields_to_collapse = ",".join([
        f"INFO/CSQ_{field}" for field in csq_string_list
    ])

    # Add MOI as INFO field
    vcf_contents.header.info.add(
        "MOI", ".", "String",
        "Mode of inheritance from PanelApp (simplified)"
    )

    return vcf_contents, sample_name, csq_fields_to_collapse


def add_filtering_flag(vcf_contents, panel_dict) -> dict:
    """
    Add the flags to each variant which will be used for filtering

    Parameters
    ----------
    vcf_contents : pysam.VariantFile object
        pysam object containing all the VCF's info
    panel_dict : dict
        default dict with gene symbol as key and gene info as val

    Returns
    -------
    gene_variant_dict : dict
        dictionary of each gene with all variants in that gene as val
        and all filtering flags added
    """
    # Add each variant in a gene/entity to a dict, with gene as key and list
    # of variants as value
    gene_variant_dict = defaultdict(list)
    for variant_record in vcf_contents:
        gene = variant_record.info['CSQ_SYMBOL'][0]
        gene_variant_dict[gene].append(variant_record)

    # For each gene, check whether certain gene info is available
    for gene, variant_list in gene_variant_dict.items():
        gene_present = gene_moi = False
        if gene in panel_dict:
            gene_present = True
            gene_moi = panel_dict[gene].get('mode_of_inheritance')

        # Iterate over all of the variants called in that gene
        # If gene present in our panel_dict and gene_moi is present,
        # add the MOI we've taken from PanelApp to variant info
        # otherwise set it to unknown
        for variant in variant_list:
            if all([gene_present, gene_moi]):
                variant.info['MOI'] = gene_moi
            else:
                variant.info['MOI'] = 'UNKNOWN'

    return gene_variant_dict


def write_out_flagged_vcf(flagged_vcf, gene_variant_dict, vcf_contents):
    """
    Write out each variant record to VCF using pysam

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


def bcftools_remove_csq_annotation(input_vcf, fields_to_collapse):
    """
    Remove expanded CSQ strings which were used to check rules, as having them
    expanded would break eggd_generate_variant workbook and they are still
    present in the VEP CSQ string

    Parameters
    ----------
    input_vcf : str
        Name of VCF file to have CSQ annotations removed
    fields_to_collapse : str
        string containing VEP split INFO fields to remove so we don't break
        eggd_generate_variant_workbooks which requires unsplit VEP CSQ fields

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
        f"bcftools annotate -x {fields_to_collapse} {input_vcf} -o {output_vcf}"
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


def add_annotation(input_vcf, panel_dict, filter_command):
    """
    Main function to take a VCF and add the flags required for filtering

    Parameters
    ----------
    input_vcf : str
        name of the input VCF
    panel_dict : dict
        default dict with each gene on panel as key and the gene info as val
    filter_command : str
        full bcftools filter command
    """
    split_vcf = f"{Path(input_vcf).stem}.split.vcf"
    filter_vcf = f"{Path(input_vcf).stem}.filter.vcf"
    flagged_vcf = f"{Path(input_vcf).stem}.flagged.vcf"

    # separate csq fields (creates split_vcf)
    bcftools_pre_process(input_vcf)

    # create pysam object of vcf for flagging
    vcf_contents, sample_name, csq_fields_to_drop = read_in_vcf(split_vcf)

    # add MOI flags from config
    gene_var_dict = add_filtering_flag(
        vcf_contents, panel_dict
    )
    write_out_flagged_vcf(flagged_vcf, gene_var_dict, vcf_contents)

    # run bcftools filter string from config (create filter_vcf)
    bcftools_filter(flagged_vcf, filter_command, filter_vcf)
    # return csq fields to standard format
    final_vcf = bcftools_remove_csq_annotation(filter_vcf, csq_fields_to_drop)
    bgzip(final_vcf)
    os.remove(split_vcf)
    os.remove(filter_vcf)
    os.remove(flagged_vcf)
    os.remove(final_vcf)
