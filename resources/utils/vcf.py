"""
Functions related to reading, processing and writing the VCF
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
    these fields to stop potential conflicts with existing INFO fields. Then
    strips the original INFO/CSQ fields

    Parameters
    ----------
    input_vcf : file
        path to VCF file to be split

    Returns
    -------
    output_vcf : str
        name of VCF split by bcftools

    Raises
    ------
    AssertionError
        Raised when non-zero exit code returned by bcftools or output VCF
        is truncated (fewer variants than were in input VCF)
    """

    print(
        f"Splitting necessary fields from {input_vcf} "
        "with bcftools +split-vep"
    )
    output_vcf = f"{Path(input_vcf).stem.split('.')[0]}.split.vcf"

    # Check total rows before splitting out columns
    pre_split_total = subprocess.run(
        f"zgrep -v '^#' {input_vcf} | wc -l",
        shell=True,
        capture_output=True
    )
    # Split out all fields from the CSQ string and name them with
    # 'CSQ_{field}' as separate INFO fields. Remove original CSQ string
    cmd = (
        f"bcftools +split-vep --columns - -a CSQ -Ou -p 'CSQ_'"
        f" -d {input_vcf} | bcftools annotate -x INFO/CSQ -o {output_vcf}"
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

    pre_split = pre_split_total.stdout.decode()
    post_split = post_split_total.stdout.decode()

    print(
        f"Total lines before splitting: {pre_split}"
        f"Total lines after splitting: {post_split}"
    )

    assert int(post_split) >= int(pre_split), (
        "Count of variants following bcftools +split-vep is fewer than the "
        "input VCF"
    )

    return output_vcf


def bcftools_filter(split_vcf, filter_command, filter_vcf):
    """
    Filter the given VCF using bcftools command provided

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
    Raises
    ------
    AssertionError
        Raised when non-zero exit code returned by bcftools or number of
        variants after filtering does not match input VCF
    """
    command = f"{filter_command} {split_vcf} -o {filter_vcf}"

    print(
        f"\nFiltering {split_vcf} with the command: \n\t{command}\n"
    )

    pre_filter_total = subprocess.run(
        f"zgrep -v '^#' {split_vcf} | wc -l",
        shell=True,
        capture_output=True
    )

    output = subprocess.run(command, shell=True, capture_output=True)

    assert output.returncode == 0, (
        f"\n\tError in filtering VCF with bcftools\n"
        f"\n\tVCF: {split_vcf}\n"
        f"\n\tExitcode:{output.returncode}\n"
        f"\n\tbcftools filter command used: {filter_command}\n"
        f"\n\t{output.stderr.decode()}"
    )

    post_filter_total = subprocess.run(
        f"zgrep -v '^#' {filter_vcf} | wc -l",
        shell=True,
        capture_output=True
    )

    pre_filter = pre_filter_total.stdout.decode()
    post_filter = post_filter_total.stdout.decode()

    print(
        f"Total lines before filtering: {pre_filter}"
        f"Total lines after filtering: {post_filter}"
    )

    assert pre_filter == post_filter, (
        "Number of variants before and after bcftools filter do not match"
    )


def read_in_vcf(vcf_file):
    """
    Read in the VCF file with the pysam package and add a new header line
    for the MOI flag

    Parameters
    ----------
    vcf_file : string
        name of the annotated VCF file to be checked against filtering
        parameters

    Returns
    -------
    vcf_contents: pysam VariantFile object
        contents of the VCF file as a pysam object
    sample_name : str
        the full name of the sample
    """
    print(f"Reading in the split VCF {vcf_file} with pysam")

    # Read in and create pysam object of the VCF
    vcf_contents = VariantFile(vcf_file, 'r')

    # Get the name of the sample from the VCF
    sample_name = list(vcf_contents.header.samples)[0]

    # Add MOI as INFO field
    vcf_contents.header.info.add(
        "MOI", ".", "String",
        "Mode of inheritance from PanelApp (simplified)"
    )

    return vcf_contents, sample_name


def add_MOI_field(vcf_contents, panel_dict) -> dict:
    """
    Add MOI INFO field to each variant which will be used for filtering

    Parameters
    ----------
    vcf_contents : pysam.VariantFile object
        pysam object containing all the VCF's info
    panel_dict : dict
        default dict with gene symbol as key and gene info as val

    Returns
    -------
    gene_variant_dict : dict
        dictionary of each gene with value as a list of all variants in that
        gene (plus additional INFO field)
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
                variant.info['MOI'] = 'NONE'

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
        # For each gene, write out each variant to VCF with extra INFO field
        for gene, variant_list in gene_variant_dict.items():
            for variant in variant_list:
                out_vcf.write(variant)


def check_written_out_vcf(
        original_vcf_contents, gene_variant_dict, flagged_vcf
    ):
    """
    Check that the VCF file written out is exactly the same as the header
    and dict of pysam variants that was meant to be written

    Parameters
    ----------
    original_vcf_contents : pysam.VariantFile object
        the pysam object which was to be written to file
    gene_variant_dict : dict
        dictionary with each gene and value as all of the variants in that gene
        as list
    flagged_vcf : file
        the VCF that was written out to be read back in with pysam
    Raises
    ------
    AssertionError
        Raised if the contents of the VCF which was to be written out does
        not match the VCF which was actually written out
    """
    # Read in the VCF (that was just written out) back in with pysam
    flagged_contents = VariantFile(flagged_vcf, 'r')

    # Get list of lines in original header that was written out
    # Get list of variant records which were to be written out
    # Make one big list for original VCF contents which were to be written out
    original_header = list(set(
        str(header) for header in original_vcf_contents.header.records
    ))
    original_records = [
        str(var) for variant_list in gene_variant_dict.values()
        for var in variant_list
    ]
    original_contents = original_header + original_records

    # Get list of lines which were written out to VCF file header
    # Get list of variant records which were in the written out VCF
    # Make one big list for written out VCF
    written_header = list(set(
        str(header) for header in flagged_contents.header.records
    ))
    written_records = [
        str(variant_record) for variant_record in flagged_contents
    ]
    written_contents = written_header + written_records

    assert original_contents == written_contents, (
        "Header and variants written to VCF not identical to those "
        "intended to be written out"
    )


def bcftools_sort(input_vcf):
    """
    Sort the VCF. As we write out variants in gene order, if a
    variant is annotated to multiple genes/transcripts and split, we can end up
    with the VCF not being sorted by position so this fixes that

    Parameters
    ----------
    input_vcf : str
        path to filtered VCF to sort
    output_vcf : str
        name of output sorted VCF

    Outputs
    -------
    {vcf}.sorted.vcf : file
        sorted vcf file output from bcftools
    """
    print(f"Sorting flagged VCF {input_vcf} with bcftools sort")
    output_vcf = f"{Path(input_vcf).stem.split('.')[0]}.sorted.vcf"

    # Check total rows before running bcftools sort
    pre_sort_total = subprocess.run(
        f"zgrep -v '^#' {input_vcf} | wc -l",
        shell=True,
        capture_output=True
    )

    cmd = f"bcftools sort {input_vcf} -o {output_vcf}"

    output = subprocess.run(cmd, shell=True, capture_output=True)

    assert output.returncode == 0, (
        f"\n\tError in sorting VCF with bcftools. VCF: {input_vcf}"
        f"\n\tExitcode:{output.returncode}"
        f"\n\t{output.stderr.decode()}"
    )

    # Check total rows after sorting
    post_sort_total = subprocess.run(
        f"zgrep -v '^#' {output_vcf} | wc -l",
        shell=True,
        capture_output=True
    )

    pre_sort = pre_sort_total.stdout.decode()
    post_sort = post_sort_total.stdout.decode()

    print(
        f"Total lines before sorting: {pre_sort}"
        f"Total lines after sorting: {post_sort}"
    )

    assert pre_sort == post_sort, (
        "Count of variants before and after bcftools sort do not match"
    )


def add_annotation(input_vcf, panel_dict, filter_command):
    """
    Main function to take a VCF and add the INFO field required for filtering

    Parameters
    ----------
    input_vcf : str
        name of the input VCF
    panel_dict : dict
        default dict with each gene on panel as key and the gene info as val
    filter_command : str
        full bcftools filter command
    """
    split_vcf = f"{Path(input_vcf).stem.split('.')[0]}.split.vcf"
    flagged_vcf = f"{Path(input_vcf).stem.split('.')[0]}.flagged.vcf"
    sorted_vcf = f"{Path(input_vcf).stem.split('.')[0]}.sorted.vcf"
    filter_vcf = f"{Path(input_vcf).stem.split('.')[0]}.optimised_filtered.vcf"

    # separate csq fields (creates split_vcf)
    bcftools_pre_process(input_vcf)

    # create pysam object of vcf for flagging
    vcf_contents, sample_name = read_in_vcf(split_vcf)

    # add MOI flags from config
    gene_var_dict = add_MOI_field(vcf_contents, panel_dict)
    write_out_flagged_vcf(flagged_vcf, gene_var_dict, vcf_contents)
    check_written_out_vcf(vcf_contents, gene_var_dict, flagged_vcf)
    bcftools_sort(flagged_vcf)

    # run bcftools filter string from config (create filter_vcf)
    bcftools_filter(sorted_vcf, filter_command, filter_vcf)

    bgzip(filter_vcf)
    os.remove(split_vcf)
    os.remove(flagged_vcf)
    os.remove(sorted_vcf)
    os.remove(filter_vcf)
