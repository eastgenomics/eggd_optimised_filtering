import pysam as ps
import subprocess
import utils as utils

from collections import defaultdict, Counter
from pathlib import Path


# Get path of parent directory
ROOT_DIR = Path(__file__).absolute().parents[1]



# class VCF():
#     def __init__(self, args) -> None:
#         self.args = args
#         self.het = ['0/1', '1/0']
#         self.hom = '1/1'
#         self.sample_name = None


#     def process(self) -> None:


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


def bcftools_pre_process(input_vcf) -> None:
    """
    Decompose multiple transcript annotations to individual records, and split VEP CSQ string fields for SYMBOL, Consequence and gnomADe_AF to individual INFO keys. Adds a 'CSQ_' prefix to these fields extracted from the CSQ string to stop potential conflicts with existing INFO fields.

    Parameters
    ----------
    input_vcf : str
        path to VCF file to be split

    Returns
    -------
    {vcf}.split.vcf : file
            vcf file output from bcftools
    """

    print(
        f"Splitting necessary fields from {input_vcf}"
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
        f"bcftools +split-vep -c SYMBOL,Consequence,gnomADe_AF -Ou -p 'CSQ_'"
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


def read_in_vcf(split_vcf_gz, flag_name):
    """
    Read in the VCF file with the pysam package

    Parameters
    ----------
    filename : string
        name of the annotated VCF file to be check against filtering parameters

    Returns
    -------
    vcf_contents: pysam VariantFile object
        contents of the VCF file as a pysam object
    """
    # Create pysam object of the VCF
    vcf_contents = ps.VariantFile(split_vcf_gz, 'r')
    sample_name = list(vcf_contents.header.samples)[0]

    # Add in our new flag to VCF header
    vcf_contents.header.info.add(flag_name, ".", "String", "Filtering flag")

    return vcf_contents, sample_name


def add_filtering_flag(sample_name, vcf_contents, panel_dict, rules, csq_types, flag_name):

    gene_variant_dict = defaultdict(list)
    for record in vcf_contents:
        gene = record.info['CSQ_SYMBOL'][0]
        gene_variant_dict[gene].append(record)


    for gene, variant_list in gene_variant_dict.items():
        gene_moi = panel_dict[gene]['mode_of_inheritance']
        af_threshold = rules[gene_moi]['af']

        variants_passing_af_csq_filters = []
        for variant in variant_list:
            consequences = variant.info['CSQ_Consequence'][0].split('&')
            allele_freq = variant.info['CSQ_gnomADe_AF'][0]
            if not allele_freq:
                allele_freq = 0.0

            if (allele_freq < af_threshold and any(x in consequences for x in csq_types)):
                variants_passing_af_csq_filters.append(variant)
            else:
                variant.info[flag_name] = 'EXCLUDE'

        # Get genotypes of all the variants in that gene which pass
        # consequence types and AF filters for that MOI
        gtypes = [
            '/'.join(
                [str(element) for element in variant.samples[sample_name]['GT']]
            ) for variant in variants_passing_af_csq_filters
        ]

        # Get number of hets and homs needed for the gene's MOI
        hets_needed = rules[gene_moi]['HET']
        homs_needed = rules[gene_moi]['HOM']

        # Count how many het and hom variants there are which pass filters
        # for that gene
        counts = Counter(gtypes)
        het_count = counts.get('1/0', 0) + counts.get('0/1', 0)
        hom_count = counts.get('1/1', 0)

        # If either enough hets or enough homs then add our PRIORITY flag
        if ((het_count >= hets_needed) or (hom_count >= homs_needed)):
            for variant in variants_passing_af_csq_filters:
                variant.info[flag_name] = 'PRIORITY'
        # If not enough hets or homs add EXCLUDE flag
        else:
            for variant in variants_passing_af_csq_filters:
                variant.info[flag_name] = 'EXCLUDE'

    return gene_variant_dict


def write_out_flagged_vcf(flagged_vcf, gene_variant_dict, vcf_contents) -> None:
    pass
    with ps.VariantFile(flagged_vcf, 'w', header=vcf_contents.header) as out_vcf:
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
    sample_name : str
        _description_

    Returns
    -------
    _type_
        _description_
    """
    output_vcf = f"{Path(input_vcf).stem}.G2P.vcf"

    # check total rows before splitting
    pre_split_total = subprocess.run(
        f"zgrep -v '^#' {input_vcf} | wc -l", shell=True,
        capture_output=True
    )

    cmd = (
        f"bcftools annotate -x INFO/CSQ_gnomADe_AF,INFO/CSQ_SYMBOL,INFO/CSQ_Consequence {input_vcf} -o {output_vcf}"
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


def add_annotation(flag_name, rules, csq_types, input_vcf, panel_dict):
    # Split and bgzip the annotated VCF
    split_vcf = f"{Path(input_vcf).stem}.split.vcf"
    split_vcf_gz = f"{Path(input_vcf).stem}.split.vcf.gz"
    flagged_vcf = f"{Path(input_vcf).stem}.flagged.vcf"

    bcftools_pre_process(input_vcf)
    bgzip(split_vcf)

    # Read in the split (and bgzipped) VCF with pysam
    vcf_contents, sample_name = read_in_vcf(split_vcf_gz, flag_name)
    gene_var_dict = add_filtering_flag(
        sample_name, vcf_contents, panel_dict, rules, csq_types, flag_name
    )
    write_out_flagged_vcf(flagged_vcf, gene_var_dict, vcf_contents)

    final_vcf = bcftools_remove_csq_annotation(flagged_vcf)
    bgzip(final_vcf)


if __name__ == "__main__":
    pass
