import argparse
import json
import pysam as ps
import subprocess

from pathlib import Path
from utils import panels as pq


# Get path of directory this script is in
ROOT_DIR = Path(__file__).absolute().parents[0]

# Set up command line args
parser = argparse.ArgumentParser(
    description='Information necessary for optimised filtering'
)

# Add CLI arg of the input annotated VCF to have filtering flags added
parser.add_argument(
    '-i',
    '--input_vcf',
    type=str,
    help="VCF file to be filtered"
)

# Add PanelApp code of panel of interest
parser.add_argument(
    '-p',
    '--panelapp_id',
    type=str,
    help="PanelApp R code for panel of interest"
)


def read_args(args):
    """
    Read in command line arguments

    Parameters
    ----------
    args : argparse.Namespace class
        argpase object

    Returns
    -------
    vcf_filename : str
        name of input annotated VCF to have filtering flags added
    panelapp_id : int
        ID of PanelApp panel
    """
    vcf_filename = args.input_vcf
    panelapp_id = args.panelapp_id

    return vcf_filename, panelapp_id


def bgzip(file) -> None:
    """
    Call bgzip on given file

    Parameters
    ----------
    file : file to compress

    Outputs
    -------
    input file, but compressed

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

    # check total rows before splitting
    pre_split_total = subprocess.run(
        f"zgrep -v '^#' {input_vcf} | wc -l",
        shell=True,
        capture_output=True
    )

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

    # check total rows after splitting
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


def read_in_vcf(input_vcf):
    """
    Read in the VCF file with the pysam package

    Parameters
    ----------
    filename : string
        name of the VCF file to be annotated

    Returns
    -------
    vcf_contents: pysam VariantFile object
        contents of the VCF file as a pysam object
    """
    vcf_filename = f"{Path(input_vcf).stem}.split.vcf"
    vcf_contents = ps.VariantFile(vcf_filename, 'r')

    return vcf_contents


def set_up_output_vcf(input_vcf, vcf_contents):
    """
    Sets up the output VCF and header to include the filter flag annotation

    Parameters
    ----------
    vcf_contents : pysam VariantFile object
        contents of the VCF file as a pysam object

    Returns
    -------
    sample_name : str
        Full name of the sample within the VCF
    vcf_out : pysam VariantFile object
        VCF file which will be written to, containing original header
        and extra header line for the INFO flag
    """
    # Get names of samples from the header
    samples = list(vcf_contents.header.samples)
    # Check only one sample is present
    assert len(samples) == 1, (
        f"More than one sample found: {samples}"
    )
    # Get the name of the single sample
    sample_name = samples[0]
    # Add a line in the VCF header for the filtering flag
    vcf_contents.header.info.add('My_Flag', ".", "String", "Filtering flag")
    # Set up an output VCF with the extra line in the header to write to
    outfile = f"{Path(input_vcf).stem}_annotated_intermediate.vcf"
    vcf_out = ps.VariantFile(
        outfile,
        'w',
        header=vcf_contents.header
    )
    return sample_name, vcf_out


def read_in_rules():
    """
    Read in the filtering rules and consequence types from JSON file

    Returns
    -------
    rules : dict
        dictionary containing filtering rules for each inheritance type
    csq_types : list
        list of relevant variant consequence types
    """
    with open(
        ROOT_DIR.joinpath(f"{ROOT_DIR}/resources/gene_rules.json"), 'r', encoding='utf8'
    ) as json_file:
        content = json.load(json_file)

    rules = content.get('rules')
    csq_types = content.get('csq_types')

    return rules, csq_types


def add_flag_to_variants(
        vcf_contents, rules, csq_types, name_of_sample, vcf_out, obesity_dict
    ):
    """
    Add a flag to the CSQ_Flag field if variant passes filters

    Parameters
    ----------
    sample_name : str
        Full name of the sample within the VCF
    vcf_out : pysam VariantFile object
        VCF file which will be written to, containing original header
        and extra header line for the INFO flag
    """
    for record in vcf_contents:
        # Get the gene the variant affects
        gene = record.info['CSQ_SYMBOL'][0]
        # Get the consequence type of the variant
        consequence = record.info['CSQ_Consequence'][0]
        # Get the gnomAD pop allele frequency
        allele_freq = record.info['CSQ_gnomADe_AF'][0]
        # If there is no gnomAD allele freq, convert None to 0.0 so can
        # compare later
        if not allele_freq:
            allele_freq = 0.0
        # Get the genotype of the sample as a '/' separated string
        # instead of a tuple
        gtype = '/'.join(
            [str(element) for element in record.samples[name_of_sample]['GT']]
        )
        # Get mode of inheritance requirement for the gene so we can check
        # the variant rules for that mode
        allelic_req = obesity_dict[gene].get('mode_of_inheritance')
        # If gene is in panel being assessed
        if gene in obesity_dict:
            # If there is a mode of inheritance for the gene/region present
            if allelic_req:
                # If variant satisfies our rules, add FILTER_FLAG
                if (
                    (allele_freq < rules[allelic_req]['af'])
                    and (gtype in rules[allelic_req]['GT'])
                    and (consequence in csq_types)
                ):
                    record.info['My_Flag'] = 'FILTER_FLAG'
                else:
                    record.info['My_Flag'] = 'NO_FILTER'
            else:
                record.info['My_Flag'] = 'MOI_NOT_PRESENT'
        else:
            record.info['My_Flag'] = 'GENE_NOT_PRESENT'

        # Write the variant to the output VCF
        vcf_out.write(record)
    vcf_out.close()
    print("Written out VCF")


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
    input_vcf = f"{Path(input_vcf).stem}_annotated_intermediate.vcf"
    output_vcf = f"{Path(input_vcf).stem}.filter_annotated.vcf"

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


def main():
    args = parser.parse_args()
    input_vcf, panelapp_id = read_args(args)
    print(input_vcf)
    rules, csq_types = read_in_rules()
    bcftools_pre_process(input_vcf)
    vcf_contents = read_in_vcf(input_vcf)
    sample_name, vcf_out = set_up_output_vcf(input_vcf, vcf_contents)
    print(vcf_out)
    obesity_dict = pq.get_formatted_dict(panelapp_id)
    add_flag_to_variants(
        vcf_contents, rules, csq_types, sample_name, vcf_out, obesity_dict
    )
    bcftools_remove_csq_annotation(input_vcf)

if __name__ == "__main__":
    main()
