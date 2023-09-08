"""
Main script which takes a file and adds INFO field flag for optimised filtering
"""
import argparse

from utils import vcf
from utils import panels
from utils import file_utils


def parse_args() -> argparse.Namespace:
    """
    Parse the command line arguments inputs given

    Returns
    -------
    args : Namespace
        Namespace of passed command line argument inputs
    """

    parser = argparse.ArgumentParser(
        description='Information necessary for optimised filtering'
    )

    # Add CLI arg of the input annotated VCF to have filtering flags added
    parser.add_argument(
        '-i',
        '--input_vcf',
        type=str,
        required=True,
        help='Annotated VCF file to have optimised filtering applied'
    )

    # Add CLI input of clinical indication string
    parser.add_argument(
        '-p',
        '--panel_string',
        type=str,
        required=True,
        help='String containing the panel(s) of interest, comma separated'
    )

    parser.add_argument(
        '-c',
        '--config',
        type=str,
        required=True,
        help='Config file containing filtering rules'
    )

    parser.add_argument(
        '-g',
        '--genepanels',
        type=str,
        required=True,
        help="genepanels file with panel IDs included"
    )

    parser.add_argument(
        '-d',
        '--panel_dump',
        type=str,
        required=True,
        help="PanelApp JSON dump"
    )

    parser.add_argument(
        '-w',
        '--whitelist',
        type=str,
        required=False,
        help='VCF of variants to always give back'
    )

    args = parser.parse_args()

    return args


def read_in_config(file_path):
    """
    Read in the info needed for filtering from JSON config

    Parameters
    ----------
    file_path : str
        path of the config file

    Returns
    -------
    flag_name : str
        name of the flag to be added
    rules : dict
        dict of the filtering rules for each gene MOI
    VEP_fields_to_split : list
        list of VEP fields to split with bcftools
    bcftools_filter_string : str
        bcftools command as a string
    """
    config_contents = file_utils.read_in_json(file_path)

    return list(map(config_contents.get, [
        'flag_name', 'filtering_rules', 'VEP_fields_to_split',
        'bcftools_filter_string'
    ]))


def main():
    args = parse_args()
    flag_name, rules, fields_to_split, filter_string = read_in_config(
        args.config
    )
    bcftools_filter_command = file_utils.unescape_bcftools_command(
        filter_string
    )
    panel_dict = panels.get_formatted_dict(
        args.panel_string, args.genepanels, args.panel_dump
    )
    vcf.add_annotation(
        flag_name, rules, fields_to_split, args.input_vcf, panel_dict,
        bcftools_filter_command
    )


if __name__ == "__main__":
    main()
