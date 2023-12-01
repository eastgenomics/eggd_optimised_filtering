"""
Main script which takes a file and adds our flag for optimised filtering
"""
import argparse
import re

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
        '-f',
        '--filter_string',
        type=str,
        required=True,
        help='BCFtools filter string to be applied'
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

    args = parser.parse_args()

    return args


def check_panel_string(panel_string):
    """
    Check that the panel string given doesn't contain multiple panels
    (HGNC IDs or one panel with extra HGNC IDs is fine)

    Parameters
    ----------
    panel_string : str
        the panel the patient is being tested for
    Raises
    ------
    AssertionError
        Raised if more than one panel is given to the tool (HGNC IDs are OK)
    """
    # Check the number of panels
    panels_from_string = re.sub(
        r'_HGNC:[\d]+(;)?', '', panel_string
    ).rstrip(';')
    panel_counts = panels_from_string.count(';')

    assert panel_counts == 0, (
        f"More than one panel given: {panels_from_string}"
    )

    return panels_from_string


def main():
    args = parse_args()
    filter_string = args.filter_string.replace("\!~", "!~")
    panels_from_string = check_panel_string(args.panel_string)
    bcftools_filter_command = file_utils.unescape_bcftools_command(
        filter_string
    )
    panel_dict = panels.get_formatted_dict(
        panels_from_string, args.genepanels, args.panel_dump
    )
    vcf.add_annotation(
        args.input_vcf, panel_dict, bcftools_filter_command
    )


if __name__ == "__main__":
    main()
