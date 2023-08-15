"""
Main script which takes a file and adds our flag for optimised filtering
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
        '--panelapp_string',
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
        '-w',
        '--whitelist',
        type=str,
        required=False,
        help='VCF of variants to always give back'
    )

    args = parser.parse_args()

    return args


def read_in_config(config_file_id):
    """
    Read in the info needed for filtering from config JSON

    Parameters
    ----------
    config_file_id : str
        proj:file ID of the config file in DNAnexus

    Returns
    -------
    flag_name : str
        name of the flag to be added
    panelapp_file : str
        file ID of the PanelApp dump in DNAnexus
    genepanels_file : str
        file ID of the genepanels file in DNAnexus
    rules : dict
        dict of the filtering rules for each gene MOI
    csq_types : list
        list of variant consequence types we want to keep
    """
    config_contents = file_utils.read_in_json_from_dnanexus(config_file_id)
    flag_name = config_contents.get('flag_name')
    panelapp_file = config_contents.get('panelapp_file_id')
    genepanels_file = config_contents.get('genepanels_file_id')
    rules = config_contents.get('filtering_rules')
    csq_types = config_contents.get('csq_types')

    return flag_name, panelapp_file, genepanels_file, rules, csq_types


def main():
    args = parse_args()
    input_vcf = args.input_vcf
    panel_string = args.panelapp_string
    config = args.config
    whitelist = args.whitelist
    flag_name, panelapp_file, genepanels_file, rules, csq_types = read_in_config(config)

    # Create dictionary from the panel
    panel_dict = panels.get_formatted_dict(
        panel_string, genepanels_file, panelapp_file
    )
    vcf.add_annotation(flag_name, rules, csq_types, input_vcf, panel_dict)


if __name__ == "__main__":
    main()
