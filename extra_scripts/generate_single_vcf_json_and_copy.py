"""
This script gets only one VCF file ID per sample (based on the most recent
VCF), output this all to JSON and copies those files into the testing project
"""
import argparse
import dxpy as dx
import os
import sys
import warnings

from collections import defaultdict
from pathlib import Path

sys.path.append(os.path.abspath(
    os.path.join(os.path.realpath(__file__), '../../')
))

from utils.file_utils import read_in_csv, read_in_json, write_out_json

warnings.simplefilter(action='ignore', category=FutureWarning)


def parse_args() -> argparse.Namespace:
    """
    Parse the command line arguments inputs given

    Returns
    -------
    args : Namespace
        Namespace of passed command line argument inputs
    """

    parser = argparse.ArgumentParser(
        description='Information necessary to query original VCFs'
    )

    parser.add_argument(
        '-s',
        '--shire_output',
        type=str,
        required=True,
        help='Path to formatted CSV output of reported cases from Shire'
    )

    parser.add_argument(
        '-j',
        '--input_json',
        type=str,
        required=True,
        help='CSV file with two columns containing mapping of X to GM no'
    )

    parser.add_argument(
        '-o',
        '--output_json',
        type=str,
        required=True,
        help='Name of output JSON with one VCF per sample'
    )

    parser.add_argument(
        '--copy-files',
        action=argparse.BooleanOptionalAction,
        help=(
            'If flag given, determined whether or not to copy files to'
            ' DNAnexus project'
        )
    )

    parser.add_argument(
        '-p',
        '--proj_id',
        type=str,
        required=False,
        help=(
            "If copy_files=True, DNAnexus project ID to copy the single VCF"
            " per sample to"
        )
    )

    parser.add_argument(
        '-f',
        '--folder_name',
        type=str,
        required=False,
        help=(
            'If copy_files=True, name of folder to copy VCF files to in'
            ' project'
        )
    )

    args = parser.parse_args()

    return args


def get_only_latest_vcf(vcf_dict):
    """
    In cases where there are multiple VCFs found, take one single VCF
    which is the most recent one. The VCFs were compared and found to be
    identical apart from the header, so this is a way of taking a single VCF
    per sample

    Parameters
    ----------
    vcf_dict : dict
        dict with each sample and the list of VCF files found

    Returns
    -------
    final_dict : dict
        dict with only 1 VCF at most per sample instead of duplicates
    """
    final_dict = defaultdict(list)

    for sample, files in vcf_dict.items():
        # If more than one VCF found
        if len(files) > 1:
            # Get just the dict for the VCF file which was created latest
            recent_vcf = [max(files, key=lambda x:x['created'])]
            # Update the value to be this instead
            final_dict[sample] = recent_vcf[0]
        elif len(files) == 1:
            final_dict[sample] = files[0]
        else:
            # Otherwise, if just one VCF or no VCF, keep this as the value
            final_dict[sample] = {}

    return final_dict


def add_x_number(vcf_dict):
    """
    Add in the X number to the dictionary with each sample and the VCF details

    Parameters
    ----------
    vcf_dict : dict
        dict with each sample and single VCF info as value

    Returns
    -------
    vcf_dict : dict
        final dict with GM number as key and dict with single VCF and sample
        info as value with added X number as nested key
    """
    for sample, vcf_file in vcf_dict.items():
        # If a VCF file has been found, get X number and add as new key
        if vcf_file:
            filename = vcf_file['filename']
            if 'GM' in filename.upper():
                x_number = filename.split('-')[0]
            else:
                x_number = filename.split('_')[0]
            vcf_dict[sample]['x_number'] = x_number

    return vcf_dict


def add_testing_outcome(vcf_dict, outcome_df):
    """
    Add in the testing outcome for the sample

    Parameters
    ----------
    vcf_dict : dict
        dict with sample as key and value is a dict of the one VCF for
        that sample
    outcome_df : pd.DataFrame
        pandas df of each sample and testing outcome and (if relevant)
        the variant(s) found

    Returns
    -------
    vcf_dict : dict
        final dict with GM number as key and dict with single VCF and sample
        info as value
    """
    outcome_df['LABNO'] = outcome_df['LABNO'].str.replace('.','')

    for sample, _ in vcf_dict.items():
        # Get the report outcome from the dataframe based on the GM no
        report_outcome = outcome_df.loc[
            outcome_df['LABNO'] == sample, 'REPORT_TYPE'
        ].item()
        vcf_dict[sample]['report_outcome'] = report_outcome

    return vcf_dict


def copy_files_to_testing_project(vcf_dict, testing_project_id, folder_name):
    """
    Copy the one file per sample to my DNAnexus testing project

    Parameters
    ----------
    vcf_dict : dict
        final dict with GM number as key and dict with single VCF and sample
        info as value
    testing_project_id : str
        project ID for the testing project
    folder_name : str
        name of folder in testing project to copy files to
    """
    for file in vcf_dict.values():
        if file:
            file_id = file['id']
            proj_id = file['project']
            file_object = dx.DXFile(file_id, project=proj_id)
            file_object.clone(testing_project_id, folder=f'/{folder_name}')


def main():
    args = parse_args()
    outcome_df = read_in_csv(args.shire_output)
    sample_vcf_dict = read_in_json(args.input_json)
    vcf_dict = get_only_latest_vcf(sample_vcf_dict)
    vcf_dict_with_x = add_x_number(vcf_dict)
    final_vcf_dict = add_testing_outcome(vcf_dict_with_x, outcome_df)
    write_out_json(args.output_json, final_vcf_dict)
    if args.copy_files:
        copy_files_to_testing_project(
            final_vcf_dict, args.proj_id, args.folder_name
        )

if __name__ == '__main__':
    main()
