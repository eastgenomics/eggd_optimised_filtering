"""
Script to unarchive the original raw VCF from Sentieon for each sample
"""
import argparse
import dxpy as dx
import os
import sys
import time

sys.path.append(os.path.abspath(
    os.path.join(os.path.realpath(__file__), '../../')
))

from utils.file_utils import read_in_json


def parse_args() -> argparse.Namespace:
    """
    Parse the command line arguments inputs given

    Returns
    -------
    args : Namespace
        Namespace of passed command line argument inputs
    """

    parser = argparse.ArgumentParser(
        description='Information for unarchiving files'
    )

    parser.add_argument(
        '-j',
        '--vcf_json',
        type=str,
        required=True,
        help='JSON containing one file for each sample to unarchive'
    )

    args = parser.parse_args()

    return args


def unarchive_files(sample_vcf_dict) -> None:
    """
    Unarchive any original VCF files which are not 'live'

    Parameters
    ----------
    sample_vcf_dict : dict
        dict containing each sample as key, with a list of dictionaries
        as the value containing info all of the original VCFs found
        for that sample
    """
    archived = {
        k: v for k, v in sample_vcf_dict.items()
        if v.get('archive') != 'live' and v.get('id')
    }

    for idx, file in enumerate(archived.values()):
        print(f"Checking file {idx + 1}/{len(archived.values())}")
        file_id = file.get('id')
        proj_id = file.get('project')
        file_object = dx.DXFile(file_id, project=proj_id)
        file_object.unarchive()
        time.sleep(5)


def main():
    args = parse_args()
    # Open the JSON containing the file IDs found for each sample
    original_vcf_ids = read_in_json(args.vcf_json)
    # Unarchive the files
    unarchive_files(original_vcf_ids)


if __name__ == "__main__":
    main()
