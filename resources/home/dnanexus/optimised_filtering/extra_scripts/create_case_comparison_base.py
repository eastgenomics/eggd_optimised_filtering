"""
This script creates a simple dataframe which acts as a base for comparing
routine and optimised filter variants
"""
import argparse
import os
import pandas as pd
import sys

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
        description='Files needed to create base df for comparisons'
    )

    parser.add_argument(
        '-i',
        '--input_json',
        type=str,
        required=True,
        help='JSON file which has each sample with one VCF found'
    )

    parser.add_argument(
        '-o',
        '--output_excel',
        type=str,
        required=True,
        help='Name of output Excel'
    )

    args = parser.parse_args()

    return args


def create_base_df(vcf_json):
    """
    Create a df to act as a base for comparisons between routine and optimised
    variants

    Parameters
    ----------
    vcf_json : dict
        dictionary containing each sample and the original VCF found

    Returns
    -------
    case_df : pd.DataFrame
        dataframe with columns gm_number, x_number and report_outcome
    """
    case_df = pd.DataFrame.from_dict(vcf_json, orient='index').reset_index()
    # Subset to only columns of interest
    case_df = case_df[["index", "x_number", "report_outcome"]]

    # Remove samples where no VCF was found
    case_df = case_df[case_df["x_number"].notna()]

    # Rename index column to gm_number
    case_df = case_df.rename(columns={"index": "gm_number"})

    sortbox = {'MUT': 1, 'UVAR': 2, 'NMD':3}
    case_df['sort_column'] = case_df.report_outcome.map(sortbox)

    case_df = case_df.sort_values('sort_column').drop(
        'sort_column', axis=1
    ).reset_index(drop=True)

    return case_df


def main():
    args = parse_args()
    case_dict = read_in_json(args.input_json)
    base_df = create_base_df(case_dict)
    base_df.to_excel(args.output_excel, index=False)


if __name__ == "__main__":
    main()
