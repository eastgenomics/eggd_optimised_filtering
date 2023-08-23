"""
This script creates a CSV for each sample describing the GM number, X number,
report outcome and the number of variants which would be returned to scientists
with current routine filtering methods
"""
import pandas as pd

from pathlib import Path

import utils.file_utils as file_utils


ROOT_DIR = Path(__file__).absolute().parents[1]


def add_in_x_number_key(sample_dict):
    """
    For each sample, get the X number and add as a new key

    Parameters
    ----------
    sample_dict : dict
        dictionary containing each sample as key with info about the
        sample as the values

    Returns
    -------
    sample_dict : dict
        dictionary containing each sample as key with extra attribute
        x_number
    """
    for gm_number, vcf_info in sample_dict.items():
        file_name = vcf_info.get('filename')
        if file_name:
            if 'GM' in file_name.upper():
                x_number = file_name.split('-')[0]
            else:
                x_number = file_name.split('_')[0]
            sample_dict[gm_number]['x_number'] = x_number
        else:
            sample_dict[gm_number]['x_number'] = ''

    return sample_dict


def get_number_of_routine_variants():
    """
    Read in file and create list of each sample with the number of
    routine variants which would be returned

    Returns
    -------
    _type_
        _description_
    """
    with open(
        ROOT_DIR.joinpath('resources', 'total_included_variants.txt')
    ) as f:
        lines = f.readlines()

    stripped = [line.replace("\n"," ").strip() for line in lines]
    # Group the two lines together which are the file name and the number of
    # variants
    grouped = [stripped[i:i+2] for i in range(0, len(stripped), 2)]

    # Split the filename to get X number with number of included variants
    x_number_list = []
    for nested_list in grouped:
        file_name = nested_list[0]
        if 'GM' in file_name.upper():
            x_number = file_name.split('-')[0]
        else:
            x_number = file_name.split('_')[0]
        updated_list = [x_number, nested_list[1]]
        x_number_list.append(updated_list)

    return x_number_list


def create_merged_df(x_number_list, sample_dict):
    """
    Create a merged dataframe containing GM number, X number, report outcome
    and number of routine variants

    Parameters
    ----------
    x_number_list : list
        list with X number and number of routine filtered variants
    """
    var_number_df = pd.DataFrame(
        x_number_list, columns=['x_number', 'routine_filter_variants']
    )

    # Create df from the big JSON with each GM number and info about the VCF
    # found
    gm_df = pd.DataFrame.from_dict(sample_dict, orient='index').reset_index()
    # Subset to only have the GM number, X number and the report outcome
    subset_gm = gm_df[["index", "x_number", "report_outcome"]]

    # Merge to get one df with the GM number, X number, report outcome and no
    # of variants left after routine filtering
    merged_df = pd.merge(subset_gm, var_number_df, on='x_number', how='left')

    # Remove NAs (the 4 samples with no VCF found)
    merged_df = merged_df[merged_df['routine_filter_variants'].notna()]
    merged_df = merged_df.rename(columns={"index": "gm_number"})

    # Write to CSV
    #merged_df.to_csv('variants_with_routine_filters.csv', index=False)

    return merged_df


def get_mean_number_of_variants(a_dataframe):
    """
    Get the mean number of variants returned

    Parameters
    ----------
    a_dataframe : pd.DataFrame
        a Pandas dataframe with the routine_filtered_variants column

    Returns
    -------
    mean_no_variants: float
        the mean of the routine_filter_variants column
    """
    a_dataframe['routine_filter_variants'] = a_dataframe[
        'routine_filter_variants'
    ].astype(int)

    mean_no_variants = a_dataframe['routine_filter_variants'].mean()

    return mean_no_variants


def print_mean_number_of_variants(merged_df):
    """
    Print out the mean number of routine filtered variants

    Parameters
    ----------
    merged_df : pd.DataFrame
        dataframe with a row for each sample with relevant in each column
    """
    positive_cases = merged_df.loc[
        merged_df.report_outcome.isin(['MUT','UVAR'])
    ]

    nmd_cases = merged_df.loc[merged_df.report_outcome == 'NMD']

    print(
        f"A mean of {get_mean_number_of_variants(positive_cases)} variants were returned to scientists"
        " in positive cases"
    )

    print(
        f"A mean of {get_mean_number_of_variants(nmd_cases)} variants were returned to scientists"
        " in negative cases"
    )


def main():
    # Read in JSON containing info about each sample
    sample_dict = file_utils.read_in_json_from_local_file(
        'resources', 'sample_file_IDs_outcome.json'
    )
    sample_dict = add_in_x_number_key(sample_dict)
    x_number_list = get_number_of_routine_variants()
    merged_df = create_merged_df(x_number_list, sample_dict)
    print_mean_number_of_variants(merged_df)
