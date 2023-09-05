"""
This script adds in filtered variants (counts and actual variants from
routine or optimised filtering) to a base excel (for comparisons between
the variants)
"""
import argparse
import pandas as pd


def parse_args() -> argparse.Namespace:
    """
    Parse the command line arguments inputs given

    Returns
    -------
    args : Namespace
        Namespace of passed command line argument inputs
    """

    parser = argparse.ArgumentParser(
        description='Information necessary for adding routine filter variants'
    )

    parser.add_argument(
        '-i',
        '--input_excel',
        type=str,
        required=True,
        help='Base Excel file with each sample and the variant(s) reported'
    )

    parser.add_argument(
        '-v',
        '--variants',
        type=str,
        required=True,
        help=(
            'txt file of each case and the variants returned in each'
            ' case'
        )
    )

    parser.add_argument(
        '-o',
        '--output_excel',
        type=str,
        required=True,
        help='Name of Excel file with each case plus variants added'
    )

    parser.add_argument(
        '-t',
        '--variant_type',
        type=str,
        required=True,
        choices=['routine', 'optimised'],
        help='Whether filter variants being added are routine or optimised'
    )

    args = parser.parse_args()

    return args


def read_in_txt_file(path_to_txt_file):
    """
    Read in a .txt file as a list of lines

    Parameters
    ----------
    path_to_txt_file : str
        Path to the txt file to read in

    Returns
    -------
    lines : list
        list containing each line from the file
    """
    with open(path_to_txt_file, encoding='utf8') as my_file:
        lines = my_file.readlines()

    return lines


def group_lines_by_case(lines_of_txt_file):
    """
    Create nested lists from one big list, with one list per case

    Parameters
    ----------
    lines_of_txt_file : list
        list containing each line of the txt file

    Returns
    -------
    grouped_list : list
        list of lists, with each list representing each case with all the
        related variants
    """
    # Split each case by \n char, creating a list with variants for each case
    # and adding to main list
    split_value = '\n'
    grouped_list = []
    temp_list = []

    for line in lines_of_txt_file[:-1]:
        if line == split_value:
            grouped_list.append(temp_list)
            temp_list = []
        else:
            temp_list.append(line)
    grouped_list.append(temp_list)

    return grouped_list


def change_filename_to_x_number(grouped_list):
    """
    Update the filename in the list to be the case's X number

    Parameters
    ----------
    grouped_list : list
        list of lists, with each list representing each case with all the
        related variants

    Returns
    -------
    grouped_x_number_list : list
        list of lists but each list contains sample X number, a string of
        all the variants prioritised and the number of variants prioritised
    """
    grouped_x_number_list = []
    for nested_list in grouped_list:
        file_name = nested_list[0]
        variants = nested_list[1:]
        stripped_variants = [
            variant.replace('\t', ' ') for variant in variants
        ]
        stripped_variants = [
            variant.replace('\n', '') for variant in stripped_variants
        ]
        joined_variants = '\n'.join(str(i) for i in stripped_variants)

        # Replace filename in list with sample X number
        if 'GM' in file_name.upper():
            x_number = file_name.split('-')[0]
        else:
            x_number = file_name.split('_')[0]

        updated_list = [x_number, joined_variants, len(stripped_variants)]
        grouped_x_number_list.append(updated_list)

    return grouped_x_number_list


def create_variant_df(grouped_x_number_list, variant_type):
    """
    Create a dataframe from list of lists with variants for each case

    Parameters
    ----------
    grouped_x_number_list : list
        list of lists but each list contains sample X number, a string of
        all the variants prioritised and the number of variants prioritised
    variant_type : str
        'routine' or 'optimised' - whether filtered variants are from
        routine or optimised filtering

    Returns
    -------
    variant_df : pd.DataFrame
        Dataframe with sample X number, variants prioritised and the number
        of variants prioritised
    """
    variant_df = pd.DataFrame(
        grouped_x_number_list, columns=[
            'x_number', f'{variant_type}_variants', f'{variant_type}_count'
        ]
    )

    return variant_df


def create_merged_df(base_df, variant_df):
    """
    Merge a base df with variants prioritised for that case

    Parameters
    ----------
    base_df : pd.DataFrame
        base dataframe to add prioritised variants to
    variant_df : pd.DataFrame
        dataframe of prioritised variants to be added to a base df

    Returns
    -------
    merged_df : pd.DataFrame
        a base df merged with variants prioritised for that case
    """

    merged_df = pd.merge(base_df, variant_df, on='x_number', how='left')

    return merged_df


def main():
    args = parse_args()
    base_df = pd.read_excel(args.input_excel)
    filter_variants = read_in_txt_file(args.variants)
    grouped_variants = group_lines_by_case(filter_variants)
    x_number_filter_vars = change_filename_to_x_number(grouped_variants)
    variant_df = create_variant_df(x_number_filter_vars, args.variant_type)
    merged_base_routine_df = create_merged_df(base_df, variant_df)
    merged_base_routine_df.to_excel(args.output_excel, index=False)


if __name__ == "__main__":
    main()
