"""
A script to do some final comparisons between variants
"""
import argparse
import pandas as pd

pd.options.mode.chained_assignment = None


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
        help=(
            'Input Excel with each case with the number of (and actual)'
            ' reported variants and number of (and actual) routine prioritised '
            'variants and number of (and actual) optimised prioritised variants'
        )
    )

    args = parser.parse_args()

    return args


def print_pos_stats(comparison_df):
    """
    Print summary statistics for positive cases

    Parameters
    ----------
    comparison_df : pd.DataFrame
        df containing each sample and any reported variants and the prioritised
        variants from routine and optimised filtering
    """
    # Select positive cases
    positive_cases = comparison_df.loc[
        (comparison_df['report_outcome'] == 'UVAR')
        | (comparison_df['report_outcome'] == 'MUT')
    ]
    positive_cases['optimised_count'] = positive_cases[
        'optimised_count'
    ].astype(int)
    total_number_of_positive_cases = len(positive_cases)

    positive_cases_only_reported_var = len(positive_cases.loc[(positive_cases['reported_variants'] == positive_cases['optimised_variants'])])

    pos_percentage = (
        positive_cases_only_reported_var / total_number_of_positive_cases
    ) * 100

    print("Positive cases:")

    print(positive_cases['classification'].value_counts())

    print(
        "Only the reported variants are returned by optimised filtering in "
        f"{positive_cases_only_reported_var}/{total_number_of_positive_cases} "
        f"({pos_percentage:.2f}%) of positive cases"
    )

    pos_mean = positive_cases['optimised_count'].mean()

    print(
        f"A mean of {pos_mean:.2f} variants would be returned for positive "
        "cases"
    )


    positive_cases['compare_vars'] = positive_cases.apply(
        lambda x: str(x.reported_variants) in str(x.optimised_variants),
        axis=1
    )

    reported_var_not_prioritised = (~positive_cases.compare_vars).sum()

    false_neg_percentage = (
        reported_var_not_prioritised / total_number_of_positive_cases
    ) * 100

    print(
        f"In {reported_var_not_prioritised}/{total_number_of_positive_cases} "
        f"({false_neg_percentage:.2f}%) of positive cases the reported "
        "variant(s) would not be prioritised (false negative)"
    )


def print_nmd_stats(comparison_df):
    """
    Print summary statistics for NMD cases

    Parameters
    ----------
    comparison_df : pd.DataFrame
        df containing each sample and any reported variants and the prioritised
        variants from routine and optimised filtering
    """

    nmd_cases = comparison_df.loc[comparison_df['report_outcome'] == 'NMD']
    nmd_cases['optimised_count'] = nmd_cases['optimised_count'].astype(int)

    total_nmd = len(nmd_cases)

    nmd_optimised = len(nmd_cases.loc[nmd_cases['optimised_count'] == 0])
    nmd_mean = nmd_cases['optimised_count'].mean()

    nmd_percentage = (nmd_optimised / total_nmd) * 100

    print("Negative cases:")
    print(
        "Zero variants are returned by optimised filtering in "
        f"{nmd_optimised}/{total_nmd} ({nmd_percentage:.2f}%) of NMD cases"
    )


    print(
        f"A mean of {nmd_mean:.2f} variants would be returned for NMD cases"
    )

def main():
    args = parse_args()
    comparison_excel = pd.read_excel(args.input_excel)
    print_pos_stats(comparison_excel)
    print_nmd_stats(comparison_excel)


if __name__ == '__main__':
    main()
