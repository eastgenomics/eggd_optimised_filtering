"""
This script compares the variants where duplicate VCFs exist so we can
decide how to choose one VCF per sample
"""
import argparse
import os
import subprocess
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
        description='Information necessary to diff VCF files per sample'
    )

    parser.add_argument(
        '-i',
        '--input_json',
        type=str,
        required=True,
        help='JSON file with all of the samples and the VCFs '
    )

    parser.add_argument(
        '-o',
        '--outputfile_name',
        type=str,
        required=True,
        help='Name of output comparison txt file'
    )

    args = parser.parse_args()

    return args


# Open file that we can append results of diff to
# For each sample, if there are 2 files found

def write_out_diff_comparison(output_name, original_vcf_ids):
    #TODO change this to read in files with dxpy rather than using
    # subprocess
    """
    Where duplicate VCF files are found for one sample, diff the
    files (excluding the header) to check the same variants are present

    Parameters
    ----------
    output_name : str
        name of comparison output txt file
    original_vcf_IDs : dict
        dict of each GM number as key and a list of the VCFs objs found as val
    """
    with open(output_name, 'a') as output_file:
        # For each sample, if there are 2 files found
        for sample, files in original_vcf_ids.items():
            if len(files) == 2:
                # List the IDs of the two files
                file_IDs = [file['id'] for file in files]
                # Build up the command with diff as the starting point
                base_command = "diff "
                # Add in the file ID of each file to the diff command
                # (excluding the header lines)
                for file_id in file_IDs:
                    vcf_string = f"<(dx cat {file_id} | zcat | grep -v '^#') "
                    base_command += vcf_string
                # Run the diff of the two files using bash
                output = subprocess.run(
                    base_command,
                    shell=True,
                    executable="/bin/bash",
                    capture_output=True
                )
                # Append the output when diffing the 2 files for each sample
                # to a new file
                output_file.write(f"{sample}\n{file_IDs}\n{output}\n\n")

def main():
    args = parse_args()
    file_ids_to_check = read_in_json(args.input_json)
    write_out_diff_comparison(args.outputfile_name, file_ids_to_check)


if __name__ == '__main__':
    main()
