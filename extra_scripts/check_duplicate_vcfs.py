"""
This script compares the variants where duplicate VCFs exist so we can
decide how to choose one VCF per sample
"""
from pathlib import Path
import subprocess

from utils import file_utils


# Get the path to the main directory
ROOT_DIR = Path(__file__).absolute().parents[1]


# Open file that we can append results of diff to
# For each sample, if there are 2 files found

def write_out_diff_comparison(original_vcf_ids):
    #TODO change this to read in files with dxpy rather than using
    # subprocess
    """
    Where duplicate VCF files are found for one sample, diff the
    files (excluding the header) to check the same variants are present

    Parameters
    ----------
    original_vcf_IDs : dict
        dict of each GM number as key and a list of the VCFs objs found as val
    """
    with open(
        ROOT_DIR.joinpath('resources', 'outcome_of_comparisons.txt'), 'a'
    ) as output_file:
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
    file_IDs_to_check = file_utils.read_in_json_from_local_file(
        "resources", "sample_VCF_IDs.json"
    )
    write_out_diff_comparison(file_IDs_to_check)


if __name__ == '__main__':
    main()
