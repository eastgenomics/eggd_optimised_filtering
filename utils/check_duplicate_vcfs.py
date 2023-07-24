import dxpy as dx
import json
import subprocess

from pathlib import Path


# Get the path to the main directory
ROOT_DIR = Path(__file__).absolute().parents[1]

with open(
    ROOT_DIR.joinpath("utils", "file_IDs_modified.json"), "r", encoding='utf8'
) as json_file:
    original_vcf_IDs = json.load(json_file)


# Open file that we can append results of diff to
# For each sample, if there are 2 files found

with open('outcome_of_comparisons_modified.txt', 'a') as output_file:
    # For each sample, if there are 2 files found
    for sample, files in original_vcf_IDs.items():
        if len(files) == 2:
            projects = [file['project'] for file in files]
            # # Only compare files where any of them aren't in these projects
            # if not 'project-GF1fV6Q4VfVz5qF9PQfZk3g8' in projects:
            #     if not 'project-G618GZQ4F21G2vg4FjKBGP2q' in projects:
            # List the IDs of the two files
            file_IDs = [file['id'] for file in files]
            base_command = "diff "
            # Build up the diff command to compare the two files (excluding the header lines) using the file IDs
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
            # Append the output for each sample to a file
            output_file.write(f"{sample}\n{file_IDs}\n{output}\n\n")
