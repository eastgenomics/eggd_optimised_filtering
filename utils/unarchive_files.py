"""
Script to unarchive the original raw VCF from Sentieon for each sample
"""

import dxpy as dx
import generate_sample_vcf_json as gsvj
import time


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
    # For each set of files for each sample
    # For each file in the set of files found
    for file in sample_vcf_dict.values():
        file_id = file.get('id')
        if file_id:
            proj_id = file['project']
            file_id = file['id']
            archive_state = file['archive']
            # If the file state isn't live, unarchive the file
            if archive_state != 'live':
                file_object = dx.DXFile(file_id, project=proj_id)
                file_object.unarchive()
                time.sleep(5)


if __name__ == '__main__':
    # Open the JSON containing the file IDs found for each sample
    original_vcf_IDs = gsvj.read_in_json_file(
        "resources", "sample_file_IDs_outcome.json"
    )
    # Unarchive the files
    unarchive_files(original_vcf_IDs)
