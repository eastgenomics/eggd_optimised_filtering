"""
Script to unarchive the original raw VCF from Sentieon for each sample
"""

import dxpy as dx
import time

from utils import file_utils

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


if __name__ == '__main__':
    # Open the JSON containing the file IDs found for each sample
    original_vcf_IDs = file_utils.read_in_json_from_local_file(
        "resources", "sample_file_IDs_outcome.json"
    )
    # Unarchive the files
    unarchive_files(original_vcf_IDs)
