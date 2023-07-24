"""
Get only one VCF file ID per sample, output to JSON and copy those files
into the testing project
"""
import dxpy as dx
import json
import pandas as pd
import query_original_vcf as qov

from collections import defaultdict
from pathlib import Path


# Get path of main dir
ROOT_DIR = Path(__file__).absolute().parents[1]

# List samples which failed QC on first instance
SAMPLES_WITH_QC_FAILS = ['GMXXXXXXX', 'GMXXXXXXX', 'GMXXXXXX']
# Define the project where the above samples failed QC on the first run
PROJECT_WHERE_FAILED = 'project-G5YPJZj4574zvxkVJqQPvvvY'


def read_in_json_file(folder, file_name):
    """
    Read in a JSON file

    Returns
    -------
    json_dict : dict
        the JSON converted to a Python dictionary
    """
    with open(
        ROOT_DIR.joinpath(folder, file_name), "r", encoding='utf8'
    ) as json_file:
        json_dict = json.load(json_file)

    return json_dict


def remove_failed_sample_vcfs(vcf_dict):
    """
    Remove three samples which failed QC the first time they were tested

    Parameters
    ----------
    vcf_dict : dict
        dict with GM number as key and list of dicts containing info on
        VCFs found for that sample

    Returns
    -------
    updated_dict : dict
        the dict as above but with the specific VCFs for the three samples
        removed if they are found in the project where they failed QC
    """
    updated_dict = defaultdict(list)

    for sample, files in vcf_dict.items():
        # If the sample is one of the three that failed
        if sample in SAMPLES_WITH_QC_FAILS:
            # Remove VCF file in specific project because
            # this is the instance of the sample that failed QC
            updated_files = [
                x for x in files if not (
                    PROJECT_WHERE_FAILED == x.get('project')
                )
            ]
            updated_dict[sample] = updated_files
        # If not sample that failed QC just take the files as the value
        else:
            updated_dict[sample] = files

    return updated_dict


def get_only_latest_vcf(modified_vcf_dict):
    """
    In cases where there are multiple VCFs found, take one single VCF
    which is the most recent one. The VCFs were compared and found to be
    identical apart from the header, so this is a way of taking a single VCF
    per sample

    Parameters
    ----------
    modified_vcf_dict : dict
        dict where the failed sample VCFs were removed

    Returns
    -------
    final_dict : dict
        dict with only 1 VCF at most per sample instead of duplicates
    """
    final_dict = defaultdict(list)

    for sample, files in modified_vcf_dict.items():
        # If more than one VCF found
        if len(files) > 1:
            # Get just the dict for the VCF file which was created latest
            recent_vcf = [max(files, key=lambda x:x['created'])]
            # Update the value to be this instead
            final_dict[sample] = recent_vcf[0]
        elif len(files) == 1:
            final_dict[sample] = files[0]
        else:
            # Otherwise, if just one VCF or no VCF, keep this as the value
            final_dict[sample] = {}

    return final_dict


def add_testing_outcome(vcf_dict, obesity_df):
    """
    Add in the testing outcome for the sample

    Parameters
    ----------
    vcf_dict : dict
        dict with sample as key and value is a dict of the one VCF for
        that sample
    obesity_df : pd.DataFrame
        pandas df of each sample and testing outcome and (if relevant)
        the variant(s) found

    Returns
    -------
    vcf_dict : dict
        final dict with GM number as key and dict with single VCF and sample
        info as value
    """
    for sample, _ in vcf_dict.items():
        # Get the report outcome from the dataframe based on the GM no
        report_outcome = obesity_df.loc[
            obesity_df['LABNO'] == sample, 'REPORT_TYPE'
        ].item()
        vcf_dict[sample]['report_outcome'] = report_outcome

    return vcf_dict


def write_out_final_json(folder, file_name, sample_vcf_dict):
    """
    Write out the final dictionary to a JSON
    Parameters
    ----------
    folder : str
        name of folder to write the file to
    file_name : str
        name of the output JSON file
    sample_vcf_dict : dict
        final dict with GM number as key and dict with single VCF and sample
        info as value
    """
    with open(ROOT_DIR.joinpath(folder, file_name), 'w') as fp:
        json.dump(sample_vcf_dict, fp, indent=4)


def copy_files_to_testing_project(vcf_dict, test_project_id):
    """
    Copy the one file per sample to my DNAnexus testing project

    Parameters
    ----------
    vcf_dict : dict
        final dict with GM number as key and dict with single VCF and sample
        info as value
    test_project_id : str
        project ID for the testing project
    """
    for file in vcf_dict.values():
        if file:
            file_id = file['id']
            proj_id = file['project']
            file_object = dx.DXFile(file_id, project=proj_id)
            file_object.clone(test_project_id, folder='/vcfs')


def main():
    obesity_df = qov.read_sample_spreadsheet(
        'resources', '230609_obesity_no_dups.csv'
    )
    sample_vcf_dict = read_in_json_file('resources', 'sample_VCF_IDs.json')
    sample_vcf_updated = remove_failed_sample_vcfs(sample_vcf_dict)
    vcf_dict = get_only_latest_vcf(sample_vcf_updated)
    final_vcf_dict = add_testing_outcome(vcf_dict, obesity_df)
    write_out_final_json(
        'resources', 'sample_file_IDs_outcome.json', final_vcf_dict
    )
    #copy_files_to_testing_project(
    #     final_vcf_dict, 'project-GVqVPk04p65vjXb2kj6FqFKf'
    # )

if __name__ == '__main__':
    main()
