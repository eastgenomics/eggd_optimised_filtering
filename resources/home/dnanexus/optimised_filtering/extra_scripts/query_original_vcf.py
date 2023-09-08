"""
Script to find the original VCF(s) outputted from Sentieon based on a
sample's GM number (or X number)
"""
import argparse
import dxpy as dx
import os
import re
import sys
import warnings

from collections import defaultdict
from pathlib import Path


sys.path.append(os.path.abspath(
    os.path.join(os.path.realpath(__file__), '../../')
))

from utils.file_utils import read_in_csv, write_out_json

warnings.simplefilter(action='ignore', category=FutureWarning)


def parse_args() -> argparse.Namespace:
    """
    Parse the command line arguments inputs given

    Returns
    -------
    args : Namespace
        Namespace of passed command line argument inputs
    """

    parser = argparse.ArgumentParser(
        description='Information necessary to query original VCFs'
    )

    parser.add_argument(
        '-s',
        '--shire_output',
        type=str,
        required=True,
        help='Path to formatted CSV output of reported cases from Shire'
    )

    parser.add_argument(
        '-m',
        '--x_gm_mapping',
        type=str,
        required=True,
        help='CSV file with two columns containing mapping of X to GM no'
    )

    parser.add_argument(
        '-o',
        '--output_json',
        type=str,
        required=True,
        help='Name of output JSON file'
    )

    parser.add_argument(
        '-p',
        '--search_projects',
        action='store',
        type=str,
        required=True,
        nargs='+',
        help='Names of 002 project suffixes to search on DNAnexus'
    )

    args = parser.parse_args()

    return args


def find_002_projects(project_ending):
    """
    Gets all the 002 diagnostic DNAnexus project IDs for projects ending with
    the relevant assay type
    Parameters
    ----------
    project_ending : str
        e.g. 'TWE'
    Returns
    -------
    project_ids : list
        list of project IDs
    """
    # Search projs with name beginning with 002, ending with assay type
    # Return only project ID field
    assay_dx_response = list(dx.find_projects(
        level='VIEW',
        name=f"002*{project_ending}",
        name_mode="glob",
        describe={
            'fields': {
                'id': True
            }
        }
    ))

    # Get a list of just project IDs to search through
    project_ids = [project['id'] for project in assay_dx_response]

    return project_ids


def get_vcfs_in_proj(project_id):
    """
    Find the raw VCF files (outputted from Sentieon app) in a specific project,
    returning the file ID, file name, its archival state and created
    time (epoch format)
    Parameters
    ----------
    project_id : str
        i.e. 'project-XXXXXX'
    Returns
    -------
    vcf_info : list
        list of dictionaries, each dict holding information about a VCF file
    """
    vcf_info = list(
        dx.find_data_objects(
            project=project_id,
            name="*recalibrated_Haplotyper.vcf.gz",
            name_mode='glob',
            classname='file',
            describe={
                'fields': {
                    'name': True,
                    'archivalState': True,
                    'created': True
                }
            }
        )
    )

    return vcf_info


def get_all_vcfs_in_projects(projects_002):
    """
    Given a list of project IDs, find all the raw VCFs in the projects

    Parameters
    ----------
    projects_002 : list
        list of project IDs to search for VCFs in

    Returns
    -------
    vcfs_found :  list
        list of dicts containing info about all the VCF files found
    """
    vcfs_found = []
    # For each 002 project, find the VCFs of interest in the project
    # Add the VCFs found to a list
    for project in projects_002:
        vcf_file_info = get_vcfs_in_proj(project)
        vcfs_found.extend(vcf_file_info)

    return vcfs_found


def modify_sample_names(data_frame):
    """
    Remove '.' from GM sample names in the LABNO column

    Parameters
    ----------
    data_frame: pd.DataFrame
        dataframe to strip '.' from the LABNO column

    Returns
    -------
    data_frame : pd.DataFrame
        dataframe where GM number has '.' removed
    """

    data_frame['LABNO'] = data_frame['LABNO'].str.replace(".", "")

    return data_frame


def get_list_of_gm_numbers(sample_df):
    """
    Read in the GM numbers of the samples of interest

    Parameters
    ----------
    sample_df : pd.DataFrame
        dataframe containing GM number, outcome and, if pos, variant(s)
        of interest

    Returns
    -------
    gm_numbers :  list
        list of all the sample GM numbers
    """

    gm_numbers = sample_df['LABNO'].tolist()

    return gm_numbers


def find_files_by_gm_number(list_of_vcf_file_objs, gm_numbers):
    #TODO compare two dicts by key (make list_of_vcf_file_objs filename as key) # rather than looping over twice
    """
    Create a dictionary of GM number key and value is list all of VCF
    objects found for that sample

    Parameters
    ----------
    list_of_vcf_file_objs : list
        list of dicts containing VCF objects
    gm_numbers : list
        list of GM numbers (strings)

    Returns
    -------
    relevant_vcf_info : dict
        dict with GM number as key and list of all VCFs found for that sample
        as value
    """
    relevant_vcf_info = defaultdict(list)

    # For each VCF file, get the file name
    for vcf_file in list_of_vcf_file_objs:
        file_name = vcf_file['describe']['name']
        # Compare the file name to GM number to see if a match
        for gm_number in gm_numbers:
            if re.search(gm_number, file_name, re.IGNORECASE):
                # If match, collect info about the file
                project = vcf_file['project']
                file_id = vcf_file['id']
                archive_state = vcf_file['describe']['archivalState']
                created = vcf_file['describe']['created']

                # Add the info for the file to dict with GM number as key
                relevant_vcf_info[gm_number].extend([
                    {
                        'project': project,
                        'id': file_id,
                        'filename': file_name,
                        'archive': archive_state,
                        'created': created
                    }
                ])

    return relevant_vcf_info


def add_in_samples_not_found(sample_vcf_dict, gm_numbers):
    #TODO Initialise dict from gm_numbers list with empty list as value above
    # so this function isn't needed
    """
    Adds empty list for samples where no VCF was found

    Parameters
    ----------
    sample_vcf_dict : dict
        dict with GM number as key and list of all VCFs found for that sample
        as value
    gm_numbers : list
        list of GM numbers (strings)

    Returns
    -------
    sample_vcf_dict : dict
        dict with GM number as key and list of all VCFs found for that sample
        as value. Samples with no VCF(s) found have value as empty list
    """
    for gm_number in gm_numbers:
        # If the GM number doesn't exist in the dict (no VCFs found)
        # then add empty string as value
        if gm_number not in sample_vcf_dict:
            sample_vcf_dict[gm_number] = []

    return sample_vcf_dict


def create_x_gm_mapping(x_gm_dataframe):
    """
    Creates a mapping of GM number to X number, as VCF file names
    did not have GM number in the name in the past but did have X number

    Parameters
    ----------
    x_gm_dataframe : pd.DataFrame
        dataframe containing GM number and corresponding X number

    Returns
    -------
    x_gm_dict : dict
        dictionary of GM number key with corresponding X number as value
    """
    # Remove the '.' from each GM number
    x_gm_dataframe['LabNumber'] = x_gm_dataframe['LabNumber'].str.replace(
        '.', ''
    )
    # Create a dict of GM number (with GM in capitals) vs X number
    x_gm_dict = dict(
        zip(
            x_gm_dataframe.LabNumber.str.upper(),
            x_gm_dataframe.ExomeNumber
        )
    )
    return x_gm_dict


def match_file_to_x_number(list_of_vcf_objects, x_number, gm_number):
    """
    Take an X number and find any files in the list of VCFs which
    have the X number in their file name

    Parameters
    ----------
    list_of_vcf_objects : list
        list of dicts containing VCF objects
    x_number : str
        the sample's X number
    gm_number : str
        the sample's GM number

    Returns
    -------
    relevant_vcf_info : list
        list of dicts containing VCF objects matching that sample
    """
    relevant_vcf_info = []
    # Loop over each VCF file object we found earlier
    for vcf_file in list_of_vcf_objects:
        file_name = vcf_file['describe']['name']
        # If the X number matches a file name
        # get the relevant file info
        if re.search(x_number, file_name, re.IGNORECASE):
            project = vcf_file['project']
            file_id = vcf_file['id']
            archive_state = vcf_file['describe']['archivalState']
            created = vcf_file['describe']['created']
            # Add this info as a dict to a list (as multiple files
            # might be found matching the X number)
            vcf_info = {
                'project': project,
                'name': gm_number,
                'filename': file_name,
                'id': file_id,
                'archive_state': archive_state,
                'created': created
            }
            relevant_vcf_info.append(vcf_info)

    return relevant_vcf_info


def find_vcf_based_on_x_number(sample_vcf_dict, x_gm_dict, all_vcf_objects):
    """
    For a given sample's GM number, get the corresponding X number
    and search through the list of VCFs to find any matching files

    Parameters
    ----------
    sample_vcf_dict : dict
        dict with GM number as key and VCF objects as value (list)
    x_gm_dict : dict
        dict of GM number to X number mapping
    all_vcf_objects : list
        list of dicts containing VCF objects

    Returns
    -------
    relevant_vcf_info : list
        list of dicts containing VCF objects matching that sample
    """
    # For each entry in the overall dictionary
    for sample, files in sample_vcf_dict.items():
        # If no files were found
        if not files:
            # Find the X number by the GM to search with that instead
            x_number = x_gm_dict.get(sample)
            # If an X number is found for that GM number
            if x_number:
                # Search for file matches based on X number instead
                file_match = match_file_to_x_number(
                    all_vcf_objects, x_number, sample
                )
                # If theres a match, append the file info to a list
                if file_match:
                    formatted_file_info = [
                        {
                            'project': x['project'],
                            'id': x['id'],
                            'filename': x['filename'],
                            'archive': x['archive_state'],
                            'created': x['created']
                        } for x in file_match
                    ]
                    # Add the info for the files found to the dict
                    sample_vcf_dict[sample].extend(formatted_file_info)
                else:
                    print("No files found:", sample)
            if not x_number:
                print("No X Number found:", sample)

    return sample_vcf_dict


def check_files_found(final_sample_vcf_dict) -> None:

    samples_where_no_vcf_found = [
        sample for sample, files in final_sample_vcf_dict.items() if not files
    ]
    print(
        f"There were {len(samples_where_no_vcf_found)} samples "
        "with no VCF found"
    )

    duplicate_vcf_samples = [
        sample for sample, files in final_sample_vcf_dict.items() if len(files) > 1
    ]
    print(
        f"There were {len(duplicate_vcf_samples)} samples "
        "with >1 VCF found"
    )

    single_vcf_found = [
        sample for sample, files in final_sample_vcf_dict.items() if len(files) == 1
    ]
    print(
        f"There were {len(single_vcf_found)} samples "
        "where one VCF was found"
    )

    more_than_two_found = [
        sample for sample, files in final_sample_vcf_dict.items() if len(files) > 2
    ]
    print(
        f"There were {len(more_than_two_found)} samples "
        "where more than 2 VCFs were found"
    )


def main():
    args = parse_args()
    all_002_project_ids = []
    for project_ending in args.search_projects:
        all_002_project_ids += find_002_projects(project_ending)
    vcfs_in_projs = get_all_vcfs_in_projects(all_002_project_ids)
    obesity_df = read_in_csv(args.shire_output)
    obesity_df = modify_sample_names(obesity_df)
    gm_nos = get_list_of_gm_numbers(obesity_df)
    my_vcfs = find_files_by_gm_number(vcfs_in_projs, gm_nos)
    all_sample_vcfs_found = add_in_samples_not_found(my_vcfs, gm_nos)
    x_gm_mapping = read_in_csv(args.x_gm_mapping)
    x_gm_dict = create_x_gm_mapping(x_gm_mapping)
    final_sample_vcf_dict = find_vcf_based_on_x_number(
        all_sample_vcfs_found, x_gm_dict, vcfs_in_projs
    )
    write_out_json(args.output_json, final_sample_vcf_dict)
    check_files_found(final_sample_vcf_dict)

if __name__ == '__main__':
    main()
