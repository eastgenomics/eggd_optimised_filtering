"""
General functions which are shared between files
"""
import dxpy as dx
import json
import pandas as pd

from pathlib import Path


ROOT_DIR = Path(__file__).absolute().parents[1]


def read_in_json_from_local_file(folder, file_name):
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


def read_in_json_from_dnanexus(file_id):
    """
    Read in a JSON file from DNAnexus

    Returns
    -------
    json_dict : dict
        the JSON converted to a Python dictionary
    """
    proj_id, file_id = file_id.split(":")

    with dx.open_dxfile(file_id, project=proj_id) as json_file:
        json_contents = json.load(json_file)

    return json_contents


def write_out_json(folder, file_name, dict_to_write_out) -> None:
    """
    Write out a dictionary to a JSON file

    Parameters
    ----------
    folder : str
        name of folder to write the file to
    file_name : str
        name of the output JSON file
    dict_to_write_out : dict
        dictionary to write to a JSON
    """
    with open(
        ROOT_DIR.joinpath(folder, file_name), 'w', encoding='utf8'
    ) as fp:
        json.dump(dict_to_write_out, fp, indent=4)


def read_in_csv(folder, file_name):
    """
    Read in spreadsheet to pandas dataframe

    Parameters
    ----------
    folder : str
        name of folder spreadsheet is saved in
    file_name : _type_
        name of spreadsheet

    Returns
    -------
    pd.DataFrame
        CSV converted to pandas dataframe table
    """

    return pd.read_csv(ROOT_DIR.joinpath(folder, file_name))
