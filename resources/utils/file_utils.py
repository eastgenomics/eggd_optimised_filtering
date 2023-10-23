"""
General functions which are shared between files
"""
import json
import pandas as pd

from pathlib import Path


def read_in_json(file_name):
    """
    Read in a JSON file to a dict

    Parameters
    ----------
    file_name : str
        name of JSON file to read in
    Returns
    -------
    json_dict : dict
        the JSON converted to a Python dictionary
    """
    with open(file_name, "r", encoding='utf8') as json_file:
        json_dict = json.load(json_file)

    return json_dict


def write_out_json(file_name, dict_to_write_out) -> None:
    """
    Write out a dictionary to a JSON file

    Parameters
    ----------
    file_name : str
        name of the output JSON file
    dict_to_write_out : dict
        dictionary to write to a JSON
    """
    with open(file_name, 'w', encoding='utf8') as fp:
        json.dump(dict_to_write_out, fp, indent=4)


def read_in_csv(file_name):
    """
    Read in spreadsheet to pandas dataframe

    Parameters
    ----------
    file_name : str
        name of file to be read in

    Returns
    -------
    pd.DataFrame
        CSV converted to pandas dataframe table
    """

    return pd.read_csv(file_name)


def unescape_bcftools_command(bcftools_filter_command):
    """
    Removes extra backslashes from bcftools filter command because
    it has to be escaped in the JSON

    Parameters
    ----------
    bcftools_filter_command : str
        full escaped bcftools filter command directly from JSON

    Returns
    -------
    unescaped_command : str
        full unescaped bcftools filter command
    """

    return json.loads(json.dumps(bcftools_filter_command))
