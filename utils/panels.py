import json

from collections import defaultdict
from panelapp import api, Panelapp, queries
from pathlib import Path

ROOT_DIR = Path(__file__).absolute().parents[1]


def read_in_mappings():
    """
    Read in the inheritance terms and their map to simpler alternatives

    Returns
    -------
    mappings : dict
        dictionary with each MOI term and its simpler alternative
    """
    # Read in inheritance mappings from JSON file
    with open(
        ROOT_DIR.joinpath(f"{ROOT_DIR}/resources/inheritance_mapping.json"), 'r', encoding='utf8'
    ) as json_file:
        content = json.load(json_file)

    mappings = content.get('mappings')

    return mappings

def get_panel_data(panel_id):
    """
    Get information on the panel by its ID

    Parameters
    ----------
    panel_id : int
        ID of the panel of interest

    Returns
    -------
    genes : 
        _description_
    """
    panel_data = Panelapp.Panel(panel_id).get_data()

    return panel_data

def get_gene_region_info(panel_data):
    genes = panel_data.get('genes')
    regions = panel_data.get('regions')
    print(type(genes), type(regions))

    return genes, regions


def format_panel_info_as_dict(genes, regions):

    panel_dict = defaultdict(dict)
    if genes:
        for gene in genes:
            gene_symbol = gene.get('gene_data').get('gene_symbol')
            moi = gene.get('mode_of_inheritance')
            conf_level = int(gene.get('confidence_level'))
            if conf_level >= 3:
                panel_dict[gene_symbol]['mode_of_inheritance'] = moi
    if regions:
        for region in regions:
            region_name = region.get('entity_name')
            conf_level = int(region.get('confidence_level'))
            moi = region.get('mode_of_inheritance')
            if conf_level >=3:
                panel_dict[region_name]['mode_of_inheritance'] = moi

    return panel_dict


def map_moi_to_simpler_terms(panel_dict,mappings):
    new_dict = defaultdict(dict)
    for gene, values in panel_dict.items():
        inheritance = values['mode_of_inheritance']
        mapped_value = mappings.get(inheritance)
        new_dict[gene]['mode_of_inheritance'] = mapped_value
    return new_dict

def get_formatted_dict(panel_id):
    mappings = read_in_mappings()
    panel_data = get_panel_data(panel_id)
    genes, regions = get_gene_region_info(panel_data)
    panel_as_dict = format_panel_info_as_dict(genes, regions)
    final_dict = map_moi_to_simpler_terms(panel_as_dict, mappings)
    return final_dict
