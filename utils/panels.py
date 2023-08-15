import dxpy as dx
import re

from collections import defaultdict
from pathlib import Path

from .file_utils import read_in_json_from_dnanexus

# Get path one directory above this file
ROOT_DIR = Path(__file__).absolute().parents[1]


def parse_genepanels(genepanels_file_id):
    """
    Parse the genepanels file to make a dict mapping clinical indication
    to PanelApp panel ID

    Parameters
    ----------
    genepanels_file_id : str
        DNAnexus file ID of the genepanels file: 'proj-XYZ:file-XYZ'

    Returns
    -------
    panel_data : dict
        dict containing each clinical ind with the panel ID as val
    """
    panel_data = {}
    proj_id, file_id = genepanels_file_id.split(":")

    with dx.open_dxfile(file_id, project=proj_id) as gp_file:
        for line in gp_file:
            panel_id, clin_ind, panel, gene = line.split('\t')
            panel_data.setdefault(clin_ind, set()).add(panel_id)

    return panel_data


def get_panel_id_from_genepanels(panel_string, genepanels_dict) -> str:
    """
    Get the panel ID from the genepanels file based on the clinical
    indication text

    Parameters
    ----------
    panel_string : str
        clinical indication e.g. 'R149.1_Severe early-onset obesity_P'
    genepanels_dict : dict
        dict containing each clinical ind with the panel ID as val

    Returns
    -------
    panel_id : str
        ID of the panel of interest, e.g. '130'
    """
    panel_set = genepanels_dict.get(panel_string)
    panel_id = [p_id for p_id in panel_set][0]

    return panel_id


def parse_panelapp_dump(panel_id, panelapp_dump):
    """
    Parse the PanelApp dump and get the dictionary for our panel

    Parameters
    ----------
    panel_id : str
        PanelApp ID of the panel
    panelapp_dump : dict
        dict of the PanelApp dump

    Returns
    -------
    _type_
        _description_
    """
    my_panel = [
        item for item in panelapp_dump if item['external_id'] == panel_id
    ][0]

    return my_panel


def format_panel_info(panel_data) -> dict:
    """
    Get the gene and region info about our panel and format into a dict
    to use for filtering

    Parameters
    ----------
    panel_data : dict
        dict of info about our panel

    Returns
    -------
    panel_dict : dict
        dict with each gene, its symbol, MOI and entity type
    """
    panel_dict = defaultdict(dict)
    genes = panel_data.get('genes')
    regions = panel_data.get('regions')

    if genes:
        for gene in genes:
            gene_symbol = gene.get('gene_symbol')
            moi = gene.get('mode_of_inheritance')
            conf_level = int(gene.get('confidence_level'))
            if conf_level >= 3:
                panel_dict[gene_symbol]['mode_of_inheritance'] = moi
                panel_dict[gene_symbol]['entity_type'] = 'gene'

    if regions:
        for region in regions:
            region_name = region.get('name')
            conf_level = int(region.get('confidence_level'))
            moi = region.get('mode_of_inheritance')
            if conf_level >=3:
                panel_dict[region_name]['mode_of_inheritance'] = moi
                panel_dict[region_name]['entity_type'] = 'region'

    return panel_dict


def simplify_MOI_terms(panel_dict) -> dict:
    """
    Simplify the mode of inheritance terms from the PanelApp panel
    because there are very slight differences in formatting or words
    between panels and genes, so one discrete mapping isn't ideal

    Parameters
    ----------
    panel_dict : dict
        default dict with each gene on the panel as key and the gene info as val

    Returns
    -------
    updated_gene_dict : dict
        default dict with each gene on the panel as key and the gene info (
        including simplified MOI) as val
    """
    updated_gene_dict = defaultdict(dict)
    for gene, moi_info in panel_dict.items():
        moi = moi_info.get('mode_of_inheritance')
        if re.search(r"^BIALLELIC", moi):
            updated_moi = 'biallelic'
        elif re.search(r"^MONOALLELIC|X-LINKED", moi):
            updated_moi = 'monoallelic'
        elif re.search(r"^BOTH", moi):
            updated_moi = 'both_monoallelic_and_biallelic'
        else:
            updated_moi = 'monoallelic'

        updated_gene_dict[gene]['mode_of_inheritance'] = updated_moi

    return updated_gene_dict


def get_formatted_dict(panel_string, genepanels_file_id, panelapp_file_id):
    """
    Main function to get a simple dictionary for each panel

    Parameters
    ----------
    panel_string : str
        The clinical indication string, including R number
    genepanels_file_id : str
        the proj:file id of the genepanels file in DNAnexus
    panelapp_file_id : str
        the proj:file id of the panelapp dump JSON in DNAnexus

    Returns
    -------
    final_panel_dict : dict
        default dict with each gene on the panel as key and the gene info as val
    """
    #test_codes = parse_R_code(panel_string)

    # Get each panel and its PanelApp ID as a dict
    genepanels_dict = parse_genepanels(genepanels_file_id)
    # Get the PanelApp ID for our panel(s) of interest
    panel_id = get_panel_id_from_genepanels(panel_string, genepanels_dict)
    # Read in the PanelApp JSON dump
    panel_dump = read_in_json_from_dnanexus(panelapp_file_id)
    # Parse the PanelApp dump to get all the info for our panel(s)
    panel_dict = parse_panelapp_dump(panel_id, panel_dump)
    # Get the gene and region info from the panel and format as dict
    panel_of_interest = format_panel_info(panel_dict)
    #final_dict = map_moi_to_simpler_terms(panel_of_interest, mappings)
    final_panel_dict = simplify_MOI_terms(panel_of_interest)

    return final_panel_dict
