import json

from collections import defaultdict
from panelapp import api, Panelapp, queries
from pathlib import Path

# Get path one directory above this file
ROOT_DIR = Path(__file__).absolute().parents[1]


def read_in_mappings():
    """
    Read in the inheritance terms and their mapping to simpler alternatives

    Returns
    -------
    mappings : dict
        dictionary with each MOI term and its simpler alternative
    """
    # Read in inheritance mappings from JSON file in /resources
    with open(
        ROOT_DIR.joinpath(f"{ROOT_DIR}/resources/inheritance_mapping.json"), 'r', encoding='utf8'
    ) as json_file:
        mappings = json.load(json_file)

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
    panel_data : dict
        dict containing information about specific panel and genes/entities on
        the panel
    Example:
    {
        'id': 123,
        'hash_id': 'XX',
        'name': 'Name of panel',
        'disease_group': 'X disorders',
        'disease_sub_group': 'X syndromes',
        'status': 'public',
        'version': 'X.X',
        'version_created': '2023-04-17T10:55:17.407924Z',
        'relevant_disorders': ['disorder'],
        'stats': {
            'number_of_genes': 44,
            'number_of_strs': 0,
            'number_of_regions': 3,
        'types': [
            {
                "name": 'Rare Disease 100K',
                'slug': 'rare-disease-100k',
                'description': 'Rare Disease 100K'
            },
            ...
        ],
        'genes': [
            {
                'gene_data': {
                    'alias': 'symbol',
                    'biotype': 'protein_coding',
                    'hgnc_id': 'HGNC:123',
                    ...
                },
                'entity_type': 'gene',
                'entity_name': 'ALMS1',
                'confidence_level': '3',
                'penetrance': 'Complete',
                'mode_of_pathogenicity': "",
                'publications': [],
                'evidence': [..],
                'phenotypes': ['X syndrome'],
                'mode_of_inheritance': "BIALLELIC, autosomal or pseudoautosomal',
                'tags': [],
                'transcript': null
            }
        ],
        'strs': [],
        'regions': [
            {
                "gene_data": null,
            "entity_type": "region",
            "entity_name": "ISCA-37404-Loss",
            "verbose_name": "15q11q13 recurrent (PWS/AS) region (BP1-BP3, Class 1) Loss",
            "confidence_level": "3",
            "penetrance": null,
            "mode_of_pathogenicity": null,
            "haploinsufficiency_score": "3",
            "triplosensitivity_score": "",
            "required_overlap_percentage": 60,
            "type_of_variants": "cnv_loss",
            "publications": [
                "22045295",
                ...
            ],
            "evidence": [
                "Expert list",
                ...
            ],
            "phenotypes": [
                "microcephaly",
                ....
            ],
            "mode_of_inheritance": "MONOALLELIC, autosomal or pseudoautosomal, imprinted status unknown",
            "chromosome": "15",
            "grch37_coordinates": null,
            "grch38_coordinates": [
                22782170,
                28134728
            }
        ]
    }
    """
    panel_data = Panelapp.Panel(panel_id).get_data()
    return panel_data

def get_gene_region_info(panel_data):
    """
    _summary_

    Parameters
    ----------
    panel_data : _type_
        _description_

    Returns
    -------
    genes : list
        list of dicts containing info about each gene
    regions : list
        list of dicts containing info about each region
    """
    genes = panel_data.get('genes')
    regions = panel_data.get('regions')

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
                panel_dict[gene_symbol]['entity_type'] = 'gene'
    if regions:
        for region in regions:
            region_name = region.get('entity_name')
            conf_level = int(region.get('confidence_level'))
            moi = region.get('mode_of_inheritance')
            if conf_level >=3:
                panel_dict[region_name]['mode_of_inheritance'] = moi
                panel_dict[region_name]['entity_type'] = 'region'

    return panel_dict


def map_moi_to_simpler_terms(panel_dict, mappings):
    new_dict = defaultdict(dict)
    # For each gene in our panel_dict, get mode of inheritance
    for gene, values in panel_dict.items():
        inheritance = values['mode_of_inheritance']
        # Get the mapped simpler value for that MOI
        mapped_value = mappings.get(inheritance)
        new_dict[gene]['mode_of_inheritance'] = mapped_value

    return new_dict

def get_formatted_dict(panel_id):
    """
    _summary_

    Parameters
    ----------
    panel_id : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """
    mappings = read_in_mappings()
    panel_data = get_panel_data(panel_id)
    genes, regions = get_gene_region_info(panel_data)
    panel_as_dict = format_panel_info_as_dict(genes, regions)
    final_dict = map_moi_to_simpler_terms(panel_as_dict, mappings)

    return final_dict
