"""
Functions related to gene panels and obtaining information from PanelApp
"""
import re

from collections import defaultdict

from .file_utils import read_in_json


def parse_genepanels(genepanels_file):
    """
    Parse the genepanels file to make a dict mapping each clinical indication
    to PanelApp panel ID. This is so we can search the JSON PanelApp dump by
    panel ID later

    Parameters
    ----------
    genepanels_file : file
        path to genepanels file which includes clinical indications with
        panel IDs

    Returns
    -------
    panel_data : dict
        dict containing each clinical ind with the panel ID as value

    Raises
    ------
    AssertionError
        Raised if multiple panel IDs are found for one clinical indication

    Example return format:
    {
        'R136.1_Primary lymphoedema_P': {'65'},
        'R130.1_Short QT syndrome_P': {'224'}
    }
    """
    panel_data = {}

    # Open the file and read the TSV into a dict with each clinical indication
    # and a set of the panel IDs associated with it
    with open(genepanels_file, encoding="utf-8") as gp_file:
        for line in gp_file:
            clin_ind, panel, gene, panel_id = line.split('\t')
            panel_data.setdefault(clin_ind, set()).add(panel_id.strip('\n'))

    # Get any panels which have more than 1 PanelApp ID in genepanels file
    duplicate_ids = {k: sorted(v) for k, v in panel_data.items() if len(v) > 1}

    # Raise error if multiple panel IDs exist for one indication, as this
    # could indicate errors across the whole genepanels TSV
    assert not duplicate_ids, (
        f"Multiple panel IDs found for clinical indications: {duplicate_ids}"
    )

    # Print warning for any panels without a panel ID (column in TSV is empty)
    panels_with_no_id = {k for k, v in panel_data.items() if v == {''}}
    if panels_with_no_id:
        print(
            "Warning: the following panels have no panel ID found: "
            f"{panels_with_no_id}"
        )

    return panel_data


def get_panel_id_from_genepanels(panel_string, genepanels_dict):
    """
    Get the panel ID from the genepanels file based on the panel string given

    Parameters
    ----------
    panel_string : str
        clinical indication e.g. 'R149.1_Severe early-onset obesity_P'
    genepanels_dict : dict
        dict containing each clinical ind with the panel ID as value

    Returns
    -------
    panel_id : str
        ID of the panel of interest, e.g. '130'
    """
    panel_id = None
    # Get the set (i.e. panel ID) for that panel string
    panel_set = genepanels_dict.get(panel_string)
    # If panel ID exists for the panel string
    if panel_set:
        panel_ids = list(panel_set)
        assert len(panel_ids) < 2, (
            f"Multiple panel IDs found for panel string: {panel_string}"
        )
        # Get the only panel ID value in the list
        panel_id = panel_ids[0]
        if not panel_id:
            print(
                f"WARNING: Panel {panel_string} has no PanelApp ID in the "
                "given genepanels file. No MOI filtering will occur."
            )

    else:
        print(
            f"WARNING: The panel string given {panel_string} was not found in "
            "the genepanels file. This is expected if only HGNCs have been "
            "entered, but no MOI-specific filtering will be performed."
        )

    return panel_id


def transform_panelapp_dump_to_dict(panel_dump):
    """
    Converts the PanelApp dump list of dicts to a dictionary with
    PanelApp IDs as keys

    Parameters
    ----------
    panel_dump : list
        PanelApp dump as list of dicts, one for each panel

    Returns
    -------
    panel_id_dict : dict
        dict where the PanelApp ID of the panel is the key

    Raises
    ------
    AssertionError
        Raised if no panels with IDs are found in the dump
    [{
        'panel_source': 'PanelApp',
        'panel_name': 'Stickler syndrome',
        'external_id': '3',
        'panel_version': '4.0',
        'genes': [
            {
                'transcript': None,
                'hgnc_id': 'HGNC:1071',
                'confidence_level': '3',
                'mode_of_inheritance': (
                    'MONOALLELIC, autosomal or pseudoautosomal, imprinted '
                    'status unknown'
                ),
                'mode_of_pathogenicity': None,
                'penetrance': None,
                'gene_justification': 'PanelApp',
                'transcript_justification': 'PanelApp',
                'alias_symbols': None,
                'gene_symbol': 'BMP4'
            },
            ..
        }]
                            |
                            |
                            ▼
    {
        '3': {
            'panel_source': 'PanelApp',
            'panel_name': 'Stickler syndrome',
            'external_id': '3',
            'panel_version': '4.0',
            'genes': [
                {
                    'transcript': None,
                    'hgnc_id': 'HGNC:1071',
                    'confidence_level': '3',
                    'mode_of_inheritance': (
                        'MONOALLELIC, autosomal or pseudoautosomal, imprinted '
                        'status unknown'
                    ),
                    'mode_of_pathogenicity': None,
                    'penetrance': None,
                    'gene_justification': 'PanelApp',
                    'transcript_justification': 'PanelApp',
                    'alias_symbols': None,
                    'gene_symbol': 'BMP4'
                }
            }, ..
    }
    """
    # Get list of all the panel IDs for zipping
    panel_ids = []
    for panel in panel_dump:
        panel_id = panel.get('external_id')
        if panel_id:
            panel_ids.append(str(panel_id))
        else:
            print(f"Panel did not have external ID key present: {panel}")

    # Create new dict with panel ID as key
    panel_id_dict = dict(zip(panel_ids, panel_dump))

    # Raise error if the panel ID dict is empty
    if not panel_id_dict:
        raise AssertionError("No panels with IDs found in PanelApp dump")

    return panel_id_dict


def parse_panelapp_dump(panel_id, panelapp_dict):
    """
    Parse the PanelApp dump and get the dictionary for our panel

    Parameters
    ----------
    panel_id : str
        PanelApp ID of the panel
    panelapp_dict : dict
        PanelApp dump dict where keys are panel IDs

    Returns
    -------
    my_panel : dict
        dict of all the info for a particular panel
    """
    my_panel = panelapp_dict.get(panel_id)

    if not my_panel:
        print(
            f"WARNING: The panel ID {panel_id} was not found in the PanelApp "
            "JSON dump. This is expected if only HGNCs were given, otherwise"
            " please check that the panel ID is correct. No MOI-specific "
            "filtering will be performed"
        )

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
    if panel_data:
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
                # conf_level 3 indicates sufficient gene-disease association
                # evidence for use in variant interpretation (3 is green, 2
                # amber, 1 red)
                if conf_level >= 3:
                    panel_dict[region_name]['mode_of_inheritance'] = moi
                    panel_dict[region_name]['entity_type'] = 'region'
    else:
        print("WARNING - panel-specific dictionary from PanelApp is empty")

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
        if moi:
            xlr = "X-LINKED: hemizygous mutation in males, biallelic"
            xld = "X-LINKED: hemizygous mutation in males, monoallelic"
            if re.search(r"^BIALLELIC", moi):
                updated_moi = 'AR'
            elif re.search(r"^MONOALLELIC", moi):
                updated_moi = 'AD'
            elif re.search(r"^"+xlr, moi):
                updated_moi = 'XLR'
            elif re.search(r"^"+xld, moi):
                updated_moi = 'XLD'
            elif re.search(r"^BOTH", moi):
                updated_moi = 'AD/AR'
            elif re.search(r"^MITOCHONDRIAL", moi):
                updated_moi = 'MITOCHONDRIAL'
            elif re.search(r"^Other", moi):
                updated_moi = 'OTHER'
            elif re.search(r"^Unknown", moi):
                updated_moi = 'UNKNOWN'
            else:
                updated_moi = 'NONE'
        else:
            updated_moi = 'NONE'

        updated_gene_dict[gene]['mode_of_inheritance'] = updated_moi

    return updated_gene_dict


def get_formatted_dict(panel_string, genepanels_file, panelapp_file):
    """
    Main function to get a simple dictionary for each panel

    Parameters
    ----------
    panel_string : str
        The clinical indication string, including R number
    genepanels_file : str
        path to genepanels file which contains panel IDs
    panelapp_file : str
        path to the PanelApp JSON dump

    Returns
    -------
    final_panel_dict : dict
        default dict with each gene on the panel as key and the gene info as val
    """
    # Call functions to get PanelApp data from dump for given panel
    # and parse out the gene and region info
    genepanels_dict = parse_genepanels(genepanels_file)
    panel_id = get_panel_id_from_genepanels(panel_string, genepanels_dict)
    panel_dump = read_in_json(panelapp_file)
    panelapp_dict = transform_panelapp_dump_to_dict(panel_dump)
    single_panel_dict = parse_panelapp_dump(panel_id, panelapp_dict)
    formatted_single_panel = format_panel_info(single_panel_dict)
    final_panel_dict = simplify_MOI_terms(formatted_single_panel)

    return final_panel_dict
