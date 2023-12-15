import os
import pytest
import re
import sys

from collections import defaultdict

sys.path.append(os.path.abspath(
    os.path.join(os.path.realpath(__file__), '../../')
))

from tests import TEST_DATA_DIR
from utils import panels

TEST_GENEPANELS_CORRECT_FORMAT = "test_genepanels_formatted_correctly.tsv"
TEST_GENEPANELS_BAD_FORMAT = "test_genepanels_bad_format.tsv"


class TestParseGenePanels():
    """
    Test the parse_genepanels function which reads a tsv into a dictionary
    containing each clinical indication as key and PanelApp panel ID set as
    value
    """
    genepanels_tsv = os.path.join(
        TEST_DATA_DIR, TEST_GENEPANELS_CORRECT_FORMAT
    )
    genepanels_tsv2 = os.path.join(TEST_DATA_DIR, TEST_GENEPANELS_BAD_FORMAT)

    def test_parse_genepanels_when_correct_format(self):
        """
        Check that the TSV is parsed into a dict correctly
        """
        assert panels.parse_genepanels(self.genepanels_tsv) == {
            'C1.1_Inherited Stroke': {''},
            'C2.1_INSR': {''},
            ('R100.3_Rare syndromic craniosynostosis or isolated multisuture'
             ' synostosis_P'): {'168'},
            ('R101.1_Ehlers Danlos syndrome with a likely monogenic '
             'cause_P'): {'53'}
        }, "Genepanels not parsed correctly into dictionary"

    def test_parse_genepanels_when_multiple_panel_ids_for_one_panel(self):
        """
        Check raises error when same panel ID present for different panels
        """
        expected_error = re.escape(
            "Multiple panel IDs found for clinical indications: {'R341."
            "1_Hereditary angioedema types I and II_G': ['', '12'], 'R367."
            "1_Inherited pancreatic cancer_P': ['14', '524']}"
        )
        with pytest.raises(AssertionError, match=expected_error):
            panels.parse_genepanels(self.genepanels_tsv2)


class TestGetPanelIDFromGenePanels():
    """
    Test the get_panel_id_from_genepanels() function which takes
    a panel string (clinical indication) and gets the PanelApp ID for the panel
    """
    test_panels_dict = {
        'R107.1_Bardet Biedl syndrome_P': {'543'},
        'R109.3_Childhood onset leukodystrophy_P': {'496'},
        'R104.3_Skeletal dysplasia_P': {'504', '112'},
        'R97.1_Thrombophilia with a likely monogenic cause_P': {'504'},
        'R97.1_Thrombophilia_P': {''},
        'R23_Test_CI_P': {'504', ''}
    }

    def test_get_panel_id_from_genepanels_when_single_id_exists(self):
        """
        Checks the correct panel ID is returned when given a
        clinical indication which exists and only has one ID in the set
        """
        panel_id = panels.get_panel_id_from_genepanels(
            'R107.1_Bardet Biedl syndrome_P', self.test_panels_dict
        )
        assert panel_id == '543', (
            "Panel ID not returned correctly when one instance of the "
            "clinical indication exists"
        )

    def test_get_panel_id_for_ci_with_multiple_panel_ids(self):
        """
        Check assertion error is raised if > 1 panel ID found for the panel
        string
        """
        expected_error = (
            "Multiple panel IDs found for panel string: R104.3_Skeletal "
            "dysplasia_P"
        )
        with pytest.raises(AssertionError, match=expected_error):
            panels.get_panel_id_from_genepanels(
                'R104.3_Skeletal dysplasia_P', self.test_panels_dict
            )

    def test_get_panel_id_when_panel_id_is_empty_string(self):
        """
        Check assertion error is raised if panel ID is empty, i.e. just ''
        """
        expected_error = (
            "The clinical indication R97.1_Thrombophilia_P does not have a "
            "panel ID found in genepanels"
        )
        with pytest.raises(AssertionError, match=expected_error):
            panels.get_panel_id_from_genepanels(
                'R97.1_Thrombophilia_P', self.test_panels_dict
            )

    def test_when_more_than_one_panel_id_exists_plus_empty_string(self):
        """
        Assertion error should be raised if >1 panel ID exists for a CI
        """
        expected_error = (
            "Multiple panel IDs found for panel string: R23_Test_CI_P"
        )
        with pytest.raises(AssertionError, match=expected_error):
            panels.get_panel_id_from_genepanels(
                'R23_Test_CI_P', self.test_panels_dict
            )

    def test_when_panel_string_doesnt_exist_in_genepanels_dict(self, capsys):
        """
        Check that warning is printed if panel/HGNCs given which are not
        in the dict
        """
        expected_warning = (
            "WARNING: The panel string given None was not found"
            " in the genepanels file. This is expected if only HGNCs have been"
            " entered, but no MOI-specific filtering will be performed"
        )
        panels.get_panel_id_from_genepanels(
            None, self.test_panels_dict
        )
        stdout = capsys.readouterr().out

        assert expected_warning in stdout, (
            "Warning not raised correctly if panel/HGNCs given which were not"
            " in genepanels file "
        )


class TestTransformPanelAppDumpToDict():
    """
    Test the transform_panelapp_dump_to_dict() function, which takes a list of
    panel dictionaries and converts it to a dict where the panel IDs are keys,
    works correctly
    """
    # minimal version of PanelApp JSON dump that would be passed
    # into the function to be converted to a dict
    test_panel_dump = [
        {
            'panel_source': 'PanelApp',
            'panel_name': 'Disorders of sex development',
            'external_id': '9',
            'panel_version': '4.0',
            'genes': [
                {
                    'transcript': None,
                    'hgnc_id': 'HGNC:464',
                    'confidence_level': '3',
                    'mode_of_inheritance': (
                        'BIALLELIC, autosomal or pseudoautosomal'
                    ),
                    'mode_of_pathogenicity': None,
                    'penetrance': None,
                    'gene_justification': 'PanelApp',
                    'transcript_justification': 'PanelApp',
                    'alias_symbols': 'MIS',
                    'gene_symbol': 'AMH'
                },
            ]
        },
        {
            'panel_source': 'PanelApp',
            'panel_name': 'Disorders of sex development',
            'panel_version': '4.0',
            'genes': [
                {
                    'transcript': None,
                    'hgnc_id': 'HGNC:464',
                    'confidence_level': '3',
                    'mode_of_inheritance': (
                        'BIALLELIC, autosomal or pseudoautosomal'
                    ),
                    'mode_of_pathogenicity': None,
                    'penetrance': None,
                    'gene_justification': 'PanelApp',
                    'transcript_justification': 'PanelApp',
                    'alias_symbols': 'MIS',
                    'gene_symbol': 'AMH'
                },
            ]
        },
        {
            'panel_source': 'PanelApp',
            'panel_name': 'Disorders of sex development',
            'external_id': '',
            'panel_version': '4.0',
            'genes': [
                {
                    'transcript': None,
                    'hgnc_id': 'HGNC:464',
                    'confidence_level': '3',
                    'mode_of_inheritance': (
                        'BIALLELIC, autosomal or pseudoautosomal'
                    ),
                    'mode_of_pathogenicity': None,
                    'penetrance': None,
                    'gene_justification': 'PanelApp',
                    'transcript_justification': 'PanelApp',
                    'alias_symbols': 'MIS',
                    'gene_symbol': 'AMH'
                },
            ]
        }
    ]

    test_panel_dump2 = [
        {
            'panel_source': 'PanelApp',
            'panel_name': 'Disorders of sex development',
            'external_id': '',
            'panel_version': '4.0',
            'genes': [
                {
                    'transcript': None,
                    'hgnc_id': 'HGNC:464',
                    'confidence_level': '3',
                    'mode_of_inheritance': (
                        'BIALLELIC, autosomal or pseudoautosomal'
                    ),
                    'mode_of_pathogenicity': None,
                    'penetrance': None,
                    'gene_justification': 'PanelApp',
                    'transcript_justification': 'PanelApp',
                    'alias_symbols': 'MIS',
                    'gene_symbol': 'AMH'
                },
            ]
        }
    ]

    def test_transform_panelapp_dump_to_dict(self):
        """
        Make sure that the two panels with no external ID are not parsed
        into my final dict
        """
        assert panels.transform_panelapp_dump_to_dict(self.test_panel_dump) == {
            '9': {
                'panel_source': 'PanelApp',
                'panel_name': 'Disorders of sex development',
                'external_id': '9',
                'panel_version': '4.0',
                'genes': [{
                    'transcript': None,
                    'hgnc_id': 'HGNC:464',
                    'confidence_level': '3',
                    'mode_of_inheritance': (
                        'BIALLELIC, autosomal or pseudoautosomal'
                    ),
                    'mode_of_pathogenicity': None,
                    'penetrance': None,
                    'gene_justification': 'PanelApp',
                    'transcript_justification': 'PanelApp',
                    'alias_symbols': 'MIS',
                    'gene_symbol': 'AMH'
                }]
            }
        }, (
                "PanelApp dict not created correctly when PanelApp external_id"
                " key is missing"
        )

    def test_transform_panelapp_dump_to_dict_when_no_panels_left(self):
        """
        Make sure that error raised if no panels are left after transforming
        """
        expected_error = "No panels with IDs found in PanelApp dump"
        with pytest.raises(AssertionError, match=expected_error):
            panels.transform_panelapp_dump_to_dict(self.test_panel_dump2)


class TestParsePanelAppDump():
    """
    Make sure the parse_panelapp_dump() function works correctly, including
    when panel IDs don't exist
    """
    test_panel_dict = {
        '9': {
                'panel_source': 'PanelApp',
                'panel_name': 'Disorders of sex development',
                'external_id': '9',
                'panel_version': '4.0',
                'genes': [
                    {
                        'transcript': None,
                        'gene_symbol': 'TEST'
                    }
                ]
        },
        '12': {
                'panel_source': 'PanelApp',
                'panel_name': 'Test panel',
                'external_id': '12',
                'panel_version': '4.0',
                'genes': [
                    {
                        'transcript': None,
                        'gene_symbol': 'TEST2'
                    }
                ]
        }
    }

    test_empty_panel_dict = {}

    def test_parse_panelapp_dump_id_exists(self):
        """
        Check panel obtained correctly if exists
        """
        assert panels.parse_panelapp_dump('9', self.test_panel_dict) == (
            {
                'panel_source': 'PanelApp',
                'panel_name': 'Disorders of sex development',
                'external_id': '9',
                'panel_version': '4.0',
                'genes': [
                    {
                        'transcript': None,
                        'gene_symbol': 'TEST'
                    }
                ]
            }
        ), "Panel not obtained correctly from panel dict when key exists"

    def test_parse_panelapp_dump_id_not_exists(self, capsys):
        """
        Check error raised correctly if panel ID not in dict
        """
        panels.parse_panelapp_dump(None, self.test_panel_dict)
        stdout = capsys.readouterr().out
        expected_warning = (
            "WARNING: The panel ID None was not found in the PanelApp "
            "JSON dump. This is expected if only HGNCs were given, otherwise"
            " please check that the panel ID is correct"
        )
        assert expected_warning in stdout, (
            "Warning not printed correctly if empty panel ID given"
        )

    def test_parse_panelapp_dump_when_dump_empty(self, capsys):
        """
        Check warning raised if panel dict empty
        """
        panels.parse_panelapp_dump('9', self.test_empty_panel_dict)
        stdout = capsys.readouterr().out
        expected_warning = (
            "WARNING: The panel ID 9 was not found in the PanelApp "
            "JSON dump. This is expected if only HGNCs were given, otherwise"
            " please check that the panel ID is correct"
        )
        assert expected_warning in stdout, (
            "Warning not printed correctly if panel ID given but PanelApp"
            " is empty"
        )


class TestFormatPanelInfo():
    """
    Test that the format_panel_info function works correctly
    """
    test_one_gene_conf_level_below_3 = {
        'panel_source': 'PanelApp',
        'panel_name': 'Severe early-onset obesity',
        'external_id': '130',
        'panel_version': '4.0',
        'genes': [{
            'transcript': None,
            'hgnc_id': 'HGNC:428',
            'confidence_level': '1',
            'mode_of_inheritance': 'BIALLELIC, autosomal or pseudoautosomal',
            'mode_of_pathogenicity': None,
            'penetrance': 'Complete',
            'gene_justification': 'PanelApp',
            'transcript_justification': 'PanelApp',
            'alias_symbols': 'KIAA0328',
            'gene_symbol': 'ALMS1'
        }]
    }

    test_format_panel_dict_gene_and_region = {
        'panel_source': 'PanelApp',
        'panel_name': 'Severe early-onset obesity',
        'external_id': '130',
        'panel_version': '4.0',
        'genes': [
            {
                'transcript': None,
                'hgnc_id': 'HGNC:13210',
                'confidence_level': '3',
                'mode_of_inheritance': (
                    'BIALLELIC, autosomal or pseudoautosomal'
                ),
                'mode_of_pathogenicity': None,
                'penetrance': 'Complete',
                'gene_justification': 'PanelApp',
                'transcript_justification': 'PanelApp',
                'alias_symbols': 'RP55',
                'gene_symbol': 'ARL6'
            },
            {
                'transcript': None,
                'hgnc_id': 'HGNC:13210',
                'confidence_level': '3',
                'mode_of_inheritance': None,
                'mode_of_pathogenicity': None,
                'penetrance': 'Complete',
                'gene_justification': 'PanelApp',
                'transcript_justification': 'PanelApp',
                'alias_symbols': 'RP55',
                'gene_symbol': 'TEST'
            }
        ],
        'regions': [{
            'confidence_level': '3',
            'mode_of_inheritance': (
                'MONOALLELIC, autosomal or pseudoautosomal, imprinted status'
                ' unknown'
            ),
            'mode_of_pathogenicity': None,
            'penetrance': None,
            'name': (
                '15q11q13 recurrent (PWS/AS) region (BP1-BP3, Class 1) Loss'
            ),
            'chrom': '15',
            'start_37': None,
            'end_37': None,
            'start_38': 22782170,
            'end_38': 28134728,
            'type': 'CNV',
            'variant_type': 'cnv_loss',
            'required_overlap': 60,
            'haploinsufficiency': '3',
            'triplosensitivity': None,
            'justification': 'PanelApp'
        }]
    }

    test_panel_dict_not_exists = None

    def test_format_panel_info_if_no_genes_above_conf_3(self):
        """
        Check empty dict created if no genes above conf 3 found
        """
        assert not panels.format_panel_info(
            self.test_one_gene_conf_level_below_3
        ), (
            "Formatted panel dict not created correctly when no genes with "
            "confidence level above 3"
        )

    def test_format_panel_info_if_genes_and_regions_above_conf_3(self):
        """
        Check dict created correctly when genes and regions above conf level
        3 are present
        """
        assert panels.format_panel_info(
            self.test_format_panel_dict_gene_and_region
        ) == {
            'ARL6': {
                'mode_of_inheritance': (
                    'BIALLELIC, autosomal or pseudoautosomal'
                ),
                'entity_type': 'gene'
            },
            'TEST': {
                'mode_of_inheritance': None,
                'entity_type': 'gene'
            },
            '15q11q13 recurrent (PWS/AS) region (BP1-BP3, Class 1) Loss': {
                'mode_of_inheritance': (
                    'MONOALLELIC, autosomal or pseudoautosomal, imprinted '
                    'status unknown'
                ),
                'entity_type': 'region'
            }
        }, (
            "Formatted panel dict not created correctly when genes and regions"
            " above conf level 3 exist"
        )

    def test_if_panel_dict_not_exists(self, capsys):
        """
        Check warning printed correctly if trying to format empty dict
        """
        panels.format_panel_info(self.test_panel_dict_not_exists)
        stdout = capsys.readouterr().out
        expected_warning = (
            "WARNING - panel-specific dictionary from PanelApp is empty"
        )
        assert expected_warning in stdout, (
            "Warning not printed correctly if empty dict given to "
            "format_panel_info"
        )


class TestSimplifyMOITerms():
    """
    Test that simplify_MOI_terms function (which converts PanelApp MOI terms to
    simpler categories to be added to VCF as INFO fields) works as expected
    """
    test_gene_dict = {
        'gene1': {
            'mode_of_inheritance': 'BIALLELIC, autosomal or pseudoautosomal',
            'entity_type': 'gene'
        },
        'gene2': {
            'mode_of_inheritance': (
                'MONOALLELIC, autosomal or pseudoautosomal, NOT imprinted'
            ),
            'entity_type': 'gene'
        },
        'gene3': {
            'mode_of_inheritance': (
                'MONOALLELIC, autosomal or pseudoautosomal, imprinted status '
                'unknown'
            ),
            'entity_type': 'gene'
        },
        'gene4': {
            'mode_of_inheritance': (
                'MONOALLELIC, autosomal or pseudoautosomal, maternally '
                'imprinted (paternal allele expressed)'
            ),
            'entity_type': 'gene'
        },
        'gene5': {
            'mode_of_inheritance': (
                'MONOALLELIC, autosomal or pseudoautosomal, paternally '
                'imprinted (maternal allele expressed)'
            ),
            'entity_type': 'gene'
        },
        'gene6': {
            'mode_of_inheritance': None,
            'entity_type': 'gene'
        },
        'gene7': {
            'mode_of_inheritance': (
                'BOTH monoallelic and biallelic (but BIALLELIC mutations '
                'cause a more SEVERE disease form), autosomal or '
                'pseudoautosomal'
            ),
            'entity_type': 'gene'
        },
        'gene8': {
            'mode_of_inheritance': (
                'BOTH monoallelic and biallelic, autosomal or pseudoautosomal'
            ),
            'entity_type': 'gene'
        },
        'gene9': {
            'mode_of_inheritance': 'MITOCHONDRIAL',
            'entity_type': 'gene'
        },
        'gene10': {
            'mode_of_inheritance': 'Other',
            'entity_type': 'gene'
        },
        'gene11': {
            'mode_of_inheritance': (
                'Other - please specifiy in evaluation comments'
            ),
            'entity_type': 'gene'
        },
        'gene12': {
            'mode_of_inheritance': (
                'Other - please specifiy in evaluation comments'
            ),
            'entity_type': 'gene'
        },
        'gene13': {
            'mode_of_inheritance': (
                'Other - please specify in evaluation comments'
            ),
            'entity_type': 'gene'
        },
        'gene14': {
            'mode_of_inheritance': 'Unknown',
            'entity_type': 'gene'
        },
        'gene15': {
            'mode_of_inheritance': (
                'X-LINKED: hemizygous mutation in males, biallelic mutations '
                'in females'
            ),
            'entity_type': 'gene'
        },
        'gene16': {
            'mode_of_inheritance': (
                'X-LINKED: hemizygous mutation in males, monoallelic mutations'
                ' in females may cause disease (may be less severe, later '
                'onset than males)'
            ),
            'entity_type': 'gene'
        },
        'gene17': {
            'mode_of_inheritance': 'MOI we havent thought of',
            'entity_type': 'gene'
        }
    }

    test_empty_dict_if_only_hgncs = defaultdict(dict)

    def test_simplify_MOI_terms(self):
        """
        Check simplify_MOI_terms function takes dict and simplifies the
        PanelApp MOI terms simplified correctly, including when MOI is None
        """
        assert panels.simplify_MOI_terms(self.test_gene_dict) == {
            'gene1': {'mode_of_inheritance': 'AR'},
            'gene2': {'mode_of_inheritance': 'AD'},
            'gene3': {'mode_of_inheritance': 'AD'},
            'gene4': {'mode_of_inheritance': 'AD'},
            'gene5': {'mode_of_inheritance': 'AD'},
            'gene6': {'mode_of_inheritance': 'NONE'},
            'gene7': {'mode_of_inheritance': 'AD/AR'},
            'gene8': {'mode_of_inheritance': 'AD/AR'},
            'gene9': {'mode_of_inheritance': 'MITOCHONDRIAL'},
            'gene10': {'mode_of_inheritance': 'OTHER'},
            'gene11': {'mode_of_inheritance': 'OTHER'},
            'gene12': {'mode_of_inheritance': 'OTHER'},
            'gene13': {'mode_of_inheritance': 'OTHER'},
            'gene14': {'mode_of_inheritance': 'UNKNOWN'},
            'gene15': {'mode_of_inheritance': 'XLR'},
            'gene16': {'mode_of_inheritance': 'XLD'},
            'gene17': {'mode_of_inheritance': 'NONE'}
        }, "MOIs not simplified correctly"

    def test_simplify_MOI_terms_if_dict_is_none(self):
        """
        Test an empty dict is given back if an empty dict is given
        """
        assert panels.simplify_MOI_terms(
            self.test_empty_dict_if_only_hgncs
        ) == defaultdict(dict), (
            "Empty dict not returned correctly after MOI simplification"
            " when empty dict given"
        )
