import os
import pytest
import sys

sys.path.append(os.path.abspath(
    os.path.join(os.path.realpath(__file__), '../../')
))

from tests import TEST_DATA_DIR
from utils import panels

TEST_GENEPANELS_CORRECT_FORMAT = "test_genepanels_formatted_correctly.tsv"
TEST_GENEPANELS_BAD_FORMAT = "test_genepanels_bad_format.tsv"
TEST_PANELDUMP = ""

class TestParseGenePanels():
    """
    Test the parse_genepanels function which reads a tsv into a dictionary
    containing each clinical indication as key and PanelApp panel ID set as
    value
    """
    genepanels_tsv = os.path.join(TEST_DATA_DIR, TEST_GENEPANELS_CORRECT_FORMAT)
    genepanels_tsv2 = os.path.join(TEST_DATA_DIR, TEST_GENEPANELS_BAD_FORMAT)

    def test_parse_genepanels_when_correct_format(self):
        """
        Check that the TSV is parsed into a dict correctly
        """
        assert panels.parse_genepanels(self.genepanels_tsv) == {
            'C1.1_Inherited Stroke': {''},
            'C2.1_INSR': {''},
            'R100.3_Rare syndromic craniosynostosis or isolated multisuture synostosis_P': {'168'},
            'R101.1_Ehlers Danlos syndrome with a likely monogenic cause_P': {'53'}
        }, "Genepanels not parsed correctly into dictionary"


    def test_parse_genepanels_when_multiple_panel_ids_for_one_panel(self):
        """
        Check raises error when same panel ID present for different panels
        """
        expected_error = (
            "Multiple panel IDs found for clinical indications: {'R341."
            "1_Hereditary angioedema types I and II_G': {'', '12'}, 'R367."
            "1_Inherited pancreatic cancer_P': {'524', '14'}}"
        )
        with pytest.raises(AssertionError, match=expected_error):
            panels.parse_genepanels(self.genepanels_tsv2)


class TestGetPanelIDFromGenePanels():
    """
    Test the get_panel_id_from_genepanels function which takes
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
        Assertion error should be raised if > 1 panel ID found
        """
        expected_error = "Multiple panel IDs found for panel string: R104.3_Skeletal dysplasia_P"
        with pytest.raises(AssertionError, match=expected_error):
            panels.get_panel_id_from_genepanels(
                'R104.3_Skeletal dysplasia_P', self.test_panels_dict
            )

    def test_get_panel_id_when_panel_id_is_empty_string(self):
        """
        Assertion error should be raised if panel ID is just '""'
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
        Assertion error should be raised if > panel ID exists for a CI
        """
        expected_error = (
            "Multiple panel IDs found for panel string: R23_Test_CI_P"
        )
        with pytest.raises(AssertionError, match=expected_error):
            panels.get_panel_id_from_genepanels(
                'R23_Test_CI_P', self.test_panels_dict
            )

    def test_when_panel_string_doesnt_exist_in_dict(self):
        """
        KeyError should be raised if panel name not found in dict
        """
        expected_error = (
            "The panel string was not found in the genepanels file: R22"
        )
        with pytest.raises(KeyError, match=expected_error):
            panels.get_panel_id_from_genepanels('R22', self.test_panels_dict)


class TestTransformPanelAppDumpToDict():
    """
    Test the transform_panelapp_dump_to_dict() function, which takes a list of
    panel dictionaries and converts it to a dict where the panel IDs are keys,
    works correctly
    """
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
                    'mode_of_inheritance': 'BIALLELIC, autosomal or pseudoautosomal',
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
                    'mode_of_inheritance': 'BIALLELIC, autosomal or pseudoautosomal',
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
                    'mode_of_inheritance': 'BIALLELIC, autosomal or pseudoautosomal',
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
                    'mode_of_inheritance': 'BIALLELIC, autosomal or pseudoautosomal',
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
                        'gene_symbol': 'TEST'
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

    def test_parse_panelapp_dump_id_not_exists(self):
        """
        Check error raised correctly if panel ID not in dict
        """
        expected_error = (
            'The panel ID 11 was not found in the PanelApp JSON dump'
        )
        with pytest.raises(KeyError, match=expected_error):
            panels.parse_panelapp_dump('11', self.test_panel_dict)

class TestSimplifyMOITerms():
    """
    _summary_
    """
    test_gene_dict = {

    }
