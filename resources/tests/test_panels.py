import os
import pytest
import sys

sys.path.append(os.path.abspath(
    os.path.join(os.path.realpath(__file__), '../../')
))

from tests import TEST_DATA_DIR
from utils import panels

TEST_GENEPANELS = "test_genepanels.tsv"
TEST_PANELDUMP = ""

class TestParseGenePanels():
    """
    Test the parse_genepanels function which reads a tsv into a dictionary
    containing each clinical indication as key and PanelApp panel ID set as
    value
    """
    genepanels_tsv = os.path.join(TEST_DATA_DIR, TEST_GENEPANELS)

    def test_parsed_correctly(self):
        """
        Check that the TSV is parsed into a dict correctly
        """
        assert panels.parse_genepanels(self.genepanels_tsv) == {
            'C1.1_Inherited Stroke': {''},
            'C2.1_INSR': {''},
            'R100.3_Rare syndromic craniosynostosis or isolated multisuture synostosis_P': {'168'},
            'R101.1_Ehlers Danlos syndrome with a likely monogenic cause_P': {'53'},
            'R341.1_Hereditary angioedema types I and II_G': {'', '12'},
            'R366.1_Inherited susceptibility to acute lymphoblastoid leukaemia (ALL)_P': {''},
            'R367.1_Inherited pancreatic cancer_P': {'524'}
        }, "Genepanels not parsed correctly into dictionary"


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
            "The clinical indication does not have a panel ID found in "
            "genepanels: R97.1_Thrombophilia_P"
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
