import os
import pytest
import sys

sys.path.append(os.path.abspath(
    os.path.join(os.path.realpath(__file__), '../../')
))

from utils import panels
from tests import TEST_DATA_DIR
from add_optimised_filtering import check_panel_string


@pytest.mark.parametrize("test_input,expected", [
    ('_HGNC:16627;_HGNC:795', ''),
    ('_HGNC:16627', ''),
    ('_HGNC:16627;R49.3_Beckwith-Wiedemann syndrome_G', 'R49.3_Beckwith-Wiedemann syndrome_G'),
    ('R49.3_Beckwith-Wiedemann syndrome_G;_HGNC:16627', 'R49.3_Beckwith-Wiedemann syndrome_G'),
    ('R49.3_Beckwith-Wiedemann syndrome_G;_HGNC:123;_HGNC:456', 'R49.3_Beckwith-Wiedemann syndrome_G'),
    ('_HGNC:123;_HGNC:456;R49.3_Beckwith-Wiedemann syndrome_G', 'R49.3_Beckwith-Wiedemann syndrome_G'),
    ('_HGNC:456;R49.3_Beckwith-Wiedemann syndrome_G;_HGNC:123;', 'R49.3_Beckwith-Wiedemann syndrome_G')
])
class TestCheckPanelString():
    """
    Tests that check that the panel string given is valid and contains either
    only HGNC IDs, multiple HGNC IDs plus one panel string or one panel string
    """

    def test_check_string_one_panel(self, test_input, expected):
        """
        Test that HGNCs are stripped from the input string as expected
        """
        assert check_panel_string(test_input) == expected

    # multiple_panels_one_hgnc_before = '_HGNC:16627;R49.3_Beckwith-Wiedemann syndrome_G;R228.1_Tuberous sclerosis_G'
    # multiple_panels_one_hgnc_after = 'R49.3_Beckwith-Wiedemann syndrome_G;R228.1_Tuberous sclerosis_G;_HGNC:16627'
    # multiple_panels_one_hgnc_inbetween = 'R49.3_Beckwith-Wiedemann syndrome_G;_HGNC:16627;R228.1_Tuberous sclerosis_G'
    # multiple_panels_multiple_hgncs_before = '_HGNC:16627;_HGNC:795;R49.3_Beckwith-Wiedemann syndrome_G;R228.1_Tuberous sclerosis_G'
    # multiple_panels_multiple_hgncs_after = 'R49.3_Beckwith-Wiedemann syndrome_G;R228.1_Tuberous sclerosis_G;_HGNC:16627;_HGNC:795'
    # multiple_panels_multiple_hgncs_inbetween = 'R49.3_Beckwith-Wiedemann syndrome_G;_HGNC:16627;R228.1_Tuberous sclerosis_G;_HGNC:795'

    # def test_check_panel_string_just_multiple_hgncs(self):
    #     """
    #     Check that error raised if
    #     """
    #     with pytest.raises(AssertionError):
    #         check_panel_string(self.multiple_panels_multiple_hgncs_before)


    # def test_check_panel_string_only_hgncs(self):
    #     assert check_panel_string(self.just_multiple_hgncs) == ''
