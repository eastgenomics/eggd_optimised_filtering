import os
import pytest
import sys

sys.path.append(os.path.abspath(
    os.path.join(os.path.realpath(__file__), '../../')
))

from add_optimised_filtering import check_panel_string


class TestCheckPanelString():
    """
    Tests that check that the panel string given is valid and contains either
    only HGNC IDs, multiple HGNC IDs plus one panel string or one panel string
    only
    """
    @pytest.mark.parametrize("test_input,expected", [
        ('_HGNC:16627;_HGNC:795', ''),
        ('_HGNC:16627', ''),
        ('_HGNC:16627;R49.3_Beckwith-Wiedemann syndrome_G',
         'R49.3_Beckwith-Wiedemann syndrome_G'),
        ('R49.3_Beckwith-Wiedemann syndrome_G;_HGNC:16627',
         'R49.3_Beckwith-Wiedemann syndrome_G'),
        ('R49.3_Beckwith-Wiedemann syndrome_G;_HGNC:123;_HGNC:456',
         'R49.3_Beckwith-Wiedemann syndrome_G'),
        ('_HGNC:123;_HGNC:456;R49.3_Beckwith-Wiedemann syndrome_G',
         'R49.3_Beckwith-Wiedemann syndrome_G'),
        ('_HGNC:456;R49.3_Beckwith-Wiedemann syndrome_G;_HGNC:123;',
         'R49.3_Beckwith-Wiedemann syndrome_G')
    ])
    def test_check_panel_string_one_panel_given(self, test_input, expected):
        """
        Test that HGNCs are stripped from the input string to give only
        the panel name or an empty string if only HGNCs, as expected
        """
        assert check_panel_string(test_input) == expected

    @pytest.mark.parametrize(
        "test_input_bad,expected_warning", [
            (('_HGNC:16627;R49.3_Beckwith-Wiedemann syndrome_G;R228.1'
            '_Tuberous sclerosis_G'),
            ('More than one panel given: R49.3_Beckwith-Wiedemann syndrome_G'
             ';R228.1_Tuberous sclerosis_G')),
            (('R49.3_Beckwith-Wiedemann syndrome_G;R228.1_Tuberous sclerosis'
              '_G;_HGNC:16627'),
             ('More than one panel given: R49.3_Beckwith-Wiedemann syndrome'
              '_G;R228.1_Tuberous sclerosis_G')),
            (('R49.3_Beckwith-Wiedemann syndrome_G;_HGNC:16627;R228.1'
              '_Tuberous sclerosis_G'),
             ('More than one panel given: R49.3_Beckwith-Wiedemann syndrome'
              '_G;R228.1_Tuberous sclerosis_G')),
            (('_HGNC:16627;_HGNC:795;R49.3_Beckwith-Wiedemann syndrome_G;'
              'R228.1_Tuberous sclerosis_G'),
             ('More than one panel given: R49.3_Beckwith-Wiedemann syndrome'
              '_G;R228.1_Tuberous sclerosis_G')),
            (('R49.3_Beckwith-Wiedemann syndrome_G;R228.1_Tuberous sclerosis'
              '_G;_HGNC:16627;_HGNC:795'),
             ('More than one panel given: R49.3_Beckwith-Wiedemann syndrome'
              '_G;R228.1_Tuberous sclerosis_G')),
            (('R49.3_Beckwith-Wiedemann syndrome_G;_HGNC:16627;R228.1_'
              'Tuberous sclerosis_G;_HGNC:795'),
             ('More than one panel given: R49.3_Beckwith-Wiedemann syndrome_G'
              ';R228.1_Tuberous sclerosis_G'))
        ]
    )
    def test_check_panel_string_multiple_panels_given(
        self, test_input_bad, expected_warning
    ):
        """
        Check warning is raised when more than one gene panel is given
        as the input string (multiple HGNCs are fine)
        """
        with pytest.raises(AssertionError, match=expected_warning):
            check_panel_string(test_input_bad)
