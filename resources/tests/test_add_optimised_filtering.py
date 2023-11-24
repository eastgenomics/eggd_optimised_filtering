import os
import pytest
import sys

sys.path.append(os.path.abspath(
    os.path.join(os.path.realpath(__file__), '../../')
))

from utils import panels
from tests import TEST_DATA_DIR
from add_optimised_filtering import check_panel_string


class TestCheckPanelString():
    """
    Tests for checking that the panel string given is checked correctly
    """
    # Test panel string examples
    just_multiple_hgncs = '_HGNC:16627;_HGNC:795'
    just_one_hgnc = '_HGNC:16627'

    one_panel_one_hgnc_before = '_HGNC:16627;R49.3_Beckwith-Wiedemann syndrome_G'
    one_panel_one_hgnc_after = 'R49.3_Beckwith-Wiedemann syndrome_G;_HGNC:16627'
    one_panel_multiple_hgncs_after = 'R49.3_Beckwith-Wiedemann syndrome_G;_HGNC:123;_HGNC:456'
    one_panel_multiple_hgncs_before = '_HGNC:123;_HGNC:456;R49.3_Beckwith-Wiedemann syndrome_G'
    one_panel_multiple_hgncs_around = '_HGNC:456;R49.3_Beckwith-Wiedemann syndrome_G;_HGNC:123;'

    multiple_panels_one_hgnc_before = '_HGNC:16627;R49.3_Beckwith-Wiedemann syndrome_G;R228.1_Tuberous sclerosis_G'
    multiple_panels_one_hgnc_after = 'R49.3_Beckwith-Wiedemann syndrome_G;R228.1_Tuberous sclerosis_G;_HGNC:16627'
    multiple_panels_one_hgnc_inbetween = 'R49.3_Beckwith-Wiedemann syndrome_G;_HGNC:16627;R228.1_Tuberous sclerosis_G'
    multiple_panels_multiple_hgncs_before = '_HGNC:16627;_HGNC:795;R49.3_Beckwith-Wiedemann syndrome_G;R228.1_Tuberous sclerosis_G'
    multiple_panels_multiple_hgncs_after = 'R49.3_Beckwith-Wiedemann syndrome_G;R228.1_Tuberous sclerosis_G;_HGNC:16627;_HGNC:795'
    multiple_panels_multiple_hgncs_inbetween = 'R49.3_Beckwith-Wiedemann syndrome_G;_HGNC:16627;R228.1_Tuberous sclerosis_G;_HGNC:795'

    def test_just_multiple_hgncs(self):
        with pytest.raises(AssertionError):
            check_panel_string(self.multiple_panels_multiple_hgncs_before)


# test_csq_fields = "SYMBOL,Consequence,gnomADe_AF,gnomADg_AF,TWE_AF,ClinVar_CLNSIG,ClinVar_CLNSIGCONF,SpliceAI_pred_DS_AG,SpliceAI_pred_DS_AL,SpliceAI_pred_DS_DG,SpliceAI_pred_DS_DL,HGMD_CLASS"
