import os
import pytest
import sys
import unittest

from pathlib import Path
from unittest.mock import Mock, patch

sys.path.append(os.path.abspath(
    os.path.join(os.path.realpath(__file__), '../../')
))

from tests import TEST_DATA_DIR
from utils import vcf


TEST_ANNOTATED_VCF = (
    "126560840-23326Q0015-23NGWES4-9526-F-103698_markdup_recalibrated_"
    "Haplotyper_annotated.vcf.gz"
)
TEST_SPLIT_VCF = (
    "126560840-23326Q0015-23NGWES4-9526-F-103698_markdup_recalibrated_"
    "Haplotyper_annotated.split_test.vcf"
)
TEST_FLAGGED_VCF = (
    "126560840-23326Q0015-23NGWES4-9526-F-103698_markdup_recalibrated_"
    "Haplotyper_annotated.flagged_test.vcf"
)

TEST_TRUNCATED_VCF = (
    "126560840-23326Q0015-23NGWES4-9526-F-103698_markdup_recalibrated_"
    "Haplotyper_annotated.flagged_truncated.vcf"
)


class TestBgzip(unittest.TestCase):
    """
    Test the function which uses subprocess to bgzip a file
    """
    annotated_split_vcf = os.path.join(TEST_DATA_DIR, TEST_SPLIT_VCF)

    def test_bgzip_output_file_exists(self):
        """
        Test that gzipped output file exists
        """
        vcf.bgzip(self.annotated_split_vcf)
        assert os.path.exists(
            os.path.join(
                TEST_DATA_DIR,
                "126560840-23326Q0015-23NGWES4-9526-F-103698_markdup_"
                "recalibrated_Haplotyper_annotated.split_test.vcf.gz"
            )
        ), "gzipped file does not exist"

        # Remove the gzipped VCF
        os.remove(os.path.join(
            TEST_DATA_DIR,
            "126560840-23326Q0015-23NGWES4-9526-F-103698_markdup_"
            "recalibrated_Haplotyper_annotated.split_test.vcf.gz"
        ))

    @patch('utils.vcf.subprocess.run')
    def test_bgzip_raises_error_if_return_code_not_zero(self, mock_vcf):
        """
        Test assertion error raised if return code of bgzip not zero
        """
        mock_vcf.return_value.returncode = 1

        with pytest.raises(AssertionError):
            vcf.bgzip(mock_vcf)


class TestBcftoolsPreProcess():
    """
    Test the function which uses subprocess to split VEP CSQ fields to
    separate INFO fields
    """
    annotated_vcf = os.path.join(TEST_DATA_DIR, TEST_ANNOTATED_VCF)

    def test_bcftools_pre_process_variant_count(self, capsys):
        """
        Test variant counts before and after bcftools +split-vep
        are printed as expected
        """
        output_vcf = vcf.bcftools_pre_process(self.annotated_vcf)
        stdout = capsys.readouterr().out

        errors = []
        if not 'Total lines before splitting: 255' in stdout:
            errors.append(
                "Variant counts pre-split not included in stdout as expected"
            )
        if not 'Total lines after splitting: 255' in stdout:
            errors.append(
                "Variant counts after split not included in stdout as expected"
            )
        if not output_vcf == (
            '126560840-23326Q0015-23NGWES4-9526-F-103698_markdup_recalibrated'
            '_Haplotyper_annotated.split.vcf'
        ):
            errors.append(
                "Output VCF from bcftools_pre_process not named as expected"
            )

        assert not errors, (
            "Errors occurred with the output of bcftools_preprocess"
            "():\n{}".format("\n".join(errors))
        )

        os.remove(
            '126560840-23326Q0015-23NGWES4-9526-F-103698_markdup_recalibrated'
            '_Haplotyper_annotated.split.vcf'
        )

    @patch('utils.vcf.subprocess.run')
    def test_bcftools_pre_process_raises_error_if_return_code_not_zero(
            self, mock_vcf
    ):
        """
        Test bcftools pre-process function which splits VEP vcf
        raises assertion error if return code not zero
        """
        mock_vcf.return_value.returncode = 124

        with pytest.raises(AssertionError):
            vcf.bcftools_pre_process(mock_vcf)

    @patch('utils.vcf.subprocess.run')
    def test_bcftools_pre_process_raises_error_if_variants_lost(
        self, mock_subprocess
    ):
        """
        Test assertion error raised if variant are lost following bcftools
        +split-vep
        """
        mock_subprocess.side_effect = [
            Mock(stdout=b'5'), Mock(returncode=0), Mock(stdout=b'14')
        ]

        with pytest.raises(AssertionError):
            vcf.bcftools_pre_process('mock_vcf')


class TestReadInVCF():
    """
    Test the VCF (with split VCF fields) is read in with pysam correctly
    """
    annotated_split_vcf = os.path.join(TEST_DATA_DIR, TEST_SPLIT_VCF)

    # Read in the control VCF with pysam and get sample name and CSQ fields
    vcf_contents, sample_name = vcf.read_in_vcf(
        annotated_split_vcf
    )

    def test_read_in_vcf_sample_parsed_correctly(self):
        """
        Check VCF read in to pysam correctly
        """

        # Check sample name parsed correctly
        assert self.sample_name == (
            '126560840-23326Q0015-23NGWES4-9526-F-103698'
        )

    def test_read_in_vcf_adds_info_header_correctly(self):
        """
        Test that the MOI INFO tag is added as a new header line as expected
        """
        # Get all of the pysam header records
        vcf_header_items = [
            record.values() for record in self.vcf_contents.header.records
        ]

        assert [
            'MOI', '.', 'String',
            '"Mode of inheritance from PanelApp (simplified)"', '131'
        ] in vcf_header_items, "MOI not added to header correctly"


class TestAddMOIFlag():
    """
    Test that filtering flag added to variants correctly
    """
    vcf_contents, _ = vcf.read_in_vcf(
        os.path.join(TEST_DATA_DIR, TEST_SPLIT_VCF)
    )

    test_panel_dict = {'POMC': {'mode_of_inheritance': 'AR'}}

    gene_variant_dict = vcf.add_MOI_field(vcf_contents, test_panel_dict)

    def test_add_MOI_check_MOI_added_correctly_for_present_gene(self):
        """
        Assert that the 2 variants present in POMC in the test VCF both
        have 'AR' as their MOI
        """
        assert [
            record.info['MOI'] for record in self.gene_variant_dict.get('POMC')
        ] == [('AR', ), ('AR',)], (
            "MOI not added correctly as AR for the two variants in POMC"
        )

    def test_MOIs_added_as_unknown_when_not_in_dict(self):
        """
        Assert that variants in all other genes not in the panel dict
        have MOI INFO field added as 'NONE'
        """
        all_mois_not_in_panel_dict = []
        for gene, variant_list in self.gene_variant_dict.items():
            if gene != 'POMC':
                all_mois_for_gene = [
                    variant.info['MOI'] for variant in variant_list
                ]
                all_mois_not_in_panel_dict.append(all_mois_for_gene)

        assert all(('NONE',) for moi in all_mois_not_in_panel_dict)


class TestWriteOutFlaggedVCF():
    """
    Test writing out the pysam object as a VCF file works as expected
    """
    vcf_contents, _ = vcf.read_in_vcf(
        os.path.join(TEST_DATA_DIR, TEST_SPLIT_VCF)
    )

    test_panel_dict = {'POMC': {'mode_of_inheritance': 'AR'}}

    flagged_vcf = (
        '126560840-23326Q0015-23NGWES4-9526-F-103698_markdup_recalibrated'
        '_Haplotyper_annotated.flagged.vcf'
    )

    gene_variant_dict = vcf.add_MOI_field(vcf_contents, test_panel_dict)

    def test_write_out_flagged_vcf(self):
        """
        Test that the write_out_flagged_vcf function creates a flagged VCF file
        as expected
        """
        vcf.write_out_flagged_vcf(
            self.flagged_vcf, self.gene_variant_dict, self.vcf_contents
        )

        assert os.path.exists(self.flagged_vcf)

        os.remove(self.flagged_vcf)


class TestCheckWrittenOutVcf():
    """
    Test that the test_check_written_out_vcf() function which checks that the
    pysam header and variants written to file match what was given to the
    function to write out (i.e. the file is not truncated)
    """
    # Read in control VCF which has had CSQ fields expanded by bcftools
    # +split-vep
    original_vcf_contents, _ = vcf.read_in_vcf(
        os.path.join(TEST_DATA_DIR, TEST_SPLIT_VCF)
    )
    # Create a test gene panel dict for obesity for adding MOI info to
    # variants
    test_panel_dict = {
        'ALMS1': {'mode_of_inheritance': 'AR'},
        'ARL6': {'mode_of_inheritance': 'AR'},
        'BBS1': {'mode_of_inheritance': 'AR'},
        'BBS10': {'mode_of_inheritance': 'AR'},
        'BBS12': {'mode_of_inheritance': 'AR'},
        'BBS2': {'mode_of_inheritance': 'AR'},
        'BBS4': {'mode_of_inheritance': 'AR'},
        'BBS5': {'mode_of_inheritance': 'AR'},
        'BBS7': {'mode_of_inheritance': 'AR'},
        'BBS9': {'mode_of_inheritance': 'AR'},
        'CEP19': {'mode_of_inheritance': 'AR'},
        'CPE': {'mode_of_inheritance': 'AR'},
        'GNAS': {'mode_of_inheritance': 'AD'},
        'KIDINS220': {'mode_of_inheritance': 'AD'},
        'LEP': {'mode_of_inheritance': 'AR'},
        'LEPR': {'mode_of_inheritance': 'AR'},
        'MC4R': {'mode_of_inheritance': 'AD/AR'},
        'MKKS': {'mode_of_inheritance': 'AR'},
        'MKS1': {'mode_of_inheritance': 'AR'},
        'MYT1L': {'mode_of_inheritance': 'AD'},
        'NTRK2': {'mode_of_inheritance': 'AD'},
        'PCSK1': {'mode_of_inheritance': 'AR'},
        'PGM2L1': {'mode_of_inheritance': 'AR'},
        'PHF6': {'mode_of_inheritance': 'XLD'},
        'PHIP': {'mode_of_inheritance': 'AD'},
        'POMC': {'mode_of_inheritance': 'AR'},
        'SDCCAG8': {'mode_of_inheritance': 'AR'},
        'SIM1': {'mode_of_inheritance': 'AD'},
        'TTC8': {'mode_of_inheritance': 'AR'},
        'VPS13B': {'mode_of_inheritance': 'AR'},
        '15q11q13 recurrent (PWS/AS) region (BP1-BP3, Class 1) Loss': {
            'mode_of_inheritance': 'AD'
        },
        '15q11q13 recurrent (PWS/AS) region (BP2-BP3, Class 2) Loss': {
            'mode_of_inheritance': 'AD'
        },
        (
            '16p11.2 recurrent region (includes SH2B1) (distal region) '
            '(BP2-BP3) Loss'
        ): {'mode_of_inheritance': 'AD'}
    }

    gene_variant_dict = vcf.add_MOI_field(
        original_vcf_contents, test_panel_dict
    )
    # This VCF has one variant removed from the end
    truncated_vcf = os.path.join(TEST_DATA_DIR, TEST_TRUNCATED_VCF)

    def test_check_written_out_vcf_raises_error(self):
        """
        Test error is raised if VCF which was written out which is different/
        truncated compared to what was supposed to be written out
        """
        with pytest.raises(AssertionError):
            vcf.check_written_out_vcf(
                self.original_vcf_contents,
                self.gene_variant_dict,
                self.truncated_vcf
            )


class TestBCftoolsSort(unittest.TestCase):
    """
    Test that bcftools sort works as expected
    """
    @patch('utils.vcf.subprocess.run')
    def test_bcftools_sort_raises_error_if_return_code_not_zero(
        self, mock_subprocess
    ):
        """
        Test assertion error raised if return code of bcftools sort not zero
        """
        mock_subprocess.return_value.returncode = 2

        with pytest.raises(AssertionError):
            vcf.bcftools_sort('input_vcf')

    @patch('utils.vcf.subprocess.run')
    def test_bcftools_sort_raises_error_if_variant_counts_not_match(
        self, mock_subprocess
    ):
        """
        Test assertion error raised if variant counts pre- and post-bcftools
        sort do not match
        """
        mock_subprocess.side_effect = [
            Mock(stdout=b'14'), Mock(returncode=0), Mock(stdout=b'5')
        ]

        with pytest.raises(AssertionError):
            vcf.bcftools_sort('flagged_vcf')


class TestBcftoolsFilter(unittest.TestCase):
    """
    Test the function which uses subprocess to run bcftools filtering
    """
    flagged_vcf = os.path.join(TEST_DATA_DIR, TEST_FLAGGED_VCF)
    filter_command = (
        "bcftools filter --soft-filter \"EXCLUDE\" -m + "
        "-e '(CSQ_Consequence~\"synonymous_variant\")'"
    )
    filter_vcf = (
        f"{Path(flagged_vcf).stem.split('.')[0]}.optimised_filtered.vcf"
    )

    def test_bcftools_filter_creates_file(self):
        """
        Test that bcftools filter output file exists
        """
        vcf.bcftools_filter(
            self.flagged_vcf, self.filter_command, self.filter_vcf
        )

        # Check exists
        assert os.path.exists(
            self.filter_vcf
        ), "bcftools filter output file does not exist"

        # Remove file so doesn't affect test if run again in future
        os.remove(self.filter_vcf)


    @patch('utils.vcf.subprocess.run')
    def test_bcftools_filter_raises_error_if_return_code_not_zero(
        self, mock_subprocess
    ):
        """
        Test assertion error raised if return code of bcftools filter not zero
        """
        mock_subprocess.return_value.returncode = 2

        with pytest.raises(AssertionError):
            vcf.bcftools_filter('flag_vcf', 'filter_command', 'filter_vcf')

    @patch('utils.vcf.subprocess.run')
    def test_bcftools_filter_raises_error_if_variant_counts_not_match(
        self, mock_subprocess
    ):
        """
        Test assertion error raised if variant counts pre- and post-bcftools
        filter do not match
        """
        mock_subprocess.side_effect = [
            Mock(stdout=b'14'), Mock(returncode=0), Mock(stdout=b'5')
        ]

        with pytest.raises(AssertionError):
            vcf.bcftools_filter('flag_vcf', 'filter_command', 'filter_vcf')
