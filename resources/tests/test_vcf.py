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
    "Haplotyper_annotated.vcf.split.vcf"
)
TEST_FLAGGED_VCF = (
    "126560840-23326Q0015-23NGWES4-9526-F-103698_markdup_recalibrated_"
    "Haplotyper_annotated.vcf.flagged.vcf"
)
TEST_FILTER_VCF = (
    "126560840-23326Q0015-23NGWES4-9526-F-103698_markdup_recalibrated_"
    "Haplotyper_annotated.vcf.filter.vcf"
)
TEST_TRUNCATED_VCF = (
    "126560840-23326Q0015-23NGWES4-9526-F-103698_markdup_recalibrated_"
    "Haplotyper_annotated.vcf.flagged_truncated.vcf"
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
                "recalibrated_Haplotyper_annotated.vcf.split.vcf.gz"
            )
        ), "gzipped file does not exist"

        # Remove the gzipped VCF
        os.remove(os.path.join(
            TEST_DATA_DIR,
            "126560840-23326Q0015-23NGWES4-9526-F-103698_markdup_"
            "recalibrated_Haplotyper_annotated.vcf.split.vcf.gz"
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
            '_Haplotyper_annotated.vcf.split.vcf'
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
            '_Haplotyper_annotated.vcf.split.vcf'
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
    vcf_contents, sample_name, csq_fields_to_collapse = vcf.read_in_vcf(
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

    def test_read_in_vcf_fields_to_collapse_parsed_correctly(self):
        """
        Check the VEP fields are parsed and converted to fields to be
        collapsed later correctly
        """
        assert self.csq_fields_to_collapse == (
            'INFO/CSQ_Allele,INFO/CSQ_SYMBOL,INFO/CSQ_HGNC_ID,INFO/'
            'CSQ_VARIANT_CLASS,INFO/CSQ_Consequence,INFO/CSQ_IMPACT,INFO/'
            'CSQ_EXON,INFO/CSQ_INTRON,INFO/CSQ_Feature,INFO/CSQ_HGVSc,INFO/'
            'CSQ_HGVSp,INFO/CSQ_HGVS_OFFSET,INFO/CSQ_Existing_variation,INFO/'
            'CSQ_STRAND,INFO/CSQ_ClinVar,INFO/CSQ_ClinVar_CLNSIG,INFO/'
            'CSQ_ClinVar_CLNSIGCONF,INFO/CSQ_ClinVar_CLNDN,INFO/'
            'CSQ_gnomADg_AC,INFO/CSQ_gnomADg_AN,INFO/CSQ_gnomADg_AF,INFO/'
            'CSQ_gnomADg_nhomalt,INFO/CSQ_gnomADg_popmax,INFO/'
            'CSQ_gnomADg_AC_popmax,INFO/CSQ_gnomADg_AN_popmax,INFO/'
            'CSQ_gnomADg_AF_popmax,INFO/CSQ_gnomADg_nhomalt_popmax,INFO/'
            'CSQ_gnomADe_AC,INFO/CSQ_gnomADe_AN,INFO/CSQ_gnomADe_AF,INFO/'
            'CSQ_gnomADe_nhomalt,INFO/CSQ_gnomADe_popmax,INFO/'
            'CSQ_gnomADe_AC_popmax,INFO/CSQ_gnomADe_AN_popmax,INFO/'
            'CSQ_gnomADe_AF_popmax,INFO/CSQ_gnomADe_nhomalt_popmax,INFO/'
            'CSQ_gnomADe_non_cancer_AC,INFO/CSQ_gnomADe_non_cancer_AN,INFO/'
            'CSQ_gnomADe_non_cancer_AF,INFO/CSQ_gnomADe_non_cancer_nhomalt,'
            'INFO/CSQ_gnomADe_non_cancer_AC_popmax,INFO/'
            'CSQ_gnomADe_non_cancer_AN_popmax,INFO/'
            'CSQ_gnomADe_non_cancer_AF_popmax,INFO/'
            'CSQ_gnomADe_non_cancer_nhomalt_popmax,INFO/'
            'CSQ_gnomADe_non_cancer_popmax,INFO/CSQ_TWE_AF,INFO/'
            'CSQ_TWE_AC_Hom,INFO/CSQ_TWE_AC_Het,INFO/CSQ_TWE_AN,'
            'INFO/CSQ_HGMD,INFO/CSQ_HGMD_PHEN,INFO/CSQ_HGMD_CLASS,'
            'INFO/CSQ_HGMD_RANKSCORE,INFO/CSQ_SpliceAI_pred_DS_AG,INFO/'
            'CSQ_SpliceAI_pred_DS_AL,INFO/CSQ_SpliceAI_pred_DS_DG,INFO/'
            'CSQ_SpliceAI_pred_DS_DL,INFO/CSQ_SpliceAI_pred_DP_AG,INFO/'
            'CSQ_SpliceAI_pred_DP_AL,INFO/CSQ_SpliceAI_pred_DP_DG,INFO/'
            'CSQ_SpliceAI_pred_DP_DL,INFO/CSQ_REVEL,INFO/CSQ_CADD_PHRED'
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
    vcf_contents, _, _ = vcf.read_in_vcf(
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
    vcf_contents, _, _ = vcf.read_in_vcf(
        os.path.join(TEST_DATA_DIR, TEST_SPLIT_VCF)
    )

    test_panel_dict = {'POMC': {'mode_of_inheritance': 'AR'}}

    flagged_vcf = (
        '126560840-23326Q0015-23NGWES4-9526-F-103698_markdup_recalibrated'
        '_Haplotyper_annotated.vcf.flagged.vcf'
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
    original_vcf_contents, _, _ = vcf.read_in_vcf(
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


class TestBcftoolsFilter(unittest.TestCase):
    """
    Test the function which uses subprocess to run bcftools filtering
    """
    flagged_vcf = os.path.join(TEST_DATA_DIR, TEST_FLAGGED_VCF)
    filter_command = (
        "bcftools filter --soft-filter \"EXCLUDE\" -m + "
        "-e '(CSQ_Consequence~\"synonymous_variant\")'"
    )
    filter_vcf = f"{Path(flagged_vcf).stem}.filter.vcf"

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


class TestBcftoolsRemoveCsqAnnotation(unittest.TestCase):
    """
    Test the function which uses subprocess to run bcftools to remove the
    split CSQ INFO fields from VEP (otherwise already split fields would
    break eggd_generate_variant_workbook)
    """
    filter_vcf = os.path.join(TEST_DATA_DIR, TEST_FILTER_VCF)
    fields_to_collapse = 'INFO/CSQ_Allele,INFO/CSQ_SYMBOL'
    resulting_vcf = f"{Path(filter_vcf).stem}.G2P.vcf"

    def test_bcftools_remove_csq_annotation_creates_file(self):
        """
        Test a file is created as expected when bcftools annotate is run
        """
        vcf.bcftools_remove_csq_annotation(
            self.filter_vcf, self.fields_to_collapse
        )
        assert os.path.exists(
                self.resulting_vcf
        ), "bcftools annotate output file does not exist"

        os.remove(self.resulting_vcf)

    @patch('utils.vcf.subprocess.run')
    def test_bcftools_remove_csq_annotation_error_if_return_code_not_zero(
        self, mock_subprocess
    ):
        """
        Test assertion error raised if return code of bcftools annotate not
        zero
        """
        mock_subprocess.return_value.returncode = 2

        with pytest.raises(AssertionError):
            vcf.bcftools_remove_csq_annotation(
                'filter_vcf', 'csq_fields_to_drop'
            )

    @patch('utils.vcf.subprocess.run')
    def test_bcftools_remove_csq_annotation_error_if_variant_counts_not_match(
        self, mock_subprocess
    ):
        """
        Test assertion error raised if variant counts pre- and post-bcftools
        annotate do not match
        """

        # This mocks the output of each of the subprocess calls in turn
        # within the bcftools_remove_csq_annotation() function
        # As the number of variants counted pre- and post- are different,
        # should raise an AssertionError
        mock_subprocess.side_effect = [
            Mock(stdout=b'12'), Mock(returncode=0), Mock(stdout=b'10')
        ]

        with pytest.raises(AssertionError):
            vcf.bcftools_remove_csq_annotation(
                'filter_vcf', 'csq_fields_to_drop'
            )
