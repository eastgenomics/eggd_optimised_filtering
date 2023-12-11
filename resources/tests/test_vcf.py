import os
import sys
import unittest

from unittest.mock import Mock, MagicMock, patch

sys.path.append(os.path.abspath(
    os.path.join(os.path.realpath(__file__), '../../')
))

from tests import TEST_DATA_DIR
from utils import vcf


TEST_ANNOTATED_VCF = (
    "126560840-23326Q0015-23NGWES4-9526-F-103698_markdup_recalibrated_"
    "Haplotyper_annotated.vcf.gz"
)

# class TestBgzip(unittest.TestCase):
#     """
#     Test the function which uses subprocess to bgzip a file
#     """
#     @patch('vcf.subprocess.run')
#     def test_bgzip_called_correctly(self, mock_vcf):
#         """
#         Test
#         """
#         # process_mock = Mock()
#         # attrs = {'communicate.return_value'}
#         # mock_stdout = MagicMock()
#         # mock_stdout.configure_mock
#         # mocker.patch.object(CheckInputs, "__init__", return_value=None)
#         mock_response = Mock()
#         mock_response.return_code = 

#         mock_vcf.return_value = mock_response

#         vcf.bgzip(mock_vcf)
#         mock_vcf.assert_called_with()


class TestBcftoolsPreProcess():
    pass


class TestBcftoolsFilter():
    pass

class TestReadInVCF():
    """
    Test the VCF is read in with pysam correctly
    """
    annotated_control_vcf = os.path.join(TEST_DATA_DIR, TEST_ANNOTATED_VCF)
    def test_read_in_vcf(self):
        """
        Check VCF read in to pysam correctly
        """
        vcf_contents, sample_name, csq_fields_to_collapse = vcf.read_in_vcf(
            self.annotated_control_vcf
        )

        # Check sample name parsed correctly
        assert sample_name == '126560840-23326Q0015-23NGWES4-9526-F-103698'
        # Check the VEP fields are parsed and converted to fields to be
        # collapsed later correctly
        assert csq_fields_to_collapse == (
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

        assert vcf_contents.is_vcf

        # Get all of the pysam header records
        vcf_header_items = [
            record.items() for record in vcf_contents.header.records
        ]

        assert [
            ('ID', 'MOI'), ('Number', '.'), ('Type', 'String'),
            ('Description', '"Mode of inheritance from PanelApp (simplified)"'),
            ('IDX', '68')
        ] in vcf_header_items, "MOI not added to header correctly"


class TestAddFilteringFlag():
    pass


class TestWriteOutFlaggedVCF():
    pass


class TestBcftoolsRemoveCSQAnnotation():
    pass
