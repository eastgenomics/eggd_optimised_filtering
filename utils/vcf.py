import subprocess
import re

from collections import defaultdict, Counter
from pathlib import Path
from pysam import VariantFile

# Get path of parent directory
ROOT_DIR = Path(__file__).absolute().parents[1]

FIELDS_TO_SPLIT = 'SYMBOL,Consequence,gnomADe_AF,gnomADg_AF,TWE_AF,ClinVar_CLNSIG,ClinVar_CLNSIGCONF,SpliceAI_pred_DS_AG,SpliceAI_pred_DS_AL,SpliceAI_pred_DS_DG,SpliceAI_pred_DS_DL'

FIELDS_TO_COLLAPSE = ",".join([
    f"INFO/CSQ_{field}" for field in FIELDS_TO_SPLIT.split(',')
])

# class VCF():
#     def __init__(self, args) -> None:
#         self.args = args
#         self.sample_name = None


#     def process(self) -> None:


def bgzip(file) -> None:
    """
    Call bgzip on a given file

    Parameters
    ----------
    file : file to compress

    Outputs
    -------
    compressed file

    Raises
    ------
    AssertionError
        Raised when non-zero exit code returned by bgzip
    """
    print(f"Calling bgzip on {file}")

    output = subprocess.run(
        f"bgzip --force {file} -c > {file}.gz",
        shell=True,
        capture_output=True
    )

    assert output.returncode == 0, (
        f"\n\tError in compressing file with bgzip. File: {file}"
        f"\n\tExitcode:{output.returncode}"
        f"\n\t{output.stderr.decode()}"
    )


def bcftools_pre_process(input_vcf) -> str:
    """
    Decompose multiple transcript annotations to individual records, and split VEP CSQ string fields for SYMBOL, Consequence and gnomADe_AF to individual INFO keys. Adds a 'CSQ_' prefix to these fields extracted from the CSQ string to stop potential conflicts with existing INFO fields.

    Parameters
    ----------
    input_vcf : file
        path to VCF file to be split

    Returns
    -------
    output_vcf : str
        name of VCF split by bcftools
    """

    print(
        f"Splitting necessary fields from {input_vcf} "
        "with bcftools +split-vep"
    )
    output_vcf = f"{Path(input_vcf).stem}.split.vcf"

    # Check total rows before splitting out columns
    pre_split_total = subprocess.run(
        f"zgrep -v '^#' {input_vcf} | wc -l",
        shell=True,
        capture_output=True
    )
    # Split out the SYMBOL, Consequence, gnomADe_AF fields from the CSQ
    # string and name them with 'CSQ_{field}' as separate INFO fields
    cmd = (
        f"bcftools +split-vep -c {FIELDS_TO_SPLIT} -Ou -p 'CSQ_'"
        f" -d {input_vcf} -o {output_vcf}"
    )

    output = subprocess.run(cmd, shell=True, capture_output=True)

    assert output.returncode == 0, (
        f"\n\tError in splitting VCF with bcftools +split-vep. VCF: {input_vcf}"
        f"\n\tExitcode:{output.returncode}"
        f"\n\t{output.stderr.decode()}"
    )

    # Check total rows after splitting
    post_split_total = subprocess.run(
        f"zgrep -v '^#' {output_vcf} | wc -l",
        shell=True,
        capture_output=True
    )

    print(
        f"Total lines before splitting: {pre_split_total.stdout.decode()}"
        f"Total lines after splitting: {post_split_total.stdout.decode()}"
    )

    return output_vcf


def read_in_vcf(split_vcf_gz, flag_name):
    """
    Read in the VCF file with the pysam package

    Parameters
    ----------
    filename : string
        name of the annotated VCF file to be check against filtering parameters

    Returns
    -------
    vcf_contents: pysam VariantFile object
        contents of the VCF file as a pysam object
    """
    print(f"Reading in the split VCF {split_vcf_gz} with pysam")

    # Read in and create pysam object of the VCF
    vcf_contents = VariantFile(split_vcf_gz, 'r')

    # Get the name of the sample from the VCF
    sample_name = list(vcf_contents.header.samples)[0]

    # Add in our new flag to VCF header
    vcf_contents.header.info.add(flag_name, ".", "String", "Filtering flag")

    return vcf_contents, sample_name


def add_filtering_flag(
        sample_name, vcf_contents, panel_dict, rules, csq_types, flag_name
) -> dict:
    """
    Add the flags which will be used for filtering to each variant

    Parameters
    ----------
    sample_name : str
        name of the sample being tested (so we can obtain the genotypes)
    vcf_contents : pysam.VariantFile object
        pysam object containing all the VCF's info
    panel_dict : dict
        default dict with gene symbol as key and gene info as val
    rules : dict
        dict of the filtering rules for each of the inheritance types
    csq_types : list
        list of consequence types that we want to keep
    flag_name : str
        name of the info field to add for the flag

    Returns
    -------
    gene_variant_dict : dict
        dictionary of each gene with all variants in that gene as val
        and all filtering flags added
    """
    # Add each variant in a gene to a dict with gene as key and list
    # of variants as value
    gene_variant_dict = defaultdict(list)
    for record in vcf_contents:
        gene = record.info['CSQ_SYMBOL'][0]
        gene_variant_dict[gene].append(record)

    # Get the MOI for that gene
    for gene, variant_list in gene_variant_dict.items():
        gene_moi = panel_dict[gene]['mode_of_inheritance']
        af_threshold = rules[gene_moi]['af']

        variants_passing_af_csq_filters = []
        for variant in variant_list:
            consequences = variant.info['CSQ_Consequence'][0].split('&')
            exome_af = variant.info['CSQ_gnomADe_AF'][0]
            genome_af = variant.info['CSQ_gnomADg_AF'][0]
            twe_af = variant.info['CSQ_TWE_AF'][0]
            clinvar_clinsig = variant.info['CSQ_ClinVar_CLNSIG'][0]
            clinvar_conf = variant.info['CSQ_ClinVar_CLNSIGCONF'][0]
            splice_ag = variant.info['CSQ_SpliceAI_pred_DS_AG'][0]
            splice_al = variant.info['CSQ_SpliceAI_pred_DS_AL'][0]
            splice_dg = variant.info['CSQ_SpliceAI_pred_DS_DG'][0]
            splice_dl = variant.info['CSQ_SpliceAI_pred_DS_DL'][0]

            splice_list = [splice_ag, splice_al, splice_dg, splice_dl]
            splice_scores = [
                0.0 if score is None else int(score) for score in splice_list
            ]
            if not exome_af:
                exome_af = 0.0
            if not genome_af:
                genome_af = 0.0
            if not twe_af:
                twe_af = 0.0
            if not clinvar_conf:
                clinvar_conf = ''

            if (((exome_af < af_threshold) and any(x in consequences for x in csq_types) and (twe_af < 0.05)) or ((re.search(r'pathogenic', clinvar_conf, re.IGNORECASE) or any(splice_scores) >= 0.2))):
                variants_passing_af_csq_filters.append(variant)
            else:
                variant.info[flag_name] = 'EXCLUDE'

        # Get genotypes of all the variants in that gene which pass
        # consequence types and AF filters for that MOI
        gtypes = [
            '/'.join(
                [str(element) for element in variant.samples[sample_name]['GT']]
            ) for variant in variants_passing_af_csq_filters
        ]

        # Get number of hets and homs needed for the gene's MOI
        hets_needed = rules[gene_moi]['HET']
        homs_needed = rules[gene_moi]['HOM']

        # Count how many het and hom variants there are which pass filters
        # for that gene
        gt_counts = Counter(gtypes)
        het_count = gt_counts.get('1/0', 0) + gt_counts.get('0/1', 0)
        hom_count = gt_counts.get('1/1', 0)

        # If either enough hets or enough homs then add our PRIORITY flag
        if ((het_count >= hets_needed) or (hom_count >= homs_needed)):
            for variant in variants_passing_af_csq_filters:
                variant.info[flag_name] = 'PRIORITY'
        # If not enough hets or homs add EXCLUDE flag
        else:
            for variant in variants_passing_af_csq_filters:
                variant.info[flag_name] = 'EXCLUDE'

    return gene_variant_dict


# def rescue_whitelist_variants(gene_variant_dict):
#     for gene, variant_list in gene_variant_dict.items():
#         for variant in variant_list:
#             print(variant.start)
#             print(variant.end)
#             print(variant.ref)
#             print(variant.alts)



def write_out_flagged_vcf(flagged_vcf, gene_variant_dict, vcf_contents) -> None:
    """
    Write out each record to VCF using pysam

    Parameters
    ----------
    flagged_vcf : str
        Name of the VCF to be written out with flags added
    gene_variant_dict : dict
        dictionary with each gene and value as all of the variants in that gene
        as list
    vcf_contents : pysam.VariantFile object
        the contents of the VCF as a pysam object
    """
    print("Writing out variants to VCF")

    with VariantFile(flagged_vcf, 'w', header=vcf_contents.header) as out_vcf:
        # For each gene, write out each of the variants to the VCF with the flag
        for gene, variant_list in gene_variant_dict.items():
            for variant in variant_list:
                out_vcf.write(variant)


def bcftools_remove_csq_annotation(input_vcf):
    """
    Remove expanded CSQ strings which were used to check rules, as having them
    expanded would break eggd_generate_variant workbook and they are still
    present in the VEP CSQ string

    Parameters
    ----------
    input_vcf : str
        Name of VCF file to have CSQ annotations removed

    Returns
    -------
    output_vcf : str
        Name of output VCF file
    """
    output_vcf = f"{Path(input_vcf).stem}.G2P.vcf"

    # check total rows before splitting
    pre_split_total = subprocess.run(
        f"zgrep -v '^#' {input_vcf} | wc -l", shell=True,
        capture_output=True
    )

    cmd = (
        f"bcftools annotate -x {FIELDS_TO_COLLAPSE} {input_vcf} -o {output_vcf}"
    )

    output = subprocess.run(cmd, shell=True, capture_output=True)

    assert output.returncode == 0, (
        f"\n\tError in annotating VCF with bcftools. VCF: {input_vcf}"
        f"\n\tExitcode:{output.returncode}"
        f"\n\t{output.stderr.decode()}"
    )

    # check total rows after splitting
    post_split_total = subprocess.run(
        f"zgrep -v '^#' {output_vcf} | wc -l", shell=True,
        capture_output=True
    )

    print(
        f"Total lines before unsplitting: {pre_split_total.stdout.decode()}"
        f"Total lines after unsplitting: {post_split_total.stdout.decode()}"
    )

    return output_vcf


def add_annotation(flag_name, rules, csq_types, input_vcf, panel_dict):
    # Split and bgzip the annotated VCF
    split_vcf = f"{Path(input_vcf).stem}.split.vcf"
    flagged_vcf = f"{Path(input_vcf).stem}.flagged.vcf"

    bcftools_pre_process(input_vcf)

    # Read in the split (and bgzipped) VCF with pysam
    vcf_contents, sample_name = read_in_vcf(split_vcf, flag_name)
    gene_var_dict = add_filtering_flag(
        sample_name, vcf_contents, panel_dict, rules, csq_types, flag_name
    )
    #rescue_whitelist_variants(gene_var_dict)
    write_out_flagged_vcf(flagged_vcf, gene_var_dict, vcf_contents)
    final_vcf = bcftools_remove_csq_annotation(flagged_vcf)
    bgzip(final_vcf)
