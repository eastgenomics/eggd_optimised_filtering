# eggd_optimised_filtering
This app takes an annotated VCF (output of eggd_VEP) and adds MOI info to VCF and filters according to config-provided bcftools string.

## Usage

To run the app:

```
dx run app-GZ9FZ78457v7qjBXPXqGByyP \
    -iinput_vcf=[annotated vcf] \
    -iconfig=[config file] \
    -ipanel_string=[panel string] \
    -igenepanels=[genepanels tsv] \
    -ipanel_dump=[panelapp dump json] \
    --destination=/path/to/output/dir -y

# example with WES vcf
dx run app-GZ9FZ78457v7qjBXPXqGByyP \
    -iinput_vcf=file-GVyyBg844vXGvyY77k9qGVyY \
    -iconfig=file-GZ9FfYj45B5ZgJQGBppGP8QZ \
    -ipanel_string="R149.1_Severe early-onset obesity_P" \
    -igenepanels=file-GY4QyKj4p65jx1xJqZKXBV79 \
    -ipanel_dump=file-GY4QxJ04p65zJf3937y01XBP \
    --destination=/output/wes_vcf -y
```

# Local Tool README

## optimised_filtering
Optimised filtering is a tool used to add MOI info, and perform filtering with bcftools.

### Description
INFO fields will be added (named 'MOI') to display the PanelApp Mode of Inheritance, simplified to the following categories:
- BIALLELIC
- MONOALLELIC
- BOTH
- XLR (X-Linked Recessive)
- XLD (X-Linked Dominant)
- MITOCHONDRIAL
- OTHER
- UNKNOWN
- NONE

Optimised filtering uses:
- [bcftools](https://samtools.github.io/bcftools/bcftools.html, "bcftools website")
- [pysam](https://pysam.readthedocs.io/en/latest/, "pysam documentation")

### Requirements
- bcftools
    - bcftools +split-vep

A config JSON file is required, which is given as an argument to the tool as a DNAnexus file ID: Example config:

```JSON
{
	"bcftools_filter_string": "bcftools filter --soft-filter \"EXCLUDE\" -m + -e '(CSQ_Consequence~\"synonymous_variant\" | CSQ_Consequence~\"intron_variant\" | CSQ_Consequence~\"upstream_gene_variant\" | CSQ_Consequence~\"downstream_gene_variant\" | CSQ_Consequence~\"intergenic_variant\" | CSQ_Consequence~\"5_prime_UTR_variant\" | CSQ_Consequence~\"3_prime_UTR_variant\" | CSQ_gnomADe_AF>0.01 | CSQ_gnomADg_AF>0.01 | CSQ_TWE_AF>0.05) & CSQ_ClinVar_CLNSIGCONF!~ \"pathogenic\\/i\" & (CSQ_SpliceAI_pred_DS_AG<0.2 | CSQ_SpliceAI_pred_DS_AG==\".\") & (CSQ_SpliceAI_pred_DS_AL<0.2 | CSQ_SpliceAI_pred_DS_AL==\".\") & (CSQ_SpliceAI_pred_DS_DG<0.2 | CSQ_SpliceAI_pred_DS_DG==\".\") & (CSQ_SpliceAI_pred_DS_DL<0.2 | CSQ_SpliceAI_pred_DS_DL==\".\") | (MOI=\"BIALLELIC\" & (CSQ_gnomadg_AF>0.005 | CSQ_gnomade_AF>0.005))'"
}
```
