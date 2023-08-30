# optimised_filtering
Optimised filtering is a tool used to add a flag to indicate variants which:
- Pass standard filtering with bcftools
- Do not exceed gnomAD AF thresholds based on the gene's MOI (from PanelApp)
- Fit the required zygosity counts (of those passing AF thresholds based on the gene's MOI)
    - E.g. a biallelic gene requires at least 1 homozygous variant or at least 2 heterozygous variants

## Description
The flag values which could be added to a variant:
- `PRIORITISED`: The variant passes filtering
- `NOT_PRIORITISED`: The variant does not pass filtering
- `NOT_ASSESSED`: The variant has not been assessed, either because:
    - The gene could not be found in the PanelApp dump
    - The gene's mode of inheritance could not be found
    - An AF threshold for the specific mode of inheritance could not be found in the config file

If the variant is not prioritised, the reason will be added to the `Filter_reason` INFO field.

Optimised filtering uses:
- [bcftools](https://samtools.github.io/bcftools/bcftools.html, "bcftools website")
- [pysam](https://pysam.readthedocs.io/en/latest/, "pysam documentation")

## Requirements
- bcftools
    - bcftools +split-vep

A config JSON file is required, which is given as an argument to the tool as a DNAnexus file ID: Example config:

```JSON
{
	"flag_name": "G2P",
	"filtering_rules": {
		"biallelic": {
			"af": 0.005,
			"HET": 2,
			"HOM": 1
		},
		"monoallelic": {
			"af": 0.0001,
			"HET": 1,
			"HOM": 1
		},
		"both_monoallelic_and_biallelic": {
			"af": 0.005,
			"HET": 1,
			"HOM": 1
		},
		"standard_filtering": {
			"af": 0.01,
			"HET": 1,
			"HOM": 1
		}
	},
	"bcftools_filter_string": "bcftools filter --soft-filter \"EXCLUDE\" -m + -e '(CSQ_Consequence~\"synonymous_variant\" | CSQ_Consequence~\"intron_variant\" | CSQ_Consequence~\"upstream_gene_variant\" | CSQ_Consequence~\"downstream_gene_variant\" | CSQ_Consequence~\"intergenic_variant\" | CSQ_Consequence~\"5_prime_UTR_variant\" | CSQ_Consequence~\"3_prime_UTR_variant\" | CSQ_gnomADe_AF>0.01 | CSQ_gnomADg_AF>0.01 | CSQ_TWE_AF>0.05) & CSQ_ClinVar_CLNSIGCONF!~ \"pathogenic\/i\" & (CSQ_SpliceAI_pred_DS_AG<0.2 | CSQ_SpliceAI_pred_DS_AG==\".\") & (CSQ_SpliceAI_pred_DS_AL<0.2 | CSQ_SpliceAI_pred_DS_AL==\".\") & (CSQ_SpliceAI_pred_DS_DG<0.2 | CSQ_SpliceAI_pred_DS_DG==\".\") & (CSQ_SpliceAI_pred_DS_DL<0.2 | CSQ_SpliceAI_pred_DS_DL==\".\")'"
}
```

## Usage
The tool takes a VCF which has been annotated with VEP and produces a VCF with the two additional INFO fields specified above. The tool can be run as follows:

```
usage: add_optimised_filtering.py {ARGUMENTS}

arguments:
 --input_vcf/-i [str]        Path to the input annotated VCF
 --config/-c [str]           DNAnexus file ID of config file containing filtering rules
 --panel_string/-p [str]     The panel(s) the patient is being tested for
 --genepanels/-g [str]       DNAnexus file ID of the genepanels file containing panel ID
 --panel_dump/-d [str]       DNAnexus file ID of a dump from PanelApp in JSON format
```


Example:
```
python3 add_optimised_filtering.py \
-i X123456_markdup_recalibrated_Haplotyper_annotated.vcf.gz \
-c project-GVqVPk04p65vjXb2kj6FqFKf:file-GYfv8V84p65bxXp1G8k5g44j \
-p 'R149.1_Severe early-onset obesity_P' \
-g project-GVqVPk04p65vjXb2kj6FqFKf:file-GY4QyKj4p65jx1xJqZKXBV79 \
-d project-GVqVPk04p65vjXb2kj6FqFKf:file-GY4QxJ04p65zJf3937y01XBP
```