# optimised_filtering
Optimised filtering is a tool used to add a flag to variants which:
- Pass standard filtering with bcftools
- Are within gnomAD AF thresholds based on the gene's MOI (from PanelApp)
- Fit the required zygosity counts (of those passing AF thresholds based on the gene's MOI)
 - E.g. a biallelic gene requires at least 1 homozygous variant or at least 2 heterozygous variants

## Description
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
`python3 add_optimised_filtering.py -i <vcf> -c <config-file-id> -p <panel_string> -g <genepanels-file-id> -d <panelapp-dump-file-id>`