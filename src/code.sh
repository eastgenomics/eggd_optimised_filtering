#!/bin/bash

set -exo pipefail

main() {
    export BCFTOOLS_PLUGINS=/usr/local/libexec/bcftools/
    export PATH=$PATH:/home/dnanexus/.local/bin  # pip installs some packages here, add to path

    mark-section "Downloading inputs"
    dx-download-all-inputs --parallel


    mark-section "Installing packages"
    sudo -H python3 -m pip install --no-index --no-deps packages/*

    mark-section "Generating VCF with optimised filtering annotations"

    /usr/bin/time -v python3 optimised_filtering/add_optimised_filtering.py --input_vcf $input_vcf_path --config $config_path --panel_string $panel_string_path --genepanels $genepanels_path --panel_dump $panel_dump_path

    final_vcf=$(find . -name "*.G2P.vcf.gz")

    mark-section "Uploading output"
    output_vcf=$(dx upload $final_vcf --brief)
    dx-jobutil-add-output output_vcf "$output_vcf" --class=file

    mark-success
}