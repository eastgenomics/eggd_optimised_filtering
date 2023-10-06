#!/bin/bash
# optimised_filtering
#
# Any code outside of main() (or any entry point you may add) is
# ALWAYS executed, followed by running the entry point itself.

# Exit at any point if there is any error and output each line as it is executed (for debugging)
set -e -x -o pipefail


main() {

    # install dependencies
    sudo -H python3 -m pip install --no-index --no-deps packages/*
    export BCFTOOLS_PLUGINS=/usr/local/libexec/bcftools

    # get inputs
    dx-download-all-inputs --parallel

    # set zygosity switch correctly for command below
    if [ $zygosity = true ]
        then zyg="-z"
    elif [ $zygosity = false ]
        then zyg=""

    # run tool
    python3 /add_optimised_filtering.py \
        -i $input_vcf_path \
        -c $config_path \
        -p "$panel_string" \
        -g $genepanels_path \
        -d $panel_dump_path \
        $zyg

    # prepare outputs
    echo "All scripts finished successfully, uploading output files to dx"
    if [ "$debug_fail_end" == 'true' ]; then exit 1; fi
    mkdir -p ~/out/result_files/
    mv *flagged*.vcf.gz ~/out/result_files/

    # Upload output files
    dx-upload-all-outputs --parallel

}