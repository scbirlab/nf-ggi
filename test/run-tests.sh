#!/usr/bin/env bash

set -x
set -e

script_dir="$(dirname $0)"

for d in "$script_dir"/*/
do
    nextflow run "$script_dir"/.. \
        -profile gh -stub \
        -c "$d"/nextflow.config \
        --test \
        --sample-sheet "$d"/inputs.sample-sheet.csv \
        --inputs "$d"/inputs \
        --outputs "$d"/outputs
done