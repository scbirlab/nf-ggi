#!/usr/bin/env bash

set -x
set -e

script_dir="$(dirname $0)"
cd "$script_dir"  # move to tests directory

for d in "$script_dir"/*/
do
    nextflow run "$script_dir"/.. \
        -profile gh -stub \
        -c "$d"/nextflow.config \
        --test \
        --inputs "$d"/inputs \
        --outputs "$d"/outputs
done