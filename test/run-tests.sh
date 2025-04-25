#!/usr/bin/env bash

set -x
set -e

script_dir="$(dirname $0)"
uniclust="/nemo/lab/johnsone/reference/hhdb/uniclust30/uniclust30_2018_08"
bfd="/nemo/lab/johnsone/reference/hhdb/bfd_metaclust/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"

# Examples without sample sheet
nextflow run "$script_dir"/.. \
    -stub -resume -profile gh \
    --test \
    --outputs custom-test \
    --organism_id 559292 \
    --filename "$script_dir"/custom/inputs/sgadata_costanzo2009_rawdata_101120-wheader_10q.txt \
    --format 'Gene_Name' \
    --column1 query_orf \
    --column2 array_orf \
    --dca --rf2t \
    --bfd "$bfd" --uniclust "$uniclust"

nextflow run "$script_dir"/.. \
    -stub -resume -profile gh \
    --test \
    --outputs bait-test \
    --organism_id 559292 \
    --bait P00931 \
    --dca --rf2t \
    --bfd "$bfd" --uniclust "$uniclust"

nextflow run "$script_dir"/.. \
    -stub -resume -profile gh \
    --test \
    --outputs bait-taxon-test \
    --organism_id 1773 \
    --bait 28369 \
    --bait_is_taxon \
    --interspecies \
    --dca --rf2t \
    --bfd "$bfd" --uniclust "$uniclust"

nextflow run "$script_dir"/.. \
    -stub -resume -profile gh \
    --test \
    --outputs self-test \
    --organism_id 243273 \
    --dca --rf2t \
    --bfd "$bfd" --uniclust "$uniclust"

# Examples with sample sheet
for d in "$script_dir"/{self,bait,bait-taxon,custom}/
do
    nextflow run "$script_dir"/.. \
        -stub -resume -profile gh \
        -c "$d"nextflow.config \
        --test \
        --plots \
        --sample_sheet "$d"inputs/sample-sheet.csv \
        --inputs "$d"inputs \
        --outputs "$d"outputs
done