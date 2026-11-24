#!/usr/bin/env bash

set -exuo pipefail

GITHUB=${1:-no}

if [ "$GITHUB" == "gh" ]
then
    export NXF_CONTAINER_ENGINE=docker
    docker_flag='-profile gh -stub'
    uniclust="uniclust30_2018_08"
    bfd="bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
    echo $uniclust > "$uniclust"_test
    echo $uniclust > "$uniclust".test
    echo $bfd > "$bfd"_test
else
    export SINGULARITY_FAKEROOT=1
    # docker_flag='-profile local -with-singularity'
    docker_flag=
    uniclust="/nemo/lab/johnsone/reference/hhdb/uniclust30/uniclust30_2018_08"
    bfd="/nemo/lab/johnsone/reference/hhdb/bfd_metaclust/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
fi

script_dir="$(dirname $0)"

# Examples without sample sheet
nextflow run "$script_dir"/.. \
    -resume $docker_flag \
    --test \
    --outputs "$script_dir"/custom-inters-test \
    --interspecies \
    --organism_id 9606 \
    --organism_id2 353152 \
    --filename "$script_dir"/custom/inputs/interactions.csv \
    --column1 protein1 \
    --column2 protein2 \
    --dca --rf2t \
    --plots \
    --bfd "$bfd" --uniclust "$uniclust"

nextflow run "$script_dir"/.. \
    -resume $docker_flag \
    --test \
    --outputs "$script_dir"/bait-taxon-test \
    --organism_id 1773 \
    --bait 28369 \
    --bait_is_taxon \
    --interspecies \
    --dca --rf2t \
    --bfd "$bfd" --uniclust "$uniclust"
    
nextflow run "$script_dir"/.. \
    -resume $docker_flag \
    --test \
    --outputs "$script_dir"/custom-test \
    --organism_id 559292 \
    --filename "$script_dir"/custom/inputs/sgadata_costanzo2009_rawdata_101120-wheader_10q.txt \
    --format 'Gene_Name' \
    --column1 query_orf \
    --column2 array_orf \
    --dca --rf2t \
    --bfd "$bfd" --uniclust "$uniclust"

nextflow run "$script_dir"/.. \
    -resume $docker_flag \
    --test \
    --outputs "$script_dir"/bait-test \
    --organism_id 559292 \
    --bait P00931 \
    --dca --rf2t \
    --bfd "$bfd" --uniclust "$uniclust"

nextflow run "$script_dir"/.. \
    -resume $docker_flag \
    --test \
    --outputs "$script_dir"/self-test \
    --organism_id 243273 \
    --dca --rf2t \
    --bfd "$bfd" --uniclust "$uniclust"

# Examples with sample sheet
for d in "$script_dir"/{self,bait,bait-taxon,custom}/
do
    nextflow run "$script_dir"/.. \
        -resume $docker_flag \
        -c "$d"nextflow.config \
        --test \
        --plots \
        --sample_sheet "$d"inputs/sample-sheet.csv \
        --inputs "$d"inputs \
        --outputs "$d"outputs \
        --bfd "$bfd" --uniclust "$uniclust"
done