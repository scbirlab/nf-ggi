#!/usr/bin/env bash

get_job_status () (
    local job_id="\$1"
    curl -i 'https://rest.uniprot.org/idmapping/status/'"\$job_id" \
        | tail -n1 \
        | jq -r '.["jobStatus"]'
)

set -e
set -x

if [[ "${filename}" == *.csv ]]
then 
    sep=,
else
    sep=\$'\\t'
fi

awk -F"\$sep" \
    -v column="${column}" \
    '
    NR == 1 { for (i = 1; i <= NF; i++) if ( \$i == column ) column_number = i } 
    NR > 1 { print \$column_number }
    ' \
    "${filename}" \
| sort -u \
| split -l 25 - name-chunk_

chunks=( name-chunk_* )
for chunk in "\${chunks[@]}"
do

    ids=\$(awk -v ORS="," '1' "\$chunk")
    job_id=\$(
        curl \
        --request POST 'https://rest.uniprot.org/idmapping/run' \
        --form 'ids="'"\$ids"'"' \
        --form 'from="${format}"'\
        --form 'to="UniProtKB"' \
        --form 'taxId="${organism_id}"' \
        | jq -r '.["jobId"]'
    )
    sleep 1
    job_status=\$(get_job_status "\$job_id")
    while [ "\$job_status" != "FINISHED" ]
    do
        if [  "\$job_status" == "ERROR" ]
        then
        echo "Got an error from UniProt for these inputs:"
        echo " - ${format}: ""\$ids"
        echo " - https://rest.uniprot.org/idmapping/status/""\$job_id""
        echo " - https://rest.uniprot.org/idmapping/stream/""\$job_id""
        exit 1
        fi
        sleep 1
        echo "Checking status of \$job_id..."
        job_status=\$(get_job_status "\$job_id")
    done

    echo "Job ready: \$job_id"
    curl -s "https://rest.uniprot.org/idmapping/stream/""\$job_id" \
    | jq -r '.["results"][]["to"]' \
    >> "uniprot-ids0.txt"

done
cat "uniprot-ids0.txt" | sort -u > "uniprot-ids.txt"
yield=\$(wc -l < "uniprot-ids.txt")
if [ "\$yield" -eq 0 ]
then
    echo "Did not find any UniProtKB IDs for filename ${filename}, column  ${column}! Example inputs:"
    echo "\$(cat "\${chunks[0]}")"
    echo "Check a lookup status here: "
    echo " - ${format}: ""\$ids"
    echo " - https://rest.uniprot.org/idmapping/status/""\$job_id""
    echo " - https://rest.uniprot.org/idmapping/stream/""\$job_id""
    exit 1
fi