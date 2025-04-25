#!/usr/bin/env bash

FILENAME="$1"
COLUMN="$2"
OUTPUT=${3:-"uniprot-ids.txt"}
TO=${4:-"UniProtKB"}
FROM=${5:-"Gene_Name"}
ORGANISM_ID=${6:-"Placeholder Organism ID"}

# Produces file called "uniprot-ids.txt"
# requires jq

tempfile="$(date)-uniprot-temp.csv"
uniprot_url='https://rest.uniprot.org/idmapping'
if [ "$FROM" == "UniProtKB" ]
then
    FROM="UniProtKB_AC-ID"
fi

get_job_status () (
    local job_id="$1"
    curl -i "$uniprot_url"/status/"$job_id" \
        | tail -n1 \
        | jq -r '.["jobStatus"]'
)

set -e
set -x

if [[ "$FILENAME" == *.csv ]]
then 
    sep=,
else
    sep=$'\\t'
fi

if [ "$(head -n1 "$FILENAME" | grep -c "$COLUMN")" -lt 1 ]
then
    echo "Column $COLUMN not in $FILENAME!"
    exit 1
fi

awk -F"$sep" \
    -v column="$COLUMN" \
    '
    NR == 1 { for (i = 1; i <= NF; i++) if ( $i == column ) column_number = i } 
    NR > 1 { print $column_number }
    ' \
    "$FILENAME" \
| sort -u \
| split -l 25 - name-chunk_

if [ "$FROM" != "UniProtKB_AC-ID" ]
then
    org_id_flag="--form taxId=""$ORGANISM_ID"
else
    org_id_flag=
fi

chunks=( name-chunk_* )
printf "$COLUMN","$TO"'\n' > "$tempfile"
for chunk in "${chunks[@]}"
do

    ids=$(awk -v ORS="," '1' "$chunk")
    curl_flags="--form ids=""$ids"" --form from=""$FROM"" --form to=""$TO"" $org_id_flag"
    job_id=$(
        curl \
        --request POST "$uniprot_url"/run \
        $curl_flags \
        | jq -r '.["jobId"]'
    )
    sleep 1
    if [  "$job_id" == "ERROR" ] ||  [ "$job_id" == "null" ]
    then
        echo "Got an error from UniProt for these inputs:"
        echo " - $FROM: $ids"
        echo " - $uniprot_url/status/$job_id"
        echo " - $uniprot_url/stream/$job_id"
        echo "Check: curl --request POST $uniprot_url/run $curl_flags"
        exit 1
    fi
    job_status=$(get_job_status "$job_id")
    while [ "$job_status" != "FINISHED" ]
    do
        if [  "$job_status" == "ERROR" ] ||  [ "$job_status" == "null" ]
        then
            echo "Got an error from UniProt for these inputs:"
            echo " - $FROM: $ids"
            echo " - $uniprot_url/status/$job_id"
            echo " - $uniprot_url/stream/$job_id"
            exit 1
        fi
        sleep 1
        echo "Checking status of $job_id..."
        job_status=$(get_job_status "$job_id")
    done

    echo "Job ready: $job_id"
    paste -d, "$chunk" \
        <(curl -s "$uniprot_url/stream/$job_id" \
            | jq -r '.["results"][]["to"]') \
    >> "$tempfile"

done

yield=$(wc -l < "$tempfile")
if [ "$yield" -eq 1 ]
then
    echo "Did not find any $TO IDs for filename $FILENAME, column  $COLUMN! Example inputs:"
    echo "$(cat "${chunks[0]}")"
    echo "Check a lookup status here: "
    echo " - $FROM: ""$ids"
    echo " - $uniprot_url/status/$job_id"
    echo " - $uniprot_url/stream/$job_id"
    exit 1
fi

grep -v '^\$' "$tempfile" > "$OUTPUT" 