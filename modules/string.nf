process fetch_string_database {

    tag "${organism_id}"
    label "some_mem"
    time '24h'

    publishDir( 
        "${params.outputs}/string", 
        mode: 'copy',
        saveAs: { "${organism_id}.string.tsv" }
    )

    input:
    tuple val( organism_id ), path( uniprot_ids )

    output:
    tuple val( organism_id ), path( "string.tsv" )

    script:
    """
    OUTFILE="string0.tsv"
    string_dl_url=https://stringdb-downloads.org/download/protein.links.full.v12.0
    string_db_url=https://string-db.org/api/tsv-no-header/get_string_ids
    printf 'string_id_1\\tstring_id_2\\tstring_cooccurence\\tstring_coexpression\\n' \
    > "\$OUTFILE"
    curl -s "\$string_dl_url/${organism_id}.protein.links.full.v12.0.txt.gz" \
    | zcat \
    | tail -n+2 \
    | tr ' ' \$'\\t' \
    | cut -f-2,6,8 \
    >> "\$OUTFILE" \
    || (
        echo "Failed to download taxonomy ID ${organism_id} from STRING <\$string_dl_url>"
        exit 1   
    )

    # Resolve UniProt IDs
    printf 'uniprot_id_1\\tstring_id_1\\n' \
    > "string-lookup.tsv"
    cat ${uniprot_ids} \
    | split -l 5 - uniprot-chunk_

    for chunk in uniprot-chunk_*
    do
        ids=\$(awk -v ORS="%0d" '1' "\$chunk")
        curl -s -X POST --data "species=${organism_id}&echo_query=1&identifiers=\$ids" \$string_db_url \
            | cut -f1,3 \
            >> "string-lookup.tsv" \
            || (
                echo "Failed to download UniProt lookup from STRING <\$string_db_url>"
                cat "\$chunk"
                exit 1  
            )
        sleep 0.1
    done

    python -c 'import pandas as pd; import sys; df = pd.read_csv("string-lookup.tsv", sep="\\t"); renamer = dict(uniprot_id_1="uniprot_id_2", string_id_1="string_id_2"); pd.read_csv("string0.tsv", sep="\\t").assign(organism_id="${organism_id}").merge(df).merge(df.rename(columns=renamer)).to_csv(sys.stdout, sep="\\t", index=False)' \
    > "string.tsv"
    """

}
