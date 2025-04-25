process fetch_rhea_database {

    tag "${rhea_url}"

    publishDir( 
        "${params.outputs}/metabolites", 
        mode: 'copy',
        saveAs: { "_${it}" }
    )

    input:
    val rhea_url

    output:
    tuple path( "rhea-reaction-smiles.tsv" ), path( "rhea2uniprot+trembl.tsv.gz" )

    script:
    """
    for f in rhea2uniprot.tsv rhea2uniprot_trembl.tsv.gz rhea-reaction-smiles.tsv
    do
        wget ${rhea_url}/tsv/\$f
    done

    tail -n+2 -q "rhea2uniprot.tsv" <(zcat rhea2uniprot_trembl.tsv.gz) \
    | sort -k4 \
    | gzip --best \
    > "rhea2uniprot+trembl.tsv.gz"
    """
}

process match_uniprot_to_reactants {

    tag "${id}"

    publishDir( 
        "${params.outputs}/metabolites", 
        mode: 'copy',
        saveAs: { "${id}.rxn-smiles.tsv" }
    )

    input:
    tuple val( id ), path( uniprot_ids ), path( rhea_smiles ), path( rhea2uniprot )

    output:
    tuple val( id ), path( "rxn-smiles.tsv" )

    script:
    """
    OUTFILE=rxn-smiles.tsv
    zcat "${rhea2uniprot}" \
    | grep -wF -f "${uniprot_ids}" \
    | awk -v OFS='\\t' '\$2 == "UN" { \$1++; print \$0; \$1++; print \$0 }; \$2 != "UN"' \
    | sort -k1 \
    | join -t\$'\\t' -j 1 - <(sort -k1 "${rhea_smiles}") \
    > rxn-smiles0.tsv

    printf 'organism_id\\trhea_reaction_id\\tdirection\\trhea_master_reaction_id\\tuniprot_id\\treaction_smiles\\treactants\\tproducts\\n' \
    | cat - \
        <(paste rxn-smiles0.tsv \
            <(cat rxn-smiles0.tsv \
                | cut -f5 \
                | awk -F '>>' -v OFS='\\t' '{ print \$1,\$2 }') \
            | awk -v OFS='\\t' -v id="${id}" '{ print id,\$0 }') \
    > \$OUTFILE

    if [ \$(cat \$OUTFILE | wc -l) -gt 1 ]
    then 
        exit 0
    else
        >&2 echo "No entries in reaction file: \$OUTFILE."
        exit 1
    fi
    """
}