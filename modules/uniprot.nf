process fetch_fastas_from_organism_id {

   tag "${id}"

   publishDir( 
      "${params.outputs}/sequences", 
      mode: 'copy',
      saveAs: { "${id}-${organism_id}.fasta.gz" }
   )

   input:
   tuple val( id ), val( organism_id )

   output:
   tuple val( id ), path( "proteome.fasta.gz" )

   script:
   """
   function get_proteome_id() {
      curl -s "https://rest.uniprot.org/proteomes/search?query=(taxonomy_id:${organism_id})&format=json" \
      | jq '.results[] \
      | select(.proteomeType == "'"\$1"' proteome").id'
   }
   QUERIES=("Reference and representative" "Reference" "Representative" "Other")
   PROTEOME_ID=
   for q in "\${QUERIES[@]}"
   do
      PROTEOME_ID=\$(get_proteome_id "\$q")
      if [ ! -z \$PROTEOME_ID ]
      then 
         break
      fi
   done

   wget "https://rest.uniprot.org/uniprotkb/stream?query=(proteome:\$PROTEOME_ID)&format=fasta&download=true&compressed=true" \
      -O proteome.fasta.gz \
      || (
         echo "Failed to download taxonomy ID ${organism_id} with proteome ID \$PROTEOME_ID from UniProt"
         exit 1   
      )
   """

}


process fetch_fasta_from_uniprot_id {

   tag "${uniprot_id}"

   publishDir( 
      "${params.outputs}/sequences", 
      mode: 'copy',
      saveAs: { "${id}-${uniprot_id}.fasta" },
   )

   input:
   tuple val( id ), val( uniprot_id )

   output:
   tuple val( id ), path( "protein.fasta" )

   script:
   """
   wget "https://rest.uniprot.org/uniprotkb/stream?query=(accession:${uniprot_id})&format=fasta&download=true&compressed=false" \
      -O protein.fasta \
   || (
      echo "Failed to download Uniprot ID ${uniprot_id} from UniProt"
      exit 1   
   )
   """

}


process fetch_fastas_from_uniprot_ids {

   tag "${uniprot_ids[0]}...${uniprot_ids[-1]}"

   publishDir( 
      "${params.outputs}/sequences", 
      mode: 'copy',
      saveAs: { "${uniprot_ids[0]}-${uniprot_ids[-1]}.fasta" },
   )

   input:
   tuple val( id ), val( uniprot_ids )

   output:
   tuple val( id ), path( "proteins.fasta" )

   script:
   """
   set -x
   curl -X GET --header 'Accept:text/x-fasta' \
      'https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=-1&accession=${uniprot_ids.join(',')}' \
   > proteins.fasta

   """

}



process map_uniprot_ids_from_file {

   tag "${id}-${column}:${from_type}"

   publishDir( 
      "${params.outputs}/uniprot_map", 
      mode: 'copy',
      saveAs: { "${id}-${column}.txt" },
   )

   input:
   tuple val( id ), val( organism_id ), path( filename ), val( from_type ), val( column )

   output:
   tuple val( id ), path( "uniprot-ids.txt" )

   script:
   """
   bash ${projectDir}/bin/uniprot/uniprot-ids.sh \
      "${filename}" "${column}" \
      "uniprot-ids0.txt" \
      UniProtKB ${from_type} \
      "${organism_id}"
   cut -f2 -d, < "uniprot-ids0.txt" \
   | tail -n+2 \
   | grep -v '^\$' \
   > "uniprot-ids.txt"
   """

}


process map_gene_names_from_file {

   tag "${id}-${column}:${out_column}"

   // if ( "${save}" == "true" ) {

   publishDir( 
      "${params.outputs}/ppi", 
      mode: 'copy',
      saveAs: { "${id}-${table.getSimpleName()}-${out_column}.tsv" },
   )

   // }

   input:
   tuple val( id ), path( table )
   val column
   val out_column
   val save

   output:
   tuple val( id ), path( "named-ids.tsv" )

   script:
   """
   set -x

   bash ${projectDir}/bin/uniprot/uniprot-ids.sh \
      "${table}" "${column}" \
      "uniprot-ids.csv" \
      Gene_Name UniProtKB

   python -c '
   import pandas as pd
   
   pd.merge(
      pd.read_csv("${table}", sep="\\t"),
      pd.read_csv("uniprot-ids.csv", sep=",").rename(columns={"Gene_Name": "${out_column}"}),
   ).drop_duplicates().to_csv("named-ids.tsv", sep="\\t", index=False)
   
   '

   """

   stub:
   """
   set -x

   cp ${projectDir}/data/559292-dca-stub.tsv input.tsv

   bash ${projectDir}/bin/uniprot/uniprot-ids.sh \
      input.tsv "${column}" \
      "uniprot-ids.csv" \
      Gene_Name UniProtKB
   
   python -c '
   import pandas as pd
   
   pd.merge(
      pd.read_csv("input.tsv", sep="\\t"),
      pd.read_csv("uniprot-ids.csv", sep=",").rename(columns={"Gene_Name": "${out_column}"}),
   ).drop_duplicates().to_csv("named-ids.tsv", sep="\\t", index=False)
   
   '

   """

}