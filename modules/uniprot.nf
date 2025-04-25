process fetch_fastas_from_organism_id {

   tag "${organism_id}"

   publishDir( 
      "${params.outputs}/sequences", 
      mode: 'copy',
      saveAs: { "${organism_id}.fasta.gz" }
   )

   input:
   tuple val( id ), val( organism_id )

   output:
   tuple val( id ), path( "proteome.fasta.gz" )

   // TODO: filter only for representative proteome if no reference proteome
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
      saveAs: { "${uniprot_id}.fasta" },
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


process map_uniprot_ids_from_file {

   tag "${id}-${column}"

   publishDir( 
      "${params.outputs}/uniprot_map", 
      mode: 'copy',
      saveAs: { "${id}-${column}.txt" },
   )

   input:
   tuple val( id ), val( organism_id ), path( filename ), val( format ), val( column )

   output:
   tuple val( id ), path( "uniprot-ids.txt" )

   script:
   """
   get_job_status () (
      local job_id="\$1"
      curl -i 'https://rest.uniprot.org/idmapping/status/'"\$job_id" \
         | tail -n1 \
         | jq -r '.["jobStatus"]'
   )
   set -x
   
   if [[ "${filename}" == *.csv ]]
   then 
      sep=,
   else
      sep=\$'\\t'
   fi

   awk -F"\$sep" \
      -v column="${column}" \
      'NR == 1 { for (i = 1; i <= NF; i++) if ( \$i == column ) column_number = i } NR > 1 { print \$column_number }' \
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
   """

}


process map_gene_names_from_file {

   tag "${id}-${column}"

   publishDir( 
      "${params.outputs}/ppi", 
      mode: 'copy',
      saveAs: { "${id}-${column}.txt" },
   )

   input:
   tuple val( id ), path( filename ), val( column )

   output:
   tuple val( id ), path( "uniprot-ids.txt" )

   script:
   """
   get_job_status () (
      local job_id="\$1"
      curl -i 'https://rest.uniprot.org/idmapping/status/'"\$job_id" \
         | tail -n1 \
         | jq -r '.["jobStatus"]'
   )
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
   """

}