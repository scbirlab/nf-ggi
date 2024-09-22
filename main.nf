#!/usr/bin/env nextflow

/*
========================================================================================
   Gene-gene interaction predicting Nextflow Workflow
========================================================================================
   Github   : https://github.com/scbirlab/nf-ggi
   Contact  : Eachan Johnson <eachan.johnson@crick.ac.uk>
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl=2

/*
========================================================================================
   Help text
========================================================================================
*/
if ( params.help ) {
   println """\
         S C B I R   G E N E - G E N E   I N T E R A C T I O N   P R E D I C T I O N   P I P E L I N E
         =============================================================================================
         Nextflow pipeline to predict gene-gene interactions based on protein-protein interaction 
         predictions and similar metabolites.

         Usage:
            nextflow run scbirlab/nf-ggi --accession <acc number>
            nextflow run scbirlab/nf-ggi -c <config-file>

         Required parameters:
            sample_sheet         UniProt accession number for organism of interest.
            uniclust, bfd        Paths to get HHblits databases.

         Optional parameters (with defaults):  
            test           Whether to run in test mode. Default: false.
            non_self       Whether to run in non_self mode. Default: false.
            bactch_size    What size to batch protein-protein interactions into. Default: 100.
            rhea_url       URL to download Rhea reaction database. Default: "https://ftp.expasy.org/databases/rhea"
            outputs        Output folder. Default: "outputs".

         The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.
   
   """.stripIndent()
   exit 0
}

/*
========================================================================================
   Check parameters
========================================================================================
*/
if ( !params.sample_sheet ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a sample sheet including UniProt proteome ID for organism of interest.")
}
if ( !params.uniclust ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to UniClust database.")
}
if ( !params.bfd ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to BFD database.")
}

working_dir = params.outputs

sequences = "sequences"
ppi = "ppi"
metabi = "metabolites"
coex = "coexpression"

log.info """\
         S C B I R   G E N E - G E N E   I N T E R A C T I O N   P R E D I C T I O N   P I P E L I N E
         =============================================================================================
         test mode               : ${params.test}
         non-self mode           : ${params.non_self}
         inputs
            sample sheet         : ${params.sample_sheet}
            UniClust database    : ${params.uniclust}
            BFD                  : ${params.bfd}
            Rhea                 : ${params.rhea_url}
         batch size              : ${params.batch_size}
         output                  : ${params.outputs}
         """
         .stripIndent()


/*
========================================================================================
   Create Channels
========================================================================================
*/

database_ch = Channel.of( tuple( params.uniclust, params.bfd ) )
rhea_url_ch = Channel.of( params.rhea_url )
sample_sheet_ch = Channel.fromPath( params.sample_sheet,
                                    checkIfExists: true )


/*
========================================================================================
   MAIN Workflow
========================================================================================
*/

workflow {

   sample_sheet_ch
      .splitCsv( header: true )
      .map { params.non_self ? tuple( it.organism_id, it.proteome_name, it.bait ) : tuple( it.organism_id, it.proteome_name ) }
      .set { sample_sheet }
   sample_sheet
      .map { it[0] }
      .unique()
      .map { sample -> [ "sequences", "msa", "ppi", 
                         "coexpression", "metabolites" ].collect { "${params.outputs}/${sample}/${it}" } }
      .subscribe { it.each { println "Creating output directory: ${it}" } }
      .map { it.each { file( it ).mkdirs() } }
   // Get FASTA sequences by UniProt ID, then split each sequence into a 
   // single FASTA file.
   sample_sheet.map { it[0] } | PULL_FASTA_SEQUENCES
   PULL_FASTA_SEQUENCES.out
      .splitFasta( elem: 1, record: [id: true, text: true] )
      .set { fastas }

   // If in test mode, only take 3 proteins, otherwise take all. 
   // Then munge to extract IDs
   ( params.test ? fastas.take(3) : fastas )
      .map { tuple( it[0], it[1].id.split('\\|')[1], it[1].id.split('\\|')[2], it[1].text ) }  // Organism ID, UniProtID, Entry Name, FASTA text
      .set { sample_fastas }

   if ( params.non_self ) {
      sample_sheet.map { it[2] }      // Bait UniProtID
         | PULL_BAIT_FASTA_SEQUENCES  // Bait UniProtID, FASTA file
      sample_sheet
         .map { tuple( it[2], it[0] ) }  // Bait UniProtID, Organism ID
         .combine( PULL_BAIT_FASTA_SEQUENCES.out, 
                   by: 0 )  // Bait UniProtID, Organism ID, FASTA file
         .set { bait_fastas }
   }
   
   // Get reactants from Rhea database using UniProtIDs
   rhea_url_ch | DOWNLOAD_RHEA_DB
   DOWNLOAD_RHEA_DB.out.set { rhea_db }
   sample_fastas
      .collectFile( newLine: true ) { [ "${it[0]}.txt", "${it[1]}\t${it[2]}" ] }
      .map { tuple( it.getSimpleName(), it ) }  // Organism ID, UniProtID list file
      .combine( rhea_db ) // Organism ID, UniProtID list file, Rhea SMILES, Rhea2UniProt, Rhea2UniProt_Trembl
      | GET_ENZYME_REACTANTS  // Organism ID, reaction SMILES file
      | CONNECT_METABOLITES 

   // Get STRING connections
   sample_fastas
      .map { it[0..1] }  // Organism ID, UniProtID
      .collectFile( newLine: true ) { [ "${it[0]}.txt", it[1] ]}  // UniProtID list
      .map { tuple( it.getSimpleName(), it ) }  // Organism ID, UniProtID list
      | DOWNLOAD_STRING_DB

   // Get MSAs and cross them within each organism
   sample_fastas
      .collectFile( newLine: true ) { [ "${it[1]}.fasta", it[-1] ] }
      .map { tuple( it.getSimpleName(), it )  }  // UniProtID, FASTA file
      .set { fasta_files }
   sample_fastas
      .map { it[1..0] }  // UniProtID, Organism ID
      .join( fasta_files, 
            by: 0,
            failOnMismatch: true, 
            failOnDuplicate: true )  // UniProtID, Organism ID, FASTA file
      .set { pre_msa }

   if ( params.non_self ) {
      pre_msa.concat( bait_fastas ).set { pre_msa }
   }

   pre_msa
      .map { tuple( it[1], it[0], it[2] ) }  // Organism ID, UniProtID, FASTA file
      .combine( database_ch )  // Organism ID, UniProtID, FASTA file, uniclust, bfd
      | GET_MSA  // Organism ID, UniProtID, MSA file
   
   if ( params.non_self ) {
      sample_sheet_ch
         .map { it[2] } // Bait Uniprot ID
         .toList()      // [Bait Uniprot ID,...]
         .set { bait_list }
      GET_MSA.out
         .filter { bait_list.contains( it[1] ) }
         .set { counter_msas }
   } else {
      GET_MSA.out
         .set { counter_msas }
   }

   counter_msas
      .map { tuple( it[0], it[2] ) }  // // Organism ID, MSA file
      .groupTuple( by: 0, sort: { a, b -> a.getSimpleName() <=> b.getSimpleName() } )  // Organism ID, [MSA file, ...]
      .map { tuple ( it[0], 
                     it[1].withIndex().collect { el, i -> Math.round(Math.floor(i / params.batch_size)) },  // add batch index
                     it[1] ) }  // Organism ID, [batch_i, ...], [MSA file, ...]
      .transpose()  // Organism ID, batch_i, MSA file
      .groupTuple( by: [0,1], sort: { a, b -> a.getSimpleName() <=> b.getSimpleName() } )  // Organism ID, batch_i, [MSA file, ...]
      .set { msas_by_proteome }  
   GET_MSA.out
      .combine( msas_by_proteome, by: 0 )  // Organism ID, UniProtID, MSA file, batch_i, [MSA file, ...]
      .map { a -> a[0..-2] + [ 
         a[-1]
         .sort { b, c -> b.getSimpleName() <=> c.getSimpleName() }
         .takeWhile { it.getSimpleName() != a[2].getSimpleName() } 
      ] }  // take lower triangle per proteome
      .filter { it[-1].size() > 0 }  // filter out trivial (size-0) elements
      .set { msas_to_process }  // Organism ID, UniProtID, MSA file, batch_i, [MSA file, ...]

   // Calculate evolutionary coupling
   msas_to_process | DIRECT_COUPLING_ANALYSIS 

   if ( params.test ) {
      msas_to_process | RF2TRACK
      msas_to_process | ALPHAFOLD2
   }
   // ALPHAFOLD2.out
   //           .map { it[2] }
   //           .collect()

}

process PULL_FASTA_SEQUENCES {

   tag "${organism_id}"

   publishDir( "${params.outputs}/${organism_id}/sequences", 
               mode: 'copy' )

   input:
   val organism_id

   output:
   tuple val( organism_id ), path( "*.fasta.gz" )

   // TODO: filter only for representative proteome if no reference proteome
   script:
   """
   function get_proteome_id() {
      curl -s "https://rest.uniprot.org/proteomes/search?query=(taxonomy_id:${organism_id})&format=json" | jq '.results[] | select(.proteomeType == "'"\$1"' proteome").id'
   }
   QUERIES=("Reference and representative" "Representative" "Other")
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
      -O ${organism_id}.fasta.gz \
      || (
         echo "Failed to download taxonomy ID ${organism_id} with proteome ID \$PROTEOME_ID from UniProt"
         exit 1   
      )
   """

}


process PULL_BAIT_FASTA_SEQUENCES {

   tag "${uniprot_id}"

   publishDir( "${params.outputs}", 
               mode: 'copy' )

   input:
   val uniprot_id

   output:
   tuple val( uniprot_id ), path( "*.fasta" )

   script:
   """
   wget "https://rest.uniprot.org/uniprotkb/stream?query=(accession:${uniprot_id})&format=fasta&download=true&compressed=false" \
      -O ${uniprot_id}.fasta \
      || (
         echo "Failed to download Uniprot ID ${uniprot_id} from UniProt"
         exit 1   
      )
   """

}


process DOWNLOAD_STRING_DB {

   tag "${organism_id}"

   publishDir( "${params.outputs}/${organism_id}/coexpression", 
               mode: 'copy' )

   input:
   tuple val( organism_id ), path( uniprot_ids )

   output:
   tuple val( organism_id ), path( "*.string.tsv" )

   script:
   """
   OUTFILE=${organism_id}.string0.tsv
   printf 'string_id_1\\tstring_id_2\\tstring_cooccurence\\tstring_coexpression\\n' > \$OUTFILE
   curl -s "https://stringdb-downloads.org/download/protein.links.full.v12.0/${organism_id}.protein.links.full.v12.0.txt.gz" \
      | zcat \
      | tail -n+2 \
      | tr ' ' \$'\\t' \
      | cut -f-2,6,8 \
      >> \$OUTFILE \
      || (
         echo "Failed to download taxonomy ID ${organism_id} from STRING <https://stringdb-downloads.org>"
         exit 1   
      )

   # Resolve UniProt IDs
   printf 'uniprot_id_1\\tstring_id_1\\n' > ${organism_id}-string-lookup.tsv
   cat ${uniprot_ids} | split -l 5 - uniprot-chunk_
   for chunk in uniprot-chunk_*
   do
      
      curl -s -X POST --data "species=${organism_id}&echo_query=1&identifiers="\$(awk -v ORS="%0d" '1' \$chunk) https://string-db.org/api/tsv-no-header/get_string_ids \
         | cut -f1,3 \
         >> ${organism_id}-string-lookup.tsv \
         || (
            echo "Failed to download UniProt lookup from STRING <https://string-db.org/api/tsv-no-header/get_string_ids>"
            cat \$chunk
            exit 1  
         )
      sleep 0.1
   done

   python -c 'import pandas as pd; import sys; lookup = pd.read_csv("${organism_id}-string-lookup.tsv", sep="\\t"); renamer = dict(uniprot_id_1="uniprot_id_2", string_id_1="string_id_2"); pd.read_csv("${organism_id}.string0.tsv", sep="\\t").assign(organism_id="${organism_id}").merge(lookup).merge(lookup.rename(columns=renamer)).to_csv(sys.stdout, sep="\\t", index=False)' \
      > ${organism_id}.string.tsv
   """

}


process DOWNLOAD_RHEA_DB {

   tag "${rhea_url}"

   input:
   val rhea_url

   output:
   tuple path( "rhea-reaction-smiles.tsv" ), path( "rhea2uniprot.tsv" ), path( "rhea2uniprot_trembl.tsv.gz" )

   script:
   """
   wget ${rhea_url}/tsv/rhea2uniprot.tsv
   wget ${rhea_url}/tsv/rhea2uniprot_trembl.tsv.gz
   wget ${rhea_url}/tsv/rhea-reaction-smiles.tsv
   """
}


process GET_ENZYME_REACTANTS {

   tag "${organism_id}"

   publishDir( "${params.outputs}/${organism_id}/metabolites", 
               mode: 'copy' )

   input:
   tuple val( organism_id ), path( uniprot_ids ), path( rhea_smiles ), path( rhea2uniprot ), path( rhea2uniprot_tr )

   output:
   tuple val( organism_id ), path( "*.rxn-smiles.tsv" )

   script:
   """
   OUTFILE=${organism_id}.rxn-smiles.tsv
   join -t\$'\\t' -1 4 -2 1 <(tail -n+2 -q "${rhea2uniprot}" <(zcat ${rhea2uniprot_tr}) | sort -k4) <(sort -k1 "${uniprot_ids}") \
      | awk -v OFS='\\t' '\$3=="UN"{ \$2++; print \$0; \$2++; print \$0 }; \$3!="UN"' \
      | sort -k2 \
      | join -t\$'\\t' -1 2 -2 1 - <(sort -k1 "${rhea_smiles}") \
      > rxn-smiles0.tsv

   cat <(printf 'organism_id\\trhea_reaction_id\\tuniprot_id\\tdirection\\trhea_master_reaction_id\\tentry_name\\treaction_smiles\\treactants\\tproducts\\n') \
      <(paste rxn-smiles0.tsv \
         <(cat rxn-smiles0.tsv | cut -f6 | awk -F '>>' -v OFS=\$'\t' '{ print \$1,\$2 }') \
         | awk -v OFS='\\t' '{ print "${organism_id}",\$0 }') \
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


process CONNECT_METABOLITES {

   tag "${organism_id}"
   label "big_time"

   publishDir( "${params.outputs}/${organism_id}/metabolites", 
               mode: 'copy' )

   input:
   tuple val( organism_id ), path( reaction_table )

   output:
   tuple val( organism_id ), path( "*.metabolism-connection.tsv" )

   script:
   """
   cat ${reaction_table} | python ${projectDir}/bin/metabolism/connect-metabolites.py > ${organism_id}.metabolism-connection.tsv
   """
}


process GET_MSA {

   tag "${organism_id} : ${uniprot_id}"
   label 'big_cpu'
   errorStrategy 'retry'
   maxRetries 2

   publishDir( "${params.outputs}/${organism_id}/msa", 
               mode: 'copy' )

   // Proteome ID, UniProtID, FASTA file, uniclust, bfd
   input:
   tuple val( organism_id ), val( uniprot_id ), file( fasta ), val( uniclust ), val( bfd )

   output:
   tuple val( organism_id ), val( uniprot_id ), path( "${organism_id}-${uniprot_id}.a3m" )

   script:
   """
   set -x
   dbs=(${uniclust} ${bfd})
   for d in \${dbs[@]}
   do
   hhblits \
      -cpu ${task.cpus} \
      -maxmem ${task.memory.getGiga()} \
      -v 2 \
      -i "${fasta}" \
      -d \$d \
      -e 0.001 \
      -o /dev/null \
      -oa3m "\$(basename \$d).a3m" \
      -cov 60 \
      -n 3 \
      -realign -realign_max 10000
   done

   outputs=( *.a3m )
   cat <(head -n2 \${outputs[0]}) <(tail -n+3 -q \${outputs[@]}) \
      > "${organism_id}-${uniprot_id}.a3m"
   for f in \${outputs[@]}
   do
      if [ \$f != "${organism_id}-${uniprot_id}.a3m" ]
      then
         rm \$f
      fi
   done
   """

   stub:
   """
   head -n2 ${fasta} > "${organism_id}-${uniprot_id}.a3m"
   """
}


process DIRECT_COUPLING_ANALYSIS {

   label 'big_time'
   tag "${uniprot_id}-batch_${batch_idx}"
   stageInMode 'link'
   errorStrategy 'retry'
   maxRetries 2

   publishDir( "${params.outputs}/${organism_id}/ppi", 
               mode: 'copy' )

   // Proteome ID, UniProtID, MSA file, [MSA file, ...]
   input:
   tuple val( organism_id ), val( uniprot_id ), path( msa1, stageAs: "ref/ref.a3m" ), val( batch_idx ), path( msa2 )

   output:
   tuple val( organism_id ), val( uniprot_id ), path( "*-batch*_dca.tsv" ), emit: main
   path "*-batch*_dca" , emit: plots

   script:
   """
   MSA_LIST=msa-list.txt
   for item in *.a3m
   do
      echo "\$item" >> \$MSA_LIST
   done
   sppid dca-single \
      <(echo "ref/ref.a3m") \
      --msa2 \$MSA_LIST \
      --list-file \
      --apc \
      --output "${organism_id}-${uniprot_id}-batch_${batch_idx}_dca.tsv" \
      --plot "${organism_id}-${uniprot_id}-batch_${batch_idx}_dca"
   """
}

process RF2TRACK {

   label 'big_time'
   tag "${uniprot_id}-batch_${batch_idx}"
   stageInMode 'link'
   errorStrategy 'retry'
   maxRetries 2

   publishDir( "${params.outputs}/${organism_id}/ppi", 
               mode: 'copy' )

   // Proteome ID, UniProtID, MSA file, [MSA file, ...]
   input:
   tuple val( organism_id ), val( uniprot_id ), path( msa1, stageAs: "ref/ref.a3m" ), val( batch_idx ), path( msa2 )

   output:
   tuple val( organism_id ), val( uniprot_id ), path( "*-batch_*.tsv" ), emit: main
   path "*-batch*_rf2t", emit: plots

   script:
   """
   MSA_LIST=msa-list.txt
   for item in *.a3m
   do
      echo "\$item" >> \$MSA_LIST
   done
   PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True sppid rf2t-single \
      <(echo "ref/ref.a3m") \
      --msa2 \$MSA_LIST \
      --list-file \
      --output "${organism_id}-${uniprot_id}-batch_${batch_idx}_rf2t.tsv" \
      --plot "${organism_id}-${uniprot_id}-batch_${batch_idx}_rf2t"
   """

   stub:
   """
   mkdir "${organism_id}-${uniprot_id}-batch_${batch_idx}_rf2t.tsv"
   touch "${organism_id}-${uniprot_id}-batch_${batch_idx}_af2/plot.png"
   echo "Skipping RF2t for stub"
   """
}


process ALPHAFOLD2 {

   label 'gpu'
   tag "${uniprot_id}-batch_${batch_idx}"
   stageInMode 'link'
   errorStrategy 'retry'
   maxRetries 2

   publishDir( "${params.outputs}/${organism_id}/ppi", 
               mode: 'copy' )

   input:
   tuple val( organism_id ), val( uniprot_id ), path( msa1, stageAs: "ref/ref.a3m" ), val( batch_idx ), path( msa2 )

   output:
   tuple val( organism_id ), val( uniprot_id ), path( "*-batch_*_af2" ), emit: main
   path "*-batch_*_af2/*.pdb", emit: pdb

   script:
   """
   MSA_LIST=msa-list.txt
   for item in *.a3m
   do
      echo "\$item" >> \$MSA_LIST
   done
   export CUDNN_PATH=\$(dirname \$(python -c "import nvidia.cudnn; print(nvidia.cudnn.__file__)"))
   export LD_LIBRARY_PATH=\${CUDNN_PATH}/lib
   >&2 echo "CuDNN path at" \$CUDNN_PATH "contains:"  # should exist and give a good path
   >&2 echo \$(ls \$CUDNN_PATH)                       # should contain stuff like a lib subdir with libcudnn .so files
   >&2 echo "LD library path at" \$LD_LIBRARY_PATH    # should exist and contain CUDNN_PATH
   >&2 echo "\$(nvcc --version)"
   >&2 python3 -c "import tensorflow as tf; print(f'Available devices:\\n{tf.config.list_physical_devices()}')"
   XLA_PYTHON_CLIENT_MEM_FRACTION=.9 sppid af2-single \
      <(echo "ref/ref.a3m") \
      --msa2 \$MSA_LIST \
      --list-file \
      --output "${organism_id}-${uniprot_id}-batch_${batch_idx}_af2" \
      --plot "${organism_id}-${uniprot_id}-batch_${batch_idx}_af2"
   """

   stub:
   """
   mkdir "${organism_id}-${uniprot_id}-batch_${batch_idx}_af2"
   touch "${organism_id}-${uniprot_id}-batch_${batch_idx}_af2/stub.pdb"
   echo "Skipping AF2 for stub"
   """
}


/*
========================================================================================
   Workflow Event Handler
========================================================================================
*/

workflow.onComplete {

   println ( workflow.success ? """
      Pipeline execution summary
      ---------------------------
      Completed at: ${workflow.complete}
      Duration    : ${workflow.duration}
      Success     : ${workflow.success}
      workDir     : ${workflow.workDir}
      exit status : ${workflow.exitStatus}
      """ : """
      Failed: ${workflow.errorReport}
      exit status : ${workflow.exitStatus}
      """
   )
}

/*
========================================================================================
   THE END
========================================================================================
*/