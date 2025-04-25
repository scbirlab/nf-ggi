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

pipeline_title = """\
   S C B I R   G E N E - G E N E   I N T E R A C T I O N   P R E D I C T I O N   P I P E L I N E
   =============================================================================================
   Nextflow pipeline to predict gene-gene interactions based on protein-protein interaction 
   predictions and similar metabolites.
   """
   .stripIndent()

/*
========================================================================================
   Help text
========================================================================================
*/
if ( params.help ) {
   println pipeline_title + """\
         Command-line usage:
            nextflow run scbirlab/nf-ggi --uniclust <path> --bfd <path> --organism_id <taxon ID>
            nextflow run scbirlab/nf-ggi --uniclust <path> --bfd <path> --organism_id <taxon ID> --bait <UniProtID>
            nextflow run scbirlab/nf-ggi --uniclust <path> --bfd <path> --organism_id <taxon ID> --bait <taxon ID> --bait_is_taxon --interspecies
            nextflow run scbirlab/nf-ggi --uniclust <path> --bfd <path> --organism_id <taxon ID> --filename <path> --column1 <gene-col1> --column2 <gene-col2> [--interspecies --organism_id2 <taxon ID>] [--format <gene-name-type>]
         Config/sample sheet usage:
            nextflow run scbirlab/nf-ggi -c <config-file>

         Command-line required parameters:
            --organism_id             Taxon ID for organism
            Bait mode:
               --bait                 UniProt ID for bait protein, or Taxon ID for bait organism
            Custom mode:
               --filename             Filename to get custom protein pairs
               --column1, --column2   Column names from --filename to get protein IDs

         Command-line optional parameters:
            --bait_is_taxon  Indicate that bait is an organism ID
            --interspecies   Run analysis between interacting species proteomes
            --organism_id2   When providing a file of pairs, if the second protein (--column2) is from another organism than the first
            --format         Type of gene identifier in --column1, --column2. Default: "Gene_Name"
            --test           Whether to run in test mode. Default: false.
            --outputs        Output folder. Default: "outputs".
            --batch_size     What size to batch protein-protein interactions into. Default: 100.
            --plots          Generate contact map plots

         Config required parameters:
            uniclust, bfd        Paths to get HHblits databases.
            *and* either:
               organism_id       Taxon ID for organism
               (And all the same named flags above for command-line)
            *or*
               sample_sheet      CSV file with columns with same names as command-line flags, one row per combination to run
               mode              "self" (all vs all), "bait" (all vs some), "custom" (some-vs-some)

         Config optional parameters (with defaults):  
            test           Whether to run in test mode. Default: false.
            batch_size     What size to batch protein-protein interactions into. Default: 100.
            rhea_url       URL to download Rhea reaction database. Default: "https://ftp.expasy.org/databases/rhea"
            outputs        Output folder. Default: "outputs".

         The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.
   
   """
   .stripIndent()
   exit 0
}

/*
========================================================================================
   Check parameters
========================================================================================
*/
if ( !params.sample_sheet ) {
   if ( !params.organism_id ) {
      throw new Exception("!!! PARAMETER MISSING: Please provide a sample sheet or at least --organism_id.")
   }
   if ( params.filename ) {
      if ( !params.format ) {
         throw new Exception("!!! PARAMETER MISSING: Please provide a --format when using --filename.")

      }
      if ( !params.column1 ) {
         throw new Exception("!!! PARAMETER MISSING: Please provide a --column1 when using --filename.")
         
      }
      if ( !params.column2 ) {
         throw new Exception("!!! PARAMETER MISSING: Please provide a --column2 when using --filename.")
         
      }
      
   }

}
if ( !params.uniclust ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to UniClust database.")
}
if ( !params.bfd ) {
   throw new Exception("!!! PARAMETER MISSING: Please provide a path to BFD database.")
}

log.info pipeline_title + """\
   test mode               : ${params.test}
   mode                    : ${params.mode}
      Bait is Taxon ID     : ${params.bait_is_taxon}
      Interspecies         : ${params.interspecies}
   inputs
      input_dir            : ${params.inputs}
      sample sheet         : ${params.sample_sheet}
      UniClust database    : ${params.uniclust}
      BFD                  : ${params.bfd}
      Rhea URL             : ${params.rhea_url}
   batch size              : ${params.batch_size}
   co-evolution analysis
      Metabolites          : ${params.metabolites}
      STRINGdb             : ${params.string}
      DCA                  : ${params.dca}
      RosettaFold-2track   : ${params.rf2t}           
      AlphaFold2           : ${params.af2}              
   output                  : ${params.outputs}
      make plots?          : ${params.plots}
   """
   .stripIndent()


/*
========================================================================================
   MAIN Workflow
========================================================================================
*/

// load modules
include {
   make_msa_from_fasta;
   make_msa_from_fasta as make_msa_from_fasta_bait;
} from './modules/msa.nf'
include {
   fetch_rhea_database;
   match_uniprot_to_reactants;
} from './modules/rhea.nf'
include {
   fetch_string_database;
} from './modules/string.nf'
include { 
   fetch_fastas_from_organism_id;
   fetch_fastas_from_organism_id as fetch_fastas_from_organism_id_bait;
   fetch_fasta_from_uniprot_id;
   map_uniprot_ids_from_file;
   map_uniprot_ids_from_file as map_uniprot_ids_from_file_bait;
   map_gene_names_from_file as map_gene_names_from_file1;
   map_gene_names_from_file as map_gene_names_from_file2;
} from './modules/uniprot.nf'
include {
   run_dca;
   run_rf2track;
   run_af2;
   stack_table as stack_dca;
   stack_table as stack_rf2t;
   stack_table as stack_af2;
} from './modules/yunta.nf'
// include { GENE_NAME_TO_UNIPROT as GENE_NAME_TO_UNIPROT } from './modules/uniprot.nf'

workflow {

   Channel.value( params.bfd ).set { bfd }
   Channel.value( params.uniclust ).set { uniclust }
   Channel.of( params.rhea_url ).set { rhea_url }

   if ( params.sample_sheet ) {

      mode = params.mode

      Channel.fromPath( 
         params.sample_sheet,
         checkIfExists: true, 
      )
         .splitCsv( header: true )
         .set { sample_rows }

   }

   else {

      Channel.of( [
         organism_id: params.organism_id,
         organism_id2: params.organism_id2,
         bait: params.bait,
         filename: params.filename,
         format: params.format,
         column1: params.column1,
         column2: params.column2,
      ] )
      .set { sample_rows }

      if ( params.bait ) {
         mode = "bait"
      }
      else if ( params.filename ) {
         mode = "custom"
      }
      else {
         mode = "self"
      }

   }
   

   if ( mode == 'self' || mode == 'bait' ) {

      sample_rows
         .map { tuple( it.organism_id, it.organism_id ) }
         .unique()
         | fetch_fastas_from_organism_id  // Organism ID, FASTAs gz
      fetch_fastas_from_organism_id.out
         .splitFasta( elem: 1, record: [id: true, text: true] )  // Organism ID, FASTA text
         .map { [ it[0] ] + it[1].id.split('\\|')[1..2] + [ it[1].text ] }  // Organism ID, UniProtID, Entry Name, FASTA text
         .set { fastas_A0 }

      ( params.test ? fastas_A0.take(3) : fastas_A0 )
         .set { fastas_A }

      if ( mode == 'self' ) {

         fastas_A.set { fastas_B }

      }

      else {

         sample_rows
            .map { tuple( it.organism_id, it.bait ) }
            .unique()
            .set { baits }

         if ( params.bait_is_taxon ) {

            baits
               | fetch_fastas_from_organism_id_bait  // Organism ID, bait FASTAs gz
            fetch_fastas_from_organism_id_bait.out
               .set { bait_fastas }

         }

         else {

            baits
               | fetch_fasta_from_uniprot_id  // Organism ID, bait FASTA file
            fetch_fasta_from_uniprot_id.out
               .set { bait_fastas }
   
         }

         bait_fastas  // Organism ID, bait FASTA file
            .splitFasta( elem: 1, record: [id: true, text: true] )   // Organism ID, bait FASTA text
            .map { 
               [ it[0] ] + it[1].id.split('\\|')[1..2] + [ it[1].text ] 
            }  // Organism ID, UniProtID, Entry Name, FASTA text
            .set { fastas_B }
         
      }

   }

   else if ( mode == 'custom' ) {

        sample_rows
            .map { tuple(
               tuple( 
                  it.organism_id,
                  ( params.interspecies ? it.organism_id2 : it.organism_id ),
               ), 
               file( 
                  "${params.sample_sheet ? params.inputs : '.'}/${it.filename}", 
                  checkIfExists: true,
               ),
               it.format,
               it.column1,
               it.column2,
            ) }
            .set { mapping_input }

         mapping_input
            .map { tuple( it[0][0], it[0][0], it[1], it[2], it[3] ) }
            | map_uniprot_ids_from_file  // Organism ID, UniProtID A file
         
         mapping_input
            .map { tuple( it[0][0], it[0][1], it[1], it[2], it[4] ) }
            | map_uniprot_ids_from_file_bait  // Organism ID, UniProtID B file

         map_uniprot_ids_from_file.out
            .splitText( elem: 1 ) { v -> v.collect { it.toString().trim() } }
            .set { uniprot_A0 }  // Organism ID, UniProtID A 
         map_uniprot_ids_from_file_bait.out
            .splitText( elem: 1 ) { v -> v.collect { it.toString().trim() } }
            .set { uniprot_B0 }  // Organism ID, UniProtID B 

         if ( params.test ) {
            uniprot_A0.take(3).set { uniprot_A }
            uniprot_B0.take(3).set { uniprot_B }
         }

         else {
            uniprot_A0.set { uniprot_A }
            uniprot_B0.set { uniprot_B }
         }

         uniprot_A
            .map { tuple( it[0], it[1] ) }
            .concat(
               uniprot_B.map { tuple( it[0], it[1] ) }
            )
            .unique()
            | fetch_fasta_from_uniprot_id   // Organism ID, FASTA file

         fetch_fasta_from_uniprot_id.out
            .splitFasta( elem: 1, record: [id: true, text: true] )   // Organism ID, bait FASTA text
            .map { [ it[0] ] + it[1].id.split('\\|')[1..2] + [ it[1].text ] }  // Organism ID, UniProtID, Entry Name, FASTA text
            .set { all_fasta_custom }
         
         all_fasta_custom
            .combine( uniprot_A, by: [0, 1] )
            .set { fastas_A }
         all_fasta_custom
            .combine( uniprot_B, by: [0, 1] )
            .set { fastas_B }

   }

   else {
      error "Unsupported mode: ${mode}. Choose from: organism, bait, custom."
   }

   if ( params.test ) {
      fastas_A.take(2).set { fastas_A2 }
      fastas_B.take(2).set { fastas_B2 }
   }

   else {
      fastas_A.set { fastas_A2 }
      fastas_B.set { fastas_B2 }
   }

   fastas_A2
      .concat( fastas_B2 ) // Organism ID, UniProtID, Entry Name, FASTA text
      .unique()
      .tap { all_fasta_info }
      .map { it[1..0] }  // UniProtID, Organism ID 
      .unique()
      .set { id_to_uniprot_map }

   if ( params.metabolites ) {
      // Get reactants from Rhea database using UniProtIDs
      fetch_rhea_database( rhea_url )
      all_fasta_info  // Organism ID, UniProtID, Entry Name, FASTA text
         .collectFile( newLine: true ) { [ "${it[0]}.txt", it[1] ] }
         .map { tuple( it.getSimpleName(), it ) }  // Organism ID, UniProtID list file
         .combine( fetch_rhea_database.out ) // Organism ID, UniProtID list file, Rhea SMILES, Rhea2UniProt
         | match_uniprot_to_reactants  // Organism ID, reaction SMILES file
         | CONNECT_METABOLITES 
   }

   if ( params.string ) {
      // Get STRING connections
      all_fasta_info  // Organism ID, UniProtID, Entry Name, FASTA text
         .map { it[0..1] }  // Organism ID, UniProtID
         .collectFile( newLine: true ) { [ "${it[0]}.txt", it[1] ] }  // UniProtID list
         .map { tuple( it.getSimpleName(), it ) }  // Organism ID, UniProtID list
         | fetch_string_database
   }

   all_fasta_info
      .map { tuple( it[1], it[3] ) }  // UniProtID, FASTA text
      .collectFile( newLine: true ) { [ "${it[0]}.fasta", "${it[1]}" ] }
      .map { tuple( it.getSimpleName(), it ) }
      .unique()
      .set { input_for_making_msas }
   make_msa_from_fasta(
      input_for_making_msas,
      uniclust,
      bfd,
   )
   make_msa_from_fasta.out  // UniProtID, MSA
      .combine(
         id_to_uniprot_map,
         by: 0,
      )  // UniProtID, MSA, Organism ID 
      .map { tuple( it[2], it[0], it[1] ) }  // Organism ID, UniProtID, MSA
      .tap { all_msas }
      .join(
         fastas_A.map{ it[0..1] },
         by: [0, 1],
         failOnDuplicate: true,
      )  // Organism ID, UniProtID, MSA
      .set { msa_A }

   all_msas
      .join(
         fastas_B2.map{ it[0..1] },
         by: [0, 1],
         failOnDuplicate: true,
      )  // Organism ID, bait UniProtID, bait MSA
      .set { msa_B }

   msa_B  // Organism ID, Bait UniProt ID, Bait MSA
      .combine(
         msa_A,
         by: 0
      )  // Organism ID, Bait UniProt ID, Bait MSA, UniProtID, MSA
      .set { crossed_msa } 

   ( 
      mode == 'self' 
      ? crossed_msa.filter { it[1] < it[-2] } 
      : crossed_msa.filter { it[1] != it[-2] } 
   )  // If self-cross, only take lower triangle
      .groupTuple( 
         by: [0, 1, 2]
      )  // Organism ID, Bait UniProt ID, Bait MSA, [UniProtID, ...], [MSA, ...]
      .map { 
         tuple(
            it[0], it[1], it[2],
            it[3].withIndex().collect { 
               el, i -> Math.round(Math.floor(i / params.batch_size)) 
            },
            it[3], it[4]
         ) 
      }  // Organism ID, Bait UniProt ID, Bait MSA, [batch_i, ...], [UniProtID, ...], [MSA, ...]
      .transpose()  // Organism ID, Bait UniProt ID, Bait MSA, batch_i, UniProtID, MSA
      .groupTuple( 
         by: [0, 1, 2, 3]
      )  // Organism ID, Bait UniProt ID, Bait MSA, batch_i, [UniProtID, ...], [MSA, ...]
      .filter { it[-1].size() > 0 }  // filter out trivial (size-0) elements
      .set { msa_pairs0 }

   if ( params.test ) {
      msa_pairs0.take(2).set { msa_pairs }
   }

   else {
      msa_pairs0.set { msa_pairs }
   }

   // Calculate evolutionary coupling
   Channel.value( params.bait_is_taxon || params.interspecies ).set { interspecies }
   Channel.value( params.plots ).set { make_coevo_plots }
   if ( params.dca ) {
      run_dca( 
         msa_pairs, 
         interspecies, 
         make_coevo_plots,
      )
      stack_dca(
         run_dca.out.main
            .groupTuple( by: 0 ),  // Organism ID, [tsv, ...]
         Channel.value( "dca" )
      ) // Organism ID, tsv
         | set { stacked_dca }
   }

   else {

      Channel.empty() 
         .set { stacked_dca }

   }

   if ( params.rf2t ) {
      run_rf2track( 
         msa_pairs, 
         interspecies, 
         make_coevo_plots,
      )
      stack_rf2t(
         run_rf2track.out.main
            .groupTuple( by: 0 ),  // Organism ID, [tsv, ...]
         Channel.value( "rf2t" )
      )  // Organism ID, tsv
         | set { stacked_rf2t }
   }

   else {

      Channel.empty() 
         .set { stacked_rf2t }

   }

   if ( params.af2 ) {
      run_af2( 
         msa_pairs, 
         interspecies, 
         make_coevo_plots,
      )
      stack_af2(
         run_af2.out.main
            .groupTuple( by: 0 ),  // Organism ID, [tsv, ...]
         Channel.value( "af2" )
      )  // Organism ID, tsv
         | set { stacked_af2 }
   }

   else {

      Channel.empty() 
         .set { stacked_af2 }

   }

   stacked_dca
      .concat(
         stacked_rf2t,
         stacked_af2,
      )
      .set { ppi_outputs }
   map_gene_names_from_file1(
      ppi_outputs,
      Channel.value( "uniprot_id_1" )
   )
   map_gene_names_from_file2(
      map_gene_names_from_file1.out,
      Channel.value( "uniprot_id_2" )
   )

}

process CONNECT_METABOLITES {

   tag "${id}"
   label "big_time"

   publishDir( 
      "${params.outputs}/metabolites", 
      mode: 'copy',
      saveAs: { "${id}.metabolism-connection.tsv" }
   )

   input:
   tuple val( id ), path( reaction_table )

   output:
   tuple val( id ), path( "metabolism-connection.tsv" )

   script:
   """
   python ${projectDir}/bin/metabolism/connect-metabolites.py \
   < "${reaction_table}" \
   > metabolism-connection.tsv
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