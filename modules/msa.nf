process make_msa_from_fasta {

   tag "${id}"
   label 'big_cpu_mem'
   
   errorStrategy 'retry'
   maxRetries 2

   publishDir( 
      "${params.outputs}/msa", 
      mode: 'copy',
      saveAs: { "${fasta.getSimpleName()}.a3m" }
   )

   // Proteome ID, UniProtID, FASTA file, uniclust, bfd
   input:
   tuple val( id ), file( fasta )
   val uniclust 
   val bfd 

   output:
   tuple val( id ), file ( 'msa.a3m' )
   // tuple val( fasta.getSimpleName() ), path( "*.a3m" )

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
   # concatenate from both databases
   cat <(head -n2 \${outputs[0]}) <(tail -n+3 -q \${outputs[@]}) \
      > "msa.a3m"
   for f in \${outputs[@]}
   do
      if [ \$f != "msa.a3m" ]
      then
         rm \$f
      fi
   done
   """

   stub:
   """
   head -n2 ${fasta} > "msa.a3m"
   """
}