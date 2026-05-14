process make_msa_from_fasta {

   tag "${id}"
   label 'big_cpu_mem'
   
   errorStrategy 'retry'  // sometimes cluster will kill the job
   maxRetries 1

   publishDir( 
      "${params.outputs}/msa", 
      mode: 'copy',
      saveAs: { "${fasta.getSimpleName()}.a3m" }
   )

   // Proteome ID, UniProtID, FASTA file, uniclust, bfd
   input:
   tuple val( id ), file( fasta )
   tuple val( uniclust_root ), path( uniclust ) 
   tuple val( bfd_root ), path( bfd )

   output:
   tuple val( id ), file ( 'msa.a3m' )
   // tuple val( fasta.getSimpleName() ), path( "*.a3m" )

   script:
   """
   set -x
   dbs=(${uniclust_root} ${bfd_root})
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


process Make_msa_from_fasta_with_MMSeqs2 {

   tag "${id}"
   label 'big_cpu_mem'
   // container 'ghcr.io/soedinglab/mmseqs2:latest'
   // label 'gpu_single_short'
   container 'ghcr.io/soedinglab/mmseqs2:master-cuda12'

   //  errorStrategy 'retry'  // sometimes cluster will kill the job
   // maxRetries 1

   publishDir(
        "${params.outputs}/msa",
        mode: 'copy',
        saveAs: { "${fasta.getSimpleName()}.a3m" }
   )

   input:
   tuple val( id ), file( fasta )
   tuple val( uniref_root ), path( uniref )
   tuple val( bfd_root ), path( bfd )
   val use_gpu
   output:
   tuple val( id ), file( 'msa.a3m' )

   script:
   """
   set -euox pipefail

   mmseqs createdb "${fasta}" queryDB
   #${use_gpu ? "mmseqs makepaddedseqdb targetDB targetDB_gpu && mmseqs rmdb targetDB && mv targetDB_gpu targetDB" : ""}

   dbs=(${uniref_root} ${bfd_root})

   for d in \${dbs[@]}
   do
       name="\$(basename \$d)"
       mkdir -p "tmp_\$name"
       # TODO: Allow GPU usage
       mmseqs search queryDB "\$d" "result_\$name" "tmp_\$name" \
           --threads ${task.cpus} \
           -e 0.001 \
           --num-iterations 3 \
           -s 8 \
           --max-seqs 10000 \
           --db-load-mode 2
       # profile expansion for MSA depth ~ HHblits iterative behaviour
       mmseqs expandaln queryDB "\$d" "result_\$name" "expanded_\$name" \
           --expansion-mode 0 \
           -e 1e-6 \
           --expand-filter-clusters 1 \
           --max-seq-id 0.95 \
           --threads ${task.cpus}
       mmseqs result2msa queryDB "\$d" "expanded_\$name" "\${name}.a3m" \
           --msa-format-mode 6 \
           --cov 0.6 --cov-mode 1 \
           --threads ${task.cpus}
       rm -rf "tmp_\$name"
   done

   outputs=( *.a3m )
   head -n2 "\${outputs[0]}" | cat - <(tail -n+3 -q \${outputs[@]}) > "msa.a3m"
   for f in \${outputs[@]}
   do
       [ \$f != "msa.a3m" ] && rm \$f
   done
   """
   stub:
   """
   head -n2 "${fasta}" > "msa.a3m"
   """
}
