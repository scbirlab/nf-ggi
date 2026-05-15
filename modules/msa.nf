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
      parallel -j ${task.cpus} \
         "hhblits \
            -cpu ${task.cpus} \
            -maxmem ${task.memory.getGiga()} \
            -v 2 \
            -i {} \
            -d \$d \
            -e 0.001 \
            -o /dev/null \
            -oa3m {.}.""\$(basename \$d)"".a3m \
            -o /dev/null \
            -cov 60 \
            -n 3 \
            -realign -realign_max 10000" \
      ::: *.fasta
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


process Colabfold_MSA {

   tag "${id}"
   label 'big_cpu_mem'
   // container 'ghcr.io/soedinglab/mmseqs2:latest'
   // label 'gpu_single_short'
   container 'ghcr.io/soedinglab/mmseqs2:master-cuda12'

   // errorStrategy 'retry'  // sometimes cluster will kill the job
   // maxRetries 1

   publishDir(
        "${params.outputs}/msa",
        mode: 'copy',
        saveAs: { "${fasta.simpleName}.a3m" }
   )

   input:
   tuple val( id ), file( fasta )
   tuple val( uniref_root ), path( uniref )
   tuple val( env_root ), path( bfd )
   val use_gpu

   output:
   tuple val( id ), file( '*.a3m' )

   script:
   // TODO: Allow GPU usage. Might need to have special databases.
   def gpu_flags = use_gpu ? "--gpu 1 --prefilter-mode 1" : "--prefilter-mode 1 --k-score 'seq:96,prof:80'"
   """
   # TODO: Allow GPU usage. Make GPU-compatible database in a separate process?
   #${use_gpu ? "mmseqs makepaddedseqdb targetDB targetDB_gpu && mmseqs rmdb targetDB && mv targetDB_gpu targetDB" : ""}
   set -euox pipefail

   BASE_FLAGS="--db-load-mode 2 --threads ${task.cpus}"

   mmseqs createdb "${fasta}" query

   dbs=("${uniref_root}" "${env_root}")
   for db in \${dbs[@]}
   do

      if [ "\$db" == "${uniref_root}" ]
      then
         QUERY=query
      else
         QUERY=prof_result_uniref
      fi

      mmseqs search "\$QUERY" "\$db" result_"\$db" tmp_"\$db" \
         --threads ${task.cpus} \
         --num-iterations 3 \
         --db-load-mode 2 \
         --prefilter-mode 0 \
         -a \
         -e 0.1 \
         -s 8 \
         --max-seqs 10000

      if [ "\$db" == "${uniref_root}" ]
      then
         # Extract profile from last iteration
         mmseqs mvdb tmp_"\$db"/latest/profile_1 prof_result_uniref
         mmseqs lndb query_h prof_result_uniref_h
         EXPAND_QUERY="query"
         EXPAND_FLAGS="--expand-filter-clusters 1 --max-seq-id 0.95"
         ALIGN_QUERY="prof_result_uniref"
      else
         EXPAND_QUERY=tmp_"\$db"/latest/profile_1
         EXPAND_FLAGS=""
         ALIGN_QUERY=tmp_"\$db"/latest/profile_1"
      fi

      # Expand: fetch all cluster members for matched representatives
      mmseqs expandaln "\$EXPAND_QUERY" "\$db.idx" result_"\$db" "\$db.idx" result_exp_"\$db" \
         \$BASE_FLAGS \
         --expansion-mode 0 \
         -e inf \$EXPAND_FLAGS

      # Realign expanded hits against the profile
      mmseqs align "\$ALIGN_QUERY" "\$db.idx" result_exp_"\$db" result_exp_realign_"\$db" \
         \$BASE_FLAGS \
         -e 10 \
         --max-accept 100000 \
         --alt-ali 10 \
         -a

      # Filter
      mmseqs filterresult query "\$db.idx" result_exp_realign_"\$db" result_exp_realign_filter_"\$db" \
         \$BASE_FLAGS \
         --qid 0 \
         --qsc 0.8 \
         --diff 0 \
         --max-seq-id 1.0 \
         --filter-min-enable 100

      # Write A3M with diversity subsampling
      mmseqs result2msa query "\$db.idx" result_exp_realign_filter_"\$db" "\$db".a3m \
         \$BASE_FLAGS \
         --msa-format-mode 6 \
         --filter-msa 1 \
         --filter-min-enable 1000 \
         --diff 3000 \
         --qid 0.0,0.2,0.4,0.6,0.8,1.0 \
         --qsc 0 \
         --max-seq-id 0.95

      # Clean up
      mmseqs rmdb result_"\$db"
      mmseqs rmdb result_exp_"\$db"
      mmseqs rmdb result_exp_realign_"\$db"
      mmseqs rmdb result_exp_realign_filter_"\$db"
   done

   # Merge and clean up
   mmseqs mergedbs query msa.a3m *.a3m
   for f in *.a3m
   do
      [[ "\$f" == "msa.a3m" ]] && continue
      mmseqs rmdb "\$f"
   done

   mmseqs unpackdb msa.a3m . --unpack-name-mode 1 --unpack-suffix .a3m
   mmseqs rmdb msa.a3m

   mmseqs rmdb prof_result_uniref
   mmseqs rmdb prof_result_uniref_h
   rm -rf tmp_*

   """
   stub:
   """
   head -n2 "${fasta}" > "msa.a3m"
   """
}
