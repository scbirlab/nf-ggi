process HHblits_MSA {

   tag "${fasta[0]}...${fasta[-1]}"
   label 'big_cpu_mem'
   
   errorStrategy 'retry'  // sometimes cluster will kill the job
   maxRetries 1

   publishDir( 
      "${params.outputs}/msa/hhblits", 
      mode: 'copy',
      // saveAs: { "${fasta.simpleName}.a3m" }
   )

   // Proteome ID, UniProtID, FASTA file, uniclust, bfd
   input:
   path( fasta )
   tuple val( uniclust_root ), path( uniclust ) 
   tuple val( bfd_root ), path( bfd )

   output:
   path( "*.a3m" )
   // tuple val( fasta.getSimpleName() ), path( "*.a3m" )

   script:
   def mem_per_job = (task.memory.getGiga() / task.cpus).toInteger()
   """
   set -euox pipefail

   dbs=("${uniclust_root}" "${bfd_root}")
   for db in \${dbs[@]}
   do
      parallel -j ${task.cpus} \
         "hhblits \
            -cpu 1 \
            -maxmem ${mem_per_job} \
            -v 2 \
            -i {} \
            -d \$db \
            -e 0.001 \
            -oa3m {.}.""\$(basename \$db)"".out.a3m \
            -o /dev/null \
            -cov 60 \
            -n 3 \
            -realign -realign_max 10000" \
      ::: *.fasta
   done

   # Per-sequence merge: uniclust header + homologs from both DBs
   for f in *.fasta
   do
      base="\${f%.fasta}"
      outputs=( "\${base}".*.out.a3m )
      # concatenate from both databases
      cat <(head -n2 \${outputs[0]}) <(tail -n+3 -q \${outputs[@]}) \
         > "\${base}.a3m"
      rm \${outputs[@]}
   done
   """

   stub:
   """
   for f in *.fasta
   do
      base="\${f%.fasta}"
      head -n2 "\$f" > "\${base}.a3m"
   done
   """
}


process Colabfold_MSA {

   tag "${fasta[0]}...${fasta[-1]}"
   label 'big_cpu_mem'
   // container 'ghcr.io/soedinglab/mmseqs2:latest'
   // label 'gpu_single_short'
   container 'ghcr.io/soedinglab/mmseqs2:master-cuda12'

   // errorStrategy 'retry'  // sometimes cluster will kill the job
   // maxRetries 1

   publishDir(
        "${params.outputs}/msa/colabfold",
        mode: 'copy',
      //   saveAs: { "${fasta.simpleName}.a3m" }
   )

   input:
   path( fasta )
   tuple val( uniref_root ), path( uniref )
   tuple val( env_root ), path( bfd )
   val use_gpu

   output:
   path( '*.a3m' )

   script:
   // TODO: Allow GPU usage. Might need to have special databases.
   def gpu_flags = use_gpu ? "--gpu 1 --prefilter-mode 1" : "--prefilter-mode 1 --k-score 'seq:96,prof:80'"
   """
   # TODO: Allow GPU usage. Make GPU-compatible database in a separate process?
   #${use_gpu ? "mmseqs makepaddedseqdb targetDB targetDB_gpu && mmseqs rmdb targetDB && mv targetDB_gpu targetDB" : ""}
   set -euox pipefail

   BASE_FLAGS="--db-load-mode 2 --threads ${task.cpus}"

   mmseqs createdb *.fasta query

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
         \$BASE_FLAGS \
         --num-iterations 3 \
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
         EXPAND_QUERY=prof_result_uniref
         EXPAND_FLAGS=""
         ALIGN_QUERY=tmp_"\$db"/latest/profile_1
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

   mmseqs unpackdb msa.a3m . --unpack-name-mode 1 --unpack-suffix .temp.a3m
   mmseqs rmdb msa.a3m
   # strip db prefix and description, keep bare accession
   for f in *.temp.a3m
   do
      header=\$(head -n1 "\$f")
      acc=\$(
         echo "\$header" \
         | awk '
            {
               id = substr(\$1, 2)
               n = split(id, a, "|")
               print (n >= 3) ? a[2] : a[1]
            }
         '
      )
      [[ -n "\$acc" ]] && mv "\$f" "\${acc}.a3m"
   done

   mmseqs rmdb prof_result_uniref
   mmseqs rmdb prof_result_uniref_h
   rm -rf tmp_*

   """
   stub:
   """
   for f in *.fasta
   do
      base="\${f%.fasta}"
      head -n2 "\$f" > "\${base}.a3m"
   done
   """
}
