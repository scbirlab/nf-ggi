process HHblits_MSA {

   tag "${fasta[0]}...${fasta[-1]}"
   label 'big_mem'
   
   errorStrategy 'retry'  // sometimes cluster will kill the job
   maxRetries 2

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
   def mem_per_job = Math.floor(0.8 * task.memory.getGiga() / task.cpus).toInteger()
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
            -realign -realign_max 10000 \
            || [[ -s {.}.\$(basename \$db).out.a3m ]]" \
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
   label 'all_cpu_mem'

   container 'ghcr.io/soedinglab/mmseqs2:15-6f452'
   // label 'gpu_single_short'
   // container 'ghcr.io/soedinglab/mmseqs2:master-cuda12'

   errorStrategy 'retry'  // sometimes cluster will kill the job, or mmseqs2 segfaults
   maxRetries 2

   publishDir(
        "${params.outputs}/msa/colabfold",
        mode: 'copy',
      //   saveAs: { "${fasta.simpleName}.a3m" }
   )

   input:
   path fasta
   tuple val( uniref_root ), path( uniref )
   tuple val( env_root ), path( bfd )
   val use_gpu

   output:
   path '*.a3m.fasta'

   script:
   def split_mem = Math.floor(task.memory.getGiga() * 0.9).toInteger()
   // TODO: Allow GPU usage. Might need to have special databases.
   def gpu_flags = use_gpu ? "--gpu 1 --prefilter-mode 1" : "--prefilter-mode 1 --k-score 'seq:96,prof:80'"
   def cmd = use_gpu ? "entrypoint" : "entrypoint"
   """
   # TODO: Allow GPU usage. Make GPU-compatible database in a separate process?
   #${use_gpu ? "mmseqs makepaddedseqdb targetDB targetDB_gpu && mmseqs rmdb targetDB && mv targetDB_gpu targetDB" : ""}
   set -euox pipefail

   BASE_FLAGS="--db-load-mode 0 --threads ${task.cpus}"
   MEM_FLAG=" --split-memory-limit ${split_mem}G --split 0"

   "${cmd}" createdb *.fasta query

   dbs=("${uniref_root}" "${env_root}")
   for db in \${dbs[@]}
   do
      DB_SEQ="\$db".idx #_seq
      DB_ALN="\$db.idx" #_aln
      if [ "\$db" == "${uniref_root}" ]
      then
         QUERY=query
      else
         QUERY=prof_result_uniref
      fi

      "${cmd}" search "\$QUERY" "\$db" result_"\$db" tmp_"\$db" \
         \$BASE_FLAGS \$MEM_FLAG \
         --num-iterations 3 \
         --prefilter-mode 0 \
         -a \
         -e 0.1 \
         -s 8 \
         --max-seqs 10000

      if [ "\$db" == "${uniref_root}" ]
      then
         # Extract profile from last iteration
         "${cmd}" mvdb tmp_"\$db"/latest/profile_1 prof_result_uniref
         "${cmd}" lndb query_h prof_result_uniref_h
         EXPAND_QUERY="query"
         EXPAND_FLAGS="--expand-filter-clusters 1 --max-seq-id 0.95"
         ALIGN_QUERY="prof_result_uniref"
      else
         EXPAND_QUERY=prof_result_uniref
         EXPAND_FLAGS=""
         ALIGN_QUERY=tmp_"\$db"/latest/profile_1
      fi

      # Expand: fetch all cluster members for matched representatives
      "${cmd}" expandaln "\$EXPAND_QUERY" "\$DB_SEQ" result_"\$db" "\$DB_ALN" result_exp_"\$db" \
         \$BASE_FLAGS \
         --expansion-mode 0 \
         -e inf \$EXPAND_FLAGS

      # Realign expanded hits against the profile
      "${cmd}" align "\$ALIGN_QUERY" "\$DB_SEQ" result_exp_"\$db" result_exp_realign_"\$db" \
         \$BASE_FLAGS \
         -e 10 \
         --max-accept 100000 \
         --alt-ali 10 \
         -a

      # Filter
      "${cmd}" filterresult query "\$DB_SEQ" result_exp_realign_"\$db" result_exp_realign_filter_"\$db" \
         \$BASE_FLAGS \
         --qid 0 \
         --qsc 0.8 \
         --diff 0 \
         --max-seq-id 1.0 \
         --filter-min-enable 100

      # Write A3M with diversity subsampling
      "${cmd}" result2msa query "\$DB_SEQ" result_exp_realign_filter_"\$db" "\$db".aln.fasta \
         \$BASE_FLAGS \
         --msa-format-mode 2 \
         --filter-msa 1 \
         --filter-min-enable 1000 \
         --diff 3000 \
         --qid 0.0,0.2,0.4,0.6,0.8,1.0 \
         --qsc 0 \
         --max-seq-id 0.95

      # Clean up
      "${cmd}" rmdb result_"\$db"
      "${cmd}" rmdb result_exp_"\$db"
      "${cmd}" rmdb result_exp_realign_"\$db"
      "${cmd}" rmdb result_exp_realign_filter_"\$db"
   done

   # Merge and clean up
   "${cmd}" mergedbs query msa.aln.fasta *.aln.fasta
   for f in *.aln.fasta
   do
      [[ "\$f" == "msa.aln.fasta" ]] && continue
      "${cmd}" rmdb "\$f"
   done

   "${cmd}" unpackdb msa.aln.fasta . --unpack-name-mode 1 --unpack-suffix .temp.fasta
   "${cmd}" rmdb msa.aln.fasta
   # strip db prefix and description, keep bare accession
   for f in *.temp.fasta
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
      [[ -n "\$acc" ]] && mv "\$f" "\${acc}.a3m.fasta"
   done

   "${cmd}" rmdb prof_result_uniref
   "${cmd}" rmdb prof_result_uniref_h
   rm -rf tmp_*

   """
   stub:
   """
   for f in *.fasta
   do
      base="\${f%.fasta}"
      head -n2 "\$f" > "\${base}.a3m.fasta"
   done
   """
}

process Convert_FASTA_to_A3M {

   tag "${fasta[0]}...${fasta[-1]}"
   // label 'all_cpu_mem'

   // errorStrategy 'retry'  // sometimes cluster will kill the job
   // maxRetries 1

   publishDir(
      "${params.outputs}/msa/colabfold-a3m",
      mode: 'copy',
   //   saveAs: { "${fasta.simpleName}.a3m" }
   )

   input:
   path fasta

   output:
   path '*.a3m'

   script:
   """
   #!/usr/bin/env python
   from glob import glob
   import os
   import sys

   for filename in glob("*.a3m.fasta"):
      seqs, order, cur = {}, [], None
      with open(filename, "r") as f:
         for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                  cur = line[1:]
                  order.append(cur)
                  seqs[cur] = ''
            else:
                  seqs[cur] += line

      query = seqs[order[0]]
      match_cols = {i for i, c in enumerate(query) if c != '-'}
      outfile = filename.replace('.a3m.fasta', '.a3m')
      with open(outfile, 'w') as out:
         for name in order:
            print('>' + name, file=out)
            s = seqs[name]
            if name == order[0]:
               # Query: defines match columns, no gaps/insertions
               print(query.replace('-', ''), file=out)
            else:
               a3m = ''
               for i, c in enumerate(s):
                  if i in match_cols:
                     a3m += c          # uppercase residue or '-' deletion
                  elif c != '-':
                     a3m += c.lower()  # insertion column, skip gaps
               print(a3m, file=out)

   """
}
