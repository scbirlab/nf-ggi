process run_dca {

   label 'big_time'
   tag "${id}-${uniprot_id_bait}:${batch_idx}"
   stageInMode 'link'
   errorStrategy 'retry'
   maxRetries 2

   publishDir( 
      "${params.outputs}/dca", 
      mode: 'copy',
      saveAs: { "${id}-${uniprot_id_bait}-batch_${batch_idx}_dca.tsv" },
      pattern: "*.tsv"
   )
   publishDir( 
      "${params.outputs}/dca", 
      mode: 'copy',
      saveAs: { "${id}-${uniprot_id_bait}-batch_${batch_idx}_${it}" },
      pattern: "*-plot"
   )

   // Proteome ID, UniProtID, MSA file, [MSA file, ...]
   input:
   tuple val( id ), val( uniprot_id_bait ), path( msa1, stageAs: "bait/bait.a3m" ), val( batch_idx ), val( uniprot_id_queries ), path( msa2, stageAs: "queries/???.a3m" )
   val interspecies
   val plots

   output:
   tuple val( id ), path( "dca.tsv" ), emit: main
   path "*-plot" , emit: plots, optional: true

   script:
   def flag = ( interspecies ? "--interspecies" : "" )
   def plots_flag = ( plots ?  "--plot dca-plot" : "" )
   """
   MSA_LIST=msa-list.txt
   for item in queries/*.a3m
   do
      echo "\$item" >> \$MSA_LIST
   done
   yunta dca-single \
      <(echo "${msa1}") \
      --msa2 \$MSA_LIST \
      --list-file  ${flag} ${plots_flag} \
      --apc \
      --output "dca0.tsv"

   awk -F'\\t' -v OFS='\\t' 'NR == 1 { print "method", \$0 } NR > 1 { print "dca", \$0 }' \
      "dca0.tsv" \
   > "dca.tsv"
   """

   stub:
   """
   touch "dca.tsv"
   mkdir "dca-plot"
   touch "dca-plot/plot.png"
   echo "Skipping DCA for stub"
   """
}

process run_rf2track {

   label 'gpu_single'
   tag "${id}-${uniprot_id_bait}:${batch_idx}"
   stageInMode 'link'
   errorStrategy 'retry'
   maxRetries 2

   publishDir( 
      "${params.outputs}/rf2t", 
      mode: 'copy',
      saveAs: { "${id}-${uniprot_id_bait}-batch_${batch_idx}_rf2t.tsv" },
      pattern: "*.tsv"
   )
   publishDir( 
      "${params.outputs}/rf2t", 
      mode: 'copy',
      saveAs: { "${id}-${uniprot_id_bait}-batch_${batch_idx}_${it}" },
      pattern: "*-plot"
   )

   // Proteome ID, UniProtID, MSA file, [MSA file, ...]
   input:
   tuple val( id ), val( uniprot_id_bait ), path( msa1, stageAs: "bait/bait.a3m" ), val( batch_idx ), val( uniprot_id_queries ), path( msa2, stageAs: "queries/???.a3m" )
   val interspecies
   val plots

   output:
   tuple val( id ), path( "rf2t.tsv" ), emit: main
   path "rf2t-plot", emit: plots, optional: true

   script:
   def flag = (interspecies ? "--interspecies" : "")
   def plots_flag = (plots ?  "--plot rf2t-plot" : "" )
   """
   MSA_LIST=msa-list.txt
   for item in queries/*.a3m
   do
      echo "\$item" >> \$MSA_LIST
   done
   PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True CUDA_LAUNCH_BLOCKING=1 \
   yunta rf2t-single \
      <(echo "${msa1}") \
      --msa2 \$MSA_LIST \
      --list-file ${flag} ${plots_flag} \
      --output "rf2t0.tsv"

   awk -F'\\t' -v OFS='\\t' 'NR == 1 { print "method", \$0 } NR > 1 { print "rf2t", \$0 }' \
      "rf2t0.tsv" \
   > "rf2t.tsv"
   """

   stub:
   """
   touch "rf2t.tsv"
   mkdir "rf2t-plot"
   touch "rf2t-plot/plot.png"
   echo "Skipping RF2t for stub"
   """
}


process run_af2 {

   label 'gpu_single'
   tag "${id}-${uniprot_id_bait}:${batch_idx}"
   stageInMode 'link'
   errorStrategy 'retry'
   maxRetries 2

   publishDir( 
      "${params.outputs}/af2", 
      mode: 'copy',
      saveAs: { "${id}-${uniprot_id_bait}-batch_${batch_idx}_af2.pdb" },
      pattern: "af2/*.pdb"
   )
   publishDir( 
      "${params.outputs}/af2", 
      mode: 'copy',
      saveAs: { "${id}-${uniprot_id_bait}-batch_${batch_idx}_${it}" },
      pattern: "*.tsv"
   )
   publishDir( 
      "${params.outputs}/af2", 
      mode: 'copy',
      saveAs: { "${id}-${uniprot_id_bait}-batch_${batch_idx}-${it}" },
      pattern: "*-plot"
   )

   input:
   tuple val( id ), val( uniprot_id_bait ), path( msa1, stageAs: "bait/bait.a3m" ), val( batch_idx ), val( uniprot_id_queries ), path( msa2, stageAs: "queries/???.a3m" )
   val interspecies
   val plots

   output:
   tuple val( id ), path( "af2.tsv" ), emit: main
   path "af2/*.pdb", emit: pdb, optional: true
   path "*-plot", emit: plots, optional: true

   script:
   def flag = (interspecies ? "--interspecies" : "")
   def plots_flag = (plots ?  "--plot af2-plot" : "" )
   """
   MSA_LIST=msa-list.txt
   for item in queries/*.a3m
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
   XLA_PYTHON_CLIENT_MEM_FRACTION=.9 yunta af2-single \
      <(echo "${msa1}") \
      --msa2 \$MSA_LIST \
      --list-file ${flag} ${plots_flag} \
      --output "af2"

   awk -F'\\t' -v OFS='\\t' 'NR == 1 { print "method", \$0 } NR > 1 { print "af2", \$0 }' \
      "af2/_all_metrics.tsv" \
   > "af2.tsv"
   """

   stub:
   """
   mkdir "af2"
   touch "af2/stub.pdb"
   touch "af2.tsv"
   echo "Skipping AF2 for stub"
   """
}

process stack_table {

   tag "${id}-${filename}"

   publishDir( 
      "${params.outputs}/ppi", 
      mode: 'copy',
      saveAs: { "${id}-${filename}.tsv" },
   )

   input:
   tuple val( id ), path( tables, stageAs: "inputs/????.tsv" )
   val filename

   output:
   tuple val( id ), path( "table.tsv" )

   script:
   """
   tables=( inputs/*.tsv )
   head -n1 "\${tables[0]}" \
   | cat - <(tail -n+2 -q "\${tables[@]}") \
   > "table0.tsv"

   head -n1 "table0.tsv" \
   | cat - <(tail -n+2 "table0.tsv" | sort -k3,4 ) \
   > "table.tsv" \
   && rm "table0.tsv"
   """
}


process stack_table_py {

   tag "${id}-${filename}"

   publishDir( 
      "${params.outputs}/ppi", 
      mode: 'copy',
      saveAs: { "${id}-${filename}.tsv" },
   )

   input:
   tuple val( id ), path( tables, stageAs: "inputs/????.tsv" )
   val filename

   output:
   tuple val( id ), path( "table.tsv" )

   script:
   """
   python -c '
   from glob import glob
   import pandas as pd
   
   files = glob("inputs/*.tsv")
   df = pd.concat([pd.read_csv(f, sep="\\t") for f in files], axis=0)
   df.sort_values(["method", "ID"]).to_csv("table.tsv", sep="\\t", index=False)
   
   '
   """
}
