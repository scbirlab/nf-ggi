process run_dca {

   label 'big_time'
   tag "${id}-${uniprot_id_bait}:${batch_idx}"
   stageInMode 'link'
   errorStrategy 'retry'
   maxRetries 2

   publishDir( 
      "${params.outputs}/interactions/dca/tables", 
      mode: 'copy',
      saveAs: { "${id}-${uniprot_id_bait}-batch_${batch_idx}.${it}" },
      pattern: "*.tsv"
   )
   publishDir( 
      "${params.outputs}/interactions/dca/plots", 
      mode: 'copy',
      saveAs: { "${id}-${uniprot_id_bait}-batch_${batch_idx}.${it}" },
      pattern: "plot-*"
   )

   // Proteome ID, UniProtID, MSA file, [MSA file, ...]
   input:
   tuple val( id ), val( uniprot_id_bait ), path( msa1, stageAs: "bait/bait.a3m" ), val( batch_idx ), val( uniprot_id_queries ), path( msa2, stageAs: "queries/???.a3m" )
   val interspecies
   val plots

   output:
   tuple val( id ), path( "dca.tsv" ), emit: main
   path "plot-dca" , emit: plots, optional: true

   script:
   """
   MSA_LIST=msa-list.txt
   for item in queries/*.a3m
   do
      echo "\$item" >> \$MSA_LIST
   done

   export MPLCONFIGDIR="\${PWD}/.mplconfig"
   export XDG_CACHE_HOME="\${PWD}/.cache"
   yunta dca-single \
      <(echo "${msa1}") \
      --msa2 \$MSA_LIST \
      --list-file  ${interspecies ? "--interspecies" : ""} ${plots ?  "--plot plot-dca" : ""} \
      --apc \
      --output "dca.tsv"
   
   """

   stub:
   """
   printf 'ID\\ndca\\tA-B\\n' > "dca.tsv"
   mkdir "plot-dca"
   touch "plot-dca/plot.png"
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
      "${params.outputs}/interactions/rf2t/tables", 
      mode: 'copy',
      saveAs: { "${id}-${uniprot_id_bait}-batch_${batch_idx}.${it}" },
      pattern: "*.tsv"
   )
   publishDir( 
      "${params.outputs}/interactions/rf2t/plots", 
      mode: 'copy',
      saveAs: { "${id}-${uniprot_id_bait}-batch_${batch_idx}.${it}" },
      pattern: "plot-*"
   )

   // Proteome ID, UniProtID, MSA file, [MSA file, ...]
   input:
   tuple val( id ), val( uniprot_id_bait ), path( msa1, stageAs: "bait/bait.a3m" ), val( batch_idx ), val( uniprot_id_queries ), path( msa2, stageAs: "queries/???.a3m" )
   val interspecies
   val plots

   output:
   tuple val( id ), path( "rf2t.tsv" ), emit: main
   path "plot-*", emit: plots, optional: true

   script:
   """
   MSA_LIST=msa-list.txt
   for item in queries/*.a3m
   do
      echo "\$item" >> \$MSA_LIST
   done

   export MPLCONFIGDIR="\${PWD}/.mplconfig"
   export XDG_CACHE_HOME="\${PWD}/.cache"
   PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True CUDA_LAUNCH_BLOCKING=1 \
   yunta rf2t-single \
      <(echo "${msa1}") \
      --msa2 \$MSA_LIST \
      --list-file ${interspecies ? "--interspecies" : ""} ${plots ?  "--plot plot-rf2t" : ""} \
      --output "rf2t.tsv"

   """

   stub:
   """
   printf 'ID\\nrf2t\\tA-B\\n' > "rf2t.tsv"
   mkdir "plot-rf2t"
   touch "plot-rf2t/plot.png"
   echo "Skipping RF2t for stub"
   """
}


process run_af2 {

   label 'gpu_single'
   tag "${id}-${uniprot_id_bait}:${batch_idx}"
   stageInMode 'link'
   // errorStrategy 'retry'
   // maxRetries 2

   publishDir( 
      "${params.outputs}/interactions/af2/pdb", 
      mode: 'copy',
      saveAs: { "${id}-${uniprot_id_bait}-batch_${batch_idx}.${it}" },
      pattern: "af2/*.pdb"
   )
   publishDir( 
      "${params.outputs}/interactions/af2/tables", 
      mode: 'copy',
      saveAs: { "${id}-${uniprot_id_bait}-batch_${batch_idx}.${it}" },
      pattern: "*.tsv"
   )
   publishDir( 
      "${params.outputs}/interactions/af2/plots", 
      mode: 'copy',
      saveAs: { "${id}-${uniprot_id_bait}-batch_${batch_idx}.${it}" },
      pattern: "plot-*"
   )

   input:
   tuple val( id ), val( uniprot_id_bait ), path( msa1, stageAs: "bait/bait.a3m" ), val( batch_idx ), val( uniprot_id_queries ), path( msa2, stageAs: "queries/???.a3m" )
   val interspecies
   val plots

   output:
   tuple val( id ), path( "af2.tsv" ), emit: main
   path "af2/*.pdb", emit: pdb, optional: true
   path "plot-*", emit: plots, optional: true

   script:
   """
   set -euox pipefail

   MSA_LIST=msa-list.txt
   for item in queries/*.a3m
   do
      echo "\$item" >> \$MSA_LIST
   done

   CUDNN_PY="\$(python - <<'PY'
   try:
      import nvidia.cudnn, os
      print(nvidia.cudnn.__file__ or "")
   except Exception:
      print("")
   PY
   )"
   if [ -n "\$CUDNN_PY" ]; then
      CUDNN_PATH="\$(dirname "\$CUDNN_PY")"
      export LD_LIBRARY_PATH="\${CUDNN_PATH}/lib:\${LD_LIBRARY_PATH:-}"
      >&2 echo "CuDNN at \$CUDNN_PATH"
   else
      >&2 echo "[cuDNN python pkg not found; skipping LD_LIBRARY_PATH edit]"
   fi

   nvcc --version 2>/dev/null || echo "[nvcc not found]" && (>&2 echo "\$(nvcc --version)")

   if command -v nvcc >/dev/null 2>&1; then
      export CUDA_HOME="\$(dirname "\$(dirname "\$(command -v nvcc)")")"   # e.g., /usr/local/cuda
      export JAX_CUDA_PATH="\$CUDA_HOME"
      export XLA_FLAGS="--xla_gpu_cuda_data_dir=\$CUDA_HOME \${XLA_FLAGS:-}"
      >&2 echo "CUDA_HOME=\$CUDA_HOME"
   fi
   
   python3 - <<'PY' || true
   try:
      import jax, os
      print("jax:", jax.__version__)
      print("CUDA_VISIBLE_DEVICES:", os.environ.get("CUDA_VISIBLE_DEVICES"))
      print("Backend:", jax.lib.xla_bridge.get_backend().platform)
   except Exception as e:
      print(f"[JAX probe skipped: {e}]")
   try:
      import tensorflow as tf
      print(f"Available devices:\\n{tf.config.list_physical_devices()}")
   except Exception as e:
      print(f"[TF probe skipped due to import error: {e}]")
   PY

   export LD_LIBRARY_PATH="/opt/conda/envs/env/lib:/usr/local/cuda/lib64:/usr/lib/x86_64-linux-gnu:\${LD_LIBRARY_PATH}"
   export MPLCONFIGDIR="\${PWD}/.mplconfig"
   export XDG_CACHE_HOME="\${PWD}/.cache"
   mkdir -p "\$MPLCONFIGDIR" "\$XDG_CACHE_HOME"

   XLA_PYTHON_CLIENT_MEM_FRACTION=.9 yunta af2-single \
      <(echo "${msa1}") \
      --msa2 \$MSA_LIST \
      --list-file ${interspecies ? "--interspecies" : ""} ${plots ?  "--plot plot-af2" : ""} \
      --output "af2.tsv"

   """

   stub:
   """
   mkdir "plot-af2" "af2"
   touch "af2/stub.pdb"
   printf 'ID\\naf2\\tA-B\\n' > "af2.tsv"
   touch plot-af2/plot.png
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
      saveAs: { "${id}.${filename}.tsv" },
   )

   input:
   tuple val( id ), path( tables, stageAs: "inputs/????.tsv" )
   val filename

   output:
   tuple val( id ), path( "table.tsv" )

   script:
   """
   #/usr/bin/env python
   from glob import glob
   import pandas as pd
   
   files = glob("inputs/*.tsv")
   df = pd.concat([pd.read_csv(f, sep="\\t") for f in files], axis=0)
   (
      df
      .pivot(
         columns=
      )
   )
   df.sort_values(["method", "ID"]).to_csv("table.tsv", sep="\\t", index=False)
   
   '
   """
}
