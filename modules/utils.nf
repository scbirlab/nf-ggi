process Stack_tables {

   tag "->${dir}/${id}.${filename}.tsv"

   publishDir( 
      "${params.outputs}/ppi", 
      mode: 'copy',
      saveAs: { "${dir}/${id}.${filename}.tsv" },
   )

   input:
   tuple val( id ), path( tables, stageAs: "inputs/????.tsv" )
   val dir
   val filename

   output:
   tuple val( id ), path( "table.tsv" )

   script:
   """
   #!/usr/bin/env python
   from glob import glob
   import pandas as pd
   
   files = glob("inputs/*.tsv")
   (
       pd.concat([
           pd.read_csv(f, sep="\\t") for f in files
       ])
       .to_csv("table.tsv", index=False, sep="\\t")
   )
   
   """
}


process Join_tables {

   tag "->${dir}/${id}.${filename}.tsv"

   publishDir( 
      "${params.outputs}/ppi", 
      mode: 'copy',
      saveAs: { "${id}.${filename}.tsv" },
   )

   input:
   tuple val( id ), path( tables, stageAs: "inputs/????.tsv" )
   val dir
   val filename

   output:
   tuple val( id ), path( "table.tsv" )

   script:
   """
   #!/usr/bin/env python
   from glob import glob

   from carabiner import print_err
   import pandas as pd
   
   files = glob("inputs/*.tsv")
   df = None
   for f in files:
      print_err(f)
      df_right = pd.read_csv(f, sep="\\t")
      if df is None:
         df = df_right
         continue
      else:
         print_err("Common columns:", set(df.columns).intersection(df_right.columns))
         df = df.merge(df_right, how="outer")
   
   df.sort_values("ID").to_csv("table.tsv", sep="\\t", index=False)
   
   """
}
