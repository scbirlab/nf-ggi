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
       .to_csv("table.tsv")
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
   import pandas as pd
   
   files = glob("inputs/*.tsv")
   df = None
   for f in files:
      if df is None:
         df = pd.read_csv(f, sep="\\t")
         continue
      else:
         df = df.merge(pd.read_csv(f, sep="\\t"), how="outer")
   
   df.sort_values("ID").to_csv("table.tsv", sep="\\t", index=False)
   
   """
}
