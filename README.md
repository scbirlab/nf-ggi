# Gene-gene interaction screening pipeline

![GitHub Workflow Status (with branch)](https://img.shields.io/github/actions/workflow/status/scbirlab/nf-ggi/nf-test.yml)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.10.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

**scbirlab/nf-ggi** is a Nextflow pipeline to screen gene-gene interactions within an organism or between organisms 
(in the case of host-pathogen or phage-bacterium interactions).

**Table of contents**

- [Processing steps](#processing-steps)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Credit](#credit)
- [Issues, problems, suggestions](#issues-problems-suggestions)
- [Further help](#further-help)

## Processing steps

1. Download Rhea DB (of metabolites) in preparation for searching.

For proteins or proteomes in the [sample sheet](#sample-sheet):

1. Download its STRING database and tidy up the data.
2. Download FASTA sequences of proteins from UniProt
    - If multiple proteomes are available, choose according to this priority: "Reference and representative", "Reference", "Representative", "Other"
3. Find reactions in Rhea DB and connect products with reactants between enzymes in the proteome.

For each FASTA sequence:

4. Generate a multiple sequence alignment with `hhblits`.

For `method == "self"`:

5. Within each organism, generate all unique pairs of proteins.

For `method == "bait"`:

5. All unique pairs of proteins between the organism and listed baits.

For `method == "custom"`:

5. All unique pairs of proteins listed.

Then for each protein pair, optionally:

6. Calculate the co-evolutionary signal with DCA, optionally generating plots of contact maps.
7. Predict the interface contact map with `yunta rf2t` (RosettaFold-2track), optionally generating plots of contact maps.
8. Predict the protein-protein complex structure map with `yunta af2` (AlphaFold2), optionally generating plots of contact maps.

## Requirements

You need access to the UniClust and BFD databases, and you need Nextflow and conda to be installed.

### Databases

To generate multiple-sequence alignments (MSAs) for co-evolutionary analysis, `hhblits` databases of 
pre-clustered sequences is required. Unfortunately, these are extremely large, so cannot be downlaoded as 
part of the pipeline. You should download the [UniClust](https://uniclust.mmseqs.com/) and [BFD](https://bfd.mmseqs.com/) 
databases, then set the `--uniclust` and `--bfd` parameters of the pipeline ([see below](#inputs)).

If you're at the Crick, these databases already reside on NEMO, and there is no need to downlaod them.

### Software

You need to have Nextflow and `conda` installed on your system.

#### First time using Nextflow?

##### Crick users

If you're at the Crick **or your shared cluster has it already installed**, try:

```bash
module load Nextflow
```

##### Others

Otherwise, if it's your first time using Nextflow on your system, you can install it using `conda`:

```bash
conda install -c bioconda nextflow 
```

You may need to set the `NXF_HOME` environment variable. For example,

```bash
mkdir -p ~/.nextflow
export NXF_HOME=~/.nextflow
```

To make this a permanent change, you can do something like the following:

```bash
mkdir -p ~/.nextflow
echo "export NXF_HOME=~/.nextflow" >> ~/.bash_profile
source ~/.bash_profile
```

## Quick start

Make a [sample sheet (see below)](#sample-sheet) and, optionally, a [`nextflow.config` file](#inputs) in the 
directory where you want the pipeline to run. Then run Nextflow.

```bash 
nextflow run scbirlab/nf-ggi
```

Each time you run the pipeline after the first time, Nextflow will use a locally-cached version which 
will not be automatically updated. If you want to ensure that you're using the very latest version of the 
pipeline, use the `-latest` flag.

```bash 
nextflow run scbirlab/nf-ggi -latest
```

If you want to run a particular tagged version of the pipeline, such as `v0.0.3`, you can do so using

```bash 
nextflow run scbirlab/nf-ggi -r v0.0.3
```

For help, use `nextflow run scbirlab/nf-ggi --help`.

The first time you run the pipeline on your system, the software dependencies in `environment.yml` will be installed. 
This may take several minutes.

## Inputs

The following parameters are **required**:

- `sample_sheet`: filename of CSV with information about the samples and FASTQ files to be processed. Must be in the `inputs` folder (see below).
- `uniclust`: Path to `hhblits` UniClust database. This is very large, so you need to have it already downlaoded on your system.
- `bfd`: Path to `hhblits` BFD database. This is very large, so you need to have it already downlaoded on your system.

The following parameters are **optional**. They have default values which can
 be overridden if necessary.

 - `mode = "self"`: Which mode to run ([see below](#sample-sheet)). Options are: "self" (all vs all), "bait" (all vs some), "custom" (some-vs-some). 
 - `interspecies = false`: Whether the proteins come from two different species.
 - `bait_is_taxon = false`: In "bait" mode, the baits can be a taxon ID instead of individual proteins, in which case the pipeline will fetch the proteome for the bait taxon.
 - `inputs = "inputs"`: Folder to look for sample sheet and any other inputs
 - `outputs = "outputs"`: Folder to put outputs from the pipeline
 - `batch_size = 100`: How many protein-protein interactions to group into one job at a time.
 - `test = false`: Whether to run in test mode. If so, only 3 proteins per organism will be analyzed.
 - `rhea_url = "https://ftp.expasy.org/databases/rhea"`: URL to download Rhea reaction database

The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.

Here is an example of the `nextflow.config` file:

```nextflow
params {
    sample_sheet = "/path/to/sample-sheet.csv"
}
```

Alternatively, you can provide the parameters on the command line:

```bash
nextflow run scbirlab/nf-ggi --sample_sheet /path/to/sample-sheet.csv 
``` 

### Sample sheet

The sample sheet is a **CSV** file indicating which organisms you want to analyze.

The file must have a header with **required** column names below, and one line per combination to be processed.

- `organism_id`: the NCBI Taxonomic IDs for your organisms. This can be found at [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy)

In `mode = "bait"`, additionally **required**:
- `bait`: Uniprot ID if `bait_is_taxon`, otherwise taxon ID

In `mode = "custom"`, additionally **required**:
- `filename`: The CSV or TSV file to get the custom combinations from
- `column1`, `column2`: The names of the columns containing protein identifiers
- `format`: The type of protein identifiers, which will be looked up on UniProt. `Gene_Name` works well and felxibly for most commonly used gene names; currently `UniprotID` doesn't work. 


You can add extra columns with human-readable annotations for your own sanity. 
We recommend:
- `proteome_name`
- (if using a bait) `bait_name`

Here is an example of the sample sheet for `mode = "self"`, to find all the mycoplasma protein-protein interactions:

| organism_id | proteome_name           |
| ----------- | ----------------------- |
| 243273      | "Mycoplasma genitalium" |

If running with `mode = "bait"`, to do a pulldown against a single bait protein, add another column with the bait UniProt ID.

| organism_id | proteome_name           | bait   | bait_name |
| ----------- | ----------------------- | ------ | --------- |
| 243273      | "Mycoplasma genitalium" | P47259 | FolD      |


If running with `mode = "custom"`, to do a pulldown against a single bait protein, add another column with the bait UniProt ID.

| organism_id | proteome_name              | format    | filename  | column1 | column2 |
| ----------- | -------------------------- | --------- | --------- | -------- | -------- |
| 559292      | "Saccharomyces cerevisiae" | Gene_Name | combos.csv | query_orf | array_orf |

In this case, `combos.csv` must be in the `inputs` folder defined above. It would look like:

| query_orf	| query_gene_name | array_orf | array_gene_name |
| --------- | --------------- | --------- | --------------- |
| YAL058W	| CNE1            | YAL068C   |	PAU8            |

Further examples are in the `test` directory of this repository.

## Outputs

Outputs are saved in the `output` folder defined above. They include these directories:

- `string`: STRING co-expression values
- `metabolites`: Reconstructed metabolic network
- `msa`: All MSA files
- `ppi`: All protein-protein interaction data
- `sequences`: Protein sequences

## Credit

The idea of using DCA, [RoseTTAFold](https://github.com/RosettaCommons/RoseTTAFold)-2track, and [AlphaFold2](https://github.com/google-deepmind/alphafold) in a cascade of increasingly expensive and specific PPI detection methods has been explored in a series of papers from David Baker's lab:

- [Cong et al., Protein interaction networks revealed by proteome coevolution. _Science_, 2019](https://doi.org/10.1126/science.aaw6718)
- [Humpreys et al., Computed structures of core eukaryotic protein complexes. _Science_, 2021](https://doi.org/10.1126/science.abm4805)
- [Humpreys et al., Protein interactions in human pathogens revealed through deep learning. _Nature Microbiology_, 2024](https://doi.org/10.1038/s41564-024-01791-x)

`scbirlab/nf-ggi` applies these algorithms in a Nextflow pipeline to allow easy scaling, and enables inter-species interactions. It also reconstructs metabolic networks, and pulls known interactions from the STRING database.

## Issues, problems, suggestions

Add to the [issue tracker](https://www.github.com/scbirlab/nf-ggi/issues).

## Further help

Here are the pages of the software and databases used by this pipeline.

Databases:

- [STRING](https://string-db.org/) for co-expression
- [Rhea](https://www.rhea-db.org/) for enzyme reactions
- [UniProt](https://www.uniprot.org/) for protein sequences
- [NCBI Genbank](https://www.ncbi.nlm.nih.gov/genbank/) for taxonomy

Software:

- [hhblits](https://github.com/soedinglab/hh-suite) for generating MSAs
- [rdkit](https://www.rdkit.org/docs/index.html) for cheminformatics of enzyme reactants and products
- [yunta](https://www.github.com/scbirlab/yunta) for running DCA, RosettaFold-2track, and AlphaFold2 on MSAs