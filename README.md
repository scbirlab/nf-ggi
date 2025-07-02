# Gene-gene interaction screening pipeline

![GitHub Workflow Status (with branch)](https://img.shields.io/github/actions/workflow/status/scbirlab/nf-ggi/nf-test.yml)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.10.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

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

**If you're at the Crick, these databases already reside on NEMO, and there is no need to downlaod them.**

### Software

You need to have Nextflow and either Singularity, Docker, of conda installed on your system.

#### First time using Nextflow?

##### Crick users

If you're at the Crick **or your shared cluster has Nextflow and Singularity already installed**, try:

```bash
module load Nextflow Singularity
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

There are three run modes for the pipeline:

- "self": run all protein-protein interactions within an organism
- "bait": run interactions between all proteins from an organism and either one protein or another organism's proteome
- "custom": run specified protein pairs from a file

The easiest way to get going is by specifying parameters on the command-line:

```bash
bfd=path/to/your/bfd
uniclust=path/to/your/uniclust
nextflow run scbirlab/nf-ggi \
    --bfd "$bfd" --uniclust "$uniclust" \
    --organism_id 243273 \
    --dca --rf2t  --plots
```

Here's what the flags mean:
- `--organism_id`: The Taxon ID of the organism, whih you can find at NCBI or UniProt
- `--dca`, `--rf2t`: Run direct-coupling analysis and RosettaFold-2track
- `--plots`: Generate amino acid contact maps. This takes about 10MB per protein-protein interaction, so be sure you have enough disk space!

You could also run `--metabolites` to get the metabolic network, and `--string` to get the STRING co-expression network.

Because only `--organism_id` is provided, the pipeline assumes "self" mode (i.e. all-vs-all within taxon 243273).

You can run bait mode by providing `--bait <UniProt ID>`:

```bash
nextflow run scbirlab/nf-ggi --bfd "$bfd" --uniclust "$uniclust" \
    --organism_id 559292 --bait P00931 \
    --dca
```

The bait can be another organisms's proteome. In this case, we need to specifiy that `--bait_is_taxon` and `--interspecies`:

```bash
nextflow run scbirlab/nf-ggi --bfd "$bfd" --uniclust "$uniclust" \
    --organism_id 559292 --bait 1773 \
    --bait_is_taxon --interspecies \
    --rf2t
```

### Running more than one query in parallel

Make a [sample sheet (see below)](#sample-sheet) with columns representing the flags above, and, optionally, a [`nextflow.config` file](#inputs) in the 
directory where you want the pipeline to run. Then simply run:

```bash 
nextflow run scbirlab/nf-ggi
```

### Pipeline versions

Each time you run the pipeline after the first time, Nextflow will use a locally-cached version which 
will not be automatically updated. If you want to ensure that you're using the very latest version of the 
pipeline, use the `-latest` flag.

```bash 
nextflow run scbirlab/nf-ggi -latest
```

If you want to run a particular tagged version of the pipeline, such as `v0.0.2`, you can do so using

```bash 
nextflow run scbirlab/nf-ggi -r v0.0.2
```

For help, use `nextflow run scbirlab/nf-ggi --help`.

The first time you run the pipeline on your system, the software dependencies in `environment.yml` will be installed. 
This may take several minutes.

## Inputs

### Command-line usage

The pipeline can be run with command-line arguments:

```bash
# intra-species all-vs-all:
nextflow run scbirlab/nf-ggi --uniclust <path> --bfd <path> \
    --organism_id <taxon ID>
# intra-species all-vs-1:
nextflow run scbirlab/nf-ggi --uniclust <path> --bfd <path> \
    --organism_id <taxon ID> --bait <UniProtID>
# inter-species all-vs-all:
nextflow run scbirlab/nf-ggi --uniclust <path> --bfd <path> \
    --organism_id <taxon ID> --bait <taxon ID> \
    --bait_is_taxon --interspecies
# custom list of pairs:
nextflow run scbirlab/nf-ggi --uniclust <path> --bfd <path> \
    --organism_id <taxon ID> \
    --filename <path> --column1 <gene-col1> --column2 <gene-col2> \
    [--interspecies --organism_id2 <taxon ID>] [--format <gene-name-type>]
```

The following parameters are **required**:

```
--organism_id             Taxon ID for organism
Bait mode:
    --bait                 UniProt ID for bait protein, or Taxon ID for bait organism
Custom mode:
    --filename             Filename to get custom protein pairs
    --column1, --column2   Column names from --filename to get protein IDs
```

The following parameters are **optional**. They have default values which can
 be overridden if necessary.

 ```bash
 --reviewed      Only pull SwissProt reviewed proteins from proteome
--isoforms       Additionally pull isoform sequences from proteome
--proteome_opts  Additonal filters for pulling from proteome. Check 
            https://www.ebi.ac.uk/proteins/api/doc/#!/proteins/search for options.
--bait_is_taxon  Indicate that bait is an organism ID
--interspecies   Run analysis between interacting species proteomes
--plots          Generate contact map plots
--organism_id2   When providing a file of pairs, if the second protein (--column2) is from another organism than the first
--format         Type of gene identifier in --column1, --column2. Default: "Gene_Name"
--test           Whether to run in test mode. Default: false.
--outputs        Output folder. Default: "outputs".
--batch_size     What size to batch protein-protein interactions into. Default: 100.
```

### Sample sheet

You can run multiple combinations in one command using a sample sheet. The sample sheet is a **CSV** file with one row per combination of parameters to run. The column headings have the same names as the required flags for command-line usage. The optional flags are still on the command line, and applied to everything in the run. Here, the mode needs to be specified with `--mode`:

```bash
nextflow run scbirlab/nf-ggi --mode self --rf2t --metabolites --sample_sheet path/to/sample-sheet.csv
```

#### Sample sheet structure

Here is an example of the sample sheet for `mode = "self"`, to find all the mycoplasma protein-protein interactions:

| organism_id | proteome_name           |
| ----------- | ----------------------- |
| 243273      | "Mycoplasma genitalium" |

The `proteome_name` column is not neccesary, but you can add extra columns with human-readable annotations for your own sanity. We recommend:
- `proteome_name`
- (if using a bait) `bait_name`

If running with `mode = "bait"`, to do a pulldown against a single bait protein, add another column with the bait UniProt ID.

| organism_id | proteome_name         | bait   | bait_name |
| ----------- | --------------------- | ------ | --------- |
| 243273      | Mycoplasma genitalium | P47259 | FolD      |


If running with `mode = "custom"`, to do a pulldown against a single bait protein, add another column with the bait UniProt ID.

| organism_id | proteome_name            | format    | filename  | column1    | column2   |
| ----------- | ------------------------ | --------- | --------- | ---------- | --------- |
| 559292      | Saccharomyces cerevisiae | Gene_Name | combos.csv | query_orf | array_orf |

In this case, `combos.csv` must be in the `inputs` folder defined above. It would look like:

| query_orf	| query_gene_name | array_orf | array_gene_name |
| --------- | --------------- | --------- | --------------- |
| YAL058W	| CNE1            | YAL068C   |	PAU8            |

Further examples are in the `test` directory of this repository.

### Config-file usage (recommended)

For reproducibility, self-documentation, and to save typing, parameters with the same names as the command line flags above can be provided in a `nextflow.config` file in the working directory. For example:

```groovy
params {
    organism_id = "243273"
    dca = true
    rf2t = true
}
```

Or with a sample sheet:

```groovy
params {
    sample_sheet = "path/to/sample-sheet.csv"
    mode = "self"
    dca = true
    rf2t = true
    metabolites = true
    string = true
}
```

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
