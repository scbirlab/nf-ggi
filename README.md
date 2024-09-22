# Gene-gene interaction screening pipeline

![GitHub Workflow Status (with branch)](https://img.shields.io/github/actions/workflow/status/scbirlab/nf-ggi/nf-test.yml)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.10.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

**scbirlab/nf-ggi** is a Nextflow pipeline to screen gene-gene interactions from a specified organism (or set of organisms). 

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

1. Download Rhea DB in preparation for searching.

For each organism ID provided:

1. Download its STRING database and tidy data.
2. Download FASTA sequences of its proteins from UniProt
    - Where possible, reference proteomes are used.
3. Find reactions in Rhea DB and connect products with reactants between enzymes in the proteome.

For each FASTA sequence in each organism:

4. Generate a multiple sequence alignment with `hhblits`.

Then within each organism:

5. Generate all unique pairs of proteins.

Then for each protein pair:

6. Calculate the co-evolutionary signal with DCA.
7. Predict the interface contact map with `yunta rf2t` (RosettaFold-2track).
8. Predict the protein-protein complex structure map with `yunta af2` (AlphaFold2).

## Requirements

### Databases

To generate multiple-sequence alignments (MSAs) for co-evolutionary analysis, `hhblits` databases of pre-clustered sequences is required. Unfortunately, these are extremely large, so cannot be downlaoded as part of the pipeline. You should download the [UniClust](https://uniclust.mmseqs.com/) and [BFD](https://bfd.mmseqs.com/) databases, then set the `--uniclust` and `--bfd` parameters of the pipeline ([see below](#inputs)).

If you're at the Crick, these databases already reside on NEMO.

### Software

You need to have Nextflow and `conda` installed on your system.

#### First time using Nextflow?

If you're at the Crick or your shared cluster has it already installed, try:

```bash
module load Nextflow
```

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

Make a [sample sheet (see below)](#sample-sheet) and, optionally, a [`nextflow.config` file](#inputs) in the directory where you want the pipeline to run. Then run Nextflow.

```bash 
nextflow run scbirlab/nf-ont-call-variants
```

Each time you run the pipeline after the first time, Nextflow will use a locally-cached version which will not be automatically updated. If you want to ensure that you're using the very latest version of the pipeline, use the `-latest` flag.

```bash 
nextflow run scbirlab/nf-ont-call-variants -latest
```
If you want to run a particular tagged version of the pipeline, such as `v0.0.1`, you can do so using

```bash 
nextflow run scbirlab/nf-ont-call-variants -r v0.0.2
```

For help, use `nextflow run scbirlab/nf-ont-call-variants --help`.

The first time you run the pipeline on your system, the software dependencies in `environment.yml` will be installed. This may take several minutes.

## Inputs

The following parameters are required:

- `sample_sheet`: path to a CSV with information about the samples and FASTQ files to be processed
- `uniclust`: Path to `hhblits` UniClust database. This is very large, so you need to have it already downlaoded on your system.
- `bfd`: Path to `hhblits` BFD database. This is very large, so you need to have it already downlaoded on your system.

The following parameters have default values which can be overridden if necessary.

 - `rhea_url = "https://ftp.expasy.org/databases/rhea"`: URL to download Rhea reaction database
 - `outputs = "outputs"`: Output folder
 - `batch_size = 100`: How many protein-protein interactions to group into one job at a time.
 - `test = false`: Whether to run in test mode. If so, only 3 proteins per organism will be analyzed.
 - `non_self = false`: Whether to run in non-self mode. This is where a whole proteome is run against a single bait protein (rather than all pairwise from the proteome).

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

The file must have a header with the column names below, and one line per organism to be processed.

- `organism_id`: the NCBI Taxonomic ID for your organism. This can be found at [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy)
- `proteome_name`: This can be anything, but should be an informative description

Here is an example of the sample sheet:

| organism_id | proteome_name           |
| ----------- | ----------------------- |
| 243273      | "Mycoplasma genitalium" |

If running with `--non-self`, to do a pulldown against a single bait protein, add another column with the bait UniProt ID.

| organism_id | proteome_name           | bait   |
| ----------- | ----------------------- | ------ |
| 243273      | "Mycoplasma genitalium" | P47259 |

## Outputs

Outputs are saved in the same directory as `sample_sheet`. They are organised under three directories:

- `coexpression`: STRING co-expression values
- `metabolites`: Reconstructed metabolic network
- `msa`: All MSA files
- `ppi`: All protein-protein interaction data
- `sequences`: Each organism's proteome sequence

## Credit

The idea of using DCA, [RoseTTAFold](https://github.com/RosettaCommons/RoseTTAFold)-2track, and [AlphaFold2](https://github.com/google-deepmind/alphafold) in a cascade of increasingly expensive and specific PPI detection methods has been explored in a series of papers from David Baker's lab:

- [Cong et al., Protein interaction networks revealed by proteome coevolution. _Science_, 2019](https://doi.org/10.1126/science.aaw6718)
- [Humpreys et al., Computed structures of core eukaryotic protein complexes. _Science_, 2021](https://doi.org/10.1126/science.abm4805)
- [Humpreys et al., Protein interactions in human pathogens revealed through deep learning. _Nature Microbiology_, 2024](https://doi.org/10.1038/s41564-024-01791-x)

`scbirlab/nf-ggi` applies these algorithms in a Nextflow pipeline to allow easy scaling. It also reconstructs metabolic networks, and pulls known interactions from the STRING database.

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