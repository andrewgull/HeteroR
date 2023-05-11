# Sequence based prediciton of Unstable and Gene-amplification Generated Heteroresistance

![hr_model](images/unstable_HR_model.png)

**Background:** Heteroresistance (HR) is a phenomenon in which a preexisting subpopulation of resistant cells can rapidly replicate in the presence of a given antibiotic, whereas the majority population of susceptible cells is killed. The mechanisms underlying HR are somewhat unclear, although unstable amplification of antibiotic resistance genes resulting in increased gene dosage is responsible for the resistant subpopulation in numerous cases

**Objective:** Predicting heteroresistance from bacterial genome data.

**Input:** 
- Nanopore and Illumina reads;
- "resistance labels"

**Metrics:**

![metrics](images/HR_workflow_features_scheme.png)

1. Presence of resistance genes (RG)
   - presence of known RG as identified from the CARD database
   - presence of efflux pumps as identified from the CARD database

2. Presence of amplifiable regions
   - For each RG all pairs of direct repeats (DR) flanking the RG. 
     - DR min length = 10 bp, max mismatch = 10%, search range 100 kb. 
     - DR pairs are scored according to their length, their level of identity, and their distance to each other. 
     - These three parameters should reflect the probability that the segment encompassing the RG will be amplified. 

3. Presence of deleterious effects
   - truncated genes (any DR inside a gene?); 
   - broken operons (any DR inside an operon?); operons will be identified via [ODB4](https://operondb.jp/)
   - co-expression as identified by [StringDB](https://string-db.org/) (Are any of the genes in the amplified segment known to be co-expressed with a gene outside of the segment?)
   - toxic essential genes: essential genes (as identified by databases TraDIS (?) and Keio(?)) being truncated may become toxic.
   - toxic if over-expressed: gene amplifications may increase the expression of the gene, which in turn may be toxic. Related data can be found in databases [EDGE](https://www.pnas.org/content/pnas/early/2013/11/05/1312361110), [ASKA](https://academic.oup.com/dnaresearch/article/12/5/291/350187) and [PandaTox](https://exploration.weizmann.ac.il/pandatox/1_0/home.html).

The pipeline is created using [Snakemake](https://snakemake.readthedocs.io/en/stable) - a Python-based workflow management system for reproducible and scalable data analysis. [The "rolling" paper reference](https://f1000research.com/articles/10-33/v2) 

### Directories description

`images/` - workflow DAGs, RGI heatmaps etc.

`resources/` is for storing retrieved/transferred data like: raw_reads, BUSCO downloads, strain lists

   - `data_raw/` - raw sequencing data (both Illumina and Nanopore reads) transferred from Argos 
   - `strain_lists/` - lists of available and processed strains, serves as input to some scripts that prepare data for the processing by the pipeline
   - `busco_downloads/` - files required by BUSCO

`results/` is for everything the pipeline produces
  
`localDB/` - local instance of the CARD database (required by RGI)

`notebooks/` - these are copies of actual notebooks that I keep on GoogleDrive. Hope these copies will be updated regularly

`workflow/` - the pipeline's actual code
 - `snakefile` - a file describing the workflow
 - `envs/` - a set of YAML files describing required conda environments
 - `scripts/` - additional scripts used by Snakemake and by me

`config.yaml` - a list of strains to be processed

## How to run the pipeline:

### Prerequisites

1. mount ARGOS
2. load a local instance of CARD db (it must be in the project dir as 'localDB' - use `rgi load`)
3. download BUSCO data base

### Steps

1. get list of strains on ARGOS `ls ~/Data/Argos/imb_sal_raw/Sequenced_reference_strains/Sequencing/Strains/ > strains_on_argos.txt`
2. get list of strains to process `bash workflow/scripts/get_new_strains_list.sh strains_on_argos.txt > strains.txt`
3. run `workflow/scripts/process_files.py -s resources/strain_lists/strains.txt -c configs/strains.yaml -l 45000000` to transfer read files from ARGOS, rename them, calculate coverage and create config file for snakemake. 
4. run the pipeline using the command `snakemake --snakefile workflow/snakefile.smk --configfile configs/strains.yaml --use-conda --cores 10`

Steps with species-specific parameters:
- PROKKA (genus, species)
- QC_assembly.py (taxonomic dataset, to find available BUSCO datasets run busco --list-datasets)
- trim_nanopore.py (parameter genlen)

NB. It might be useful to change memory limitations in FastQC script which is used by quality check steps of the pipeline.
Current settings are: 

`memory = 1250 * threads, stack = 1000 * threads; if no threads specified, max memory 10000m, min 5000m`

## Installation

The basic requirement is `snakemake`, install it using `conda` or `mamba`:

```
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

Additionally install FastQC, path to the fastqc executable file must be provided to quality check scripts in the pipeline.

### RGI Database installation

RGI requires its database to be in the project's directory and named `localDB`. 
No other way to set the database location is available.
That's why `localDB` must be in the project's directory even though it is not listed here.

### BUSCO installation

BUSCO tools and databases (needed for searching BUSCO genes) -- works in Linux only!
Download the database using the command below
```
quast-download-busco
```

## Dependencies

All dependencies are installed by `snakemake` itself in isolated environments using `mamba`. 
The environments are described using YAML files that can be found in ``workflow/envs``

### List of tools used in the pipeline

1. FastQC
2. fastp
3. filtlong
4. Unicycler
5. BUSCO
6. QUAST
7. BWA
8. SAMtools
9. SPAdes
10. PROKKA
11. tRANscan-SE
12. RGI
13. BEDtools
14. GRF

## Current workflow's DAG

![dag](images/dag.png)

## Data analysis

All the data analysis code including machine learning part can be found in ``notebooks/``