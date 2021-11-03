# Sequence based prediciton of Unstable and Gene-amplification Generated Heteroresistance

**Background:** Heteroresistance (HR) is a phenomenon in which a preexisting subpopulation of resistant cells can rapidly replicate in the presence of a given antibiotic, whereas the majority population of susceptible cells is killed. The mechanisms underlying HR are somewhat unclear, although unstable amplification of antibiotic resistance genes resulting in increased gene dosage is responsible for the resistant subpopulation in numerous cases

**Objective:** Predicting heteroresistance from bacterial genome data.

**Input:** 
- Nanopore and Illumina reads;
- "resistance labels"

**Metrics:**

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

## The project's structure:

The project home directory is `/home/andrei/Data/HetroR`
 which contains the following directories:

`data_raw/` - raw sequencing data (both Illumina and Nanopore reads) copied from Argos 

`data_filtered/` - sequencing data after filtering with filtlong and fastp

`qualcheck_reads/` - quality control data

`assemblies/` - hybrid assemblies made with Unicycler

`qualcheck_assembly/` - assembly quality control made with BUSCO and QUAST

`mapping/` - mapping of short reads onto genome assembly to collect unmapped reads

`plasmids/` - assembly of unmapped reads with SPAdes ('plasmid' mode) to assemble plasmids missed by Unicycler

`assemblies_joined/` - keeps merged hybrid and plasmid assembly

`annotations/` - assembly annotations made with PROKKA and tRNA-ScanSE

`resistance_genes/` - resistance genes tables identified by RGI (requires local CARD database)

`logs/` - tool's logs

`strain_lists/` - lists of available and processed strains, serves as input to some scripts that prepare data for the processing by the pipeline

`test_dir/` - various test of tools used in the pipeline

`tools/` - tools required by the pipeline but not available through conda package manager

`localDB/` - local instance of the CARD database (required by RGI)

`busco_downloads/` - files required by BUSCO

`coverage/` - a bunch of tables with Nanopore coverage

`final/` - "regulatory" dir created by Snakemake, nothing important there

`dags/` - directed acyclic graphs created by Snakemake representing the workflow (see below for the most recent example)

`notebooks/` - these are copies of actual notebooks that I keep on GoogleDrive. Hope these copies will be updated regularly

`workflow/` - the pipeline's actual code
 - `snakefile` - a file describing the workflow
 - `envs` - a set of YAML files describing required conda environments
 - `scripts` - additional scripts used by Snakemake and by me

`config.yaml` - a list of strains to be processed

## Workflow:

1. mount ARGOS
2. run `workflow/scripts/process_files.py` to transfer read files from ARGOS, rename them, calculate coverage and create config file.
3. load a local instance of CARD db (it must be in the project dir as 'localDB' - use `rgi load`)
4. run the pipeline using the command `snakemake --use-conda --cores 14 --resources mem_mb=12000`
5. run the following command to produce a nice heatmap of resistance hits in your strains:
   ```
   cd resistance genes; 
   for D in DA*; do ln -s "/home/andrei/Data/HeteroR/resistance_genes/"$D"/rgi_table.json" "/home/andrei/Data/HeteroR/resistance_genes/linked/"$D"_rgi_table.json"; done && 
   rgi heatmap -i linked -o heatmap -cat gene_family -clus samples
   ```

## Installation

The basic requiremnt is `snakemake`, install it using `conda` or `mamba`:

```
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

### SPADE installation

Right now the original version of SPADE can not be installed using conda and is not working properly.

I [forked](https://github.com/andrewgull/SPADE) the original SPADE [repo](https://github.com/yachielab/SPADE) to fix the code
and make it usable.

SPADE dependencies for the pipeline are described in `workflow/envs/spade-env.yaml`

My version of SPADE should be cloned and then executable files should be copied to
`~/minicinda3/envs/snakemake/bin` which is suboptimal but it's the only method available so far

## Dependencies

All dependencies are installed by `snakemake` itself in isolated environments using `conda`. 
The environments are described using YAML files that can be found in ``workflow/envs``

### List of tools used in the pipeline

1. Unicycler

## Current workflow's DAG

![dag](figures/dag_full.png)

## Heatmap example

the most recent version

The most recent version of RGI heatmap
![resistance genes heatmap](figures/heatmap54.png)

AMR genes categorised by AMR Gene Family and samples have been clustered hierarchically (see SciPy documentation). 
Yellow represents a perfect hit, teal represents a strict hit, purple represents no hit. 
Genes with asterisks (*) appear multiple times because they belong to more than one AMR Gene Family category in the antibiotic resistance ontology (ARO).