# HeteroR
code for heteroresistance project

Example DAG

![dag](figures/dag_full.png)

project structure:

`data_raw/`

`data_filtered/`

`qualcheck_reads/`

`assemblies/`

`qualcheck_assemblies/`

`mapping/`

`plasmids/`

`annotations/`

`resistance_genes/`

`workflow/`
 - `.snakemake`
 - `envs`
 - `scripts`

`strain_list.txt`

`config.yaml`

`snakefile`

workflow:

1. mount ARGOS
2. copy files
3. prepare files
4. get coverage
5. create config
6. run the pipeline on these files (incl. spades --plasmid)