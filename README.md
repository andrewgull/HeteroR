# HeteroR
code for the heteroresistance project

## The project structure:

`data_raw/`

`data_filtered/`

`qualcheck_reads/`

`assemblies/`

`qualcheck_assembly/`

`mapping/`

`plasmids/`

`assemblies_joined/`

`annotations/`

`resistance_genes/`

`logs/`

`strain_lists/`

`test_dir/`

`tools/`

`localDB/`

`busco_downloads/`

`coverage/`

`final/`

`dags/`

`card_database/`

`notebooks/` - these are copies of actual notebooks that I keep on GoogleDrive. Hope these copies will be updated regularly

`workflow/`
 - `snakefile`
 - `envs`
 - `scripts`

`config.yaml`

## Workflow:

1. mount ARGOS
2. copy files
3. prepare files
4. get coverage
5. create config
6. load a local instance of CARD db (it must be in the project dir as 'localDB' - `rgi load`)
7. run the pipeline on these files: `snakemake --use-conda --cores 12 --resources mem_mb=12000`
8. run the following command to get an overview of resistance hits in your strains
   ```
   cd resistance genes; 
   for D in DA*; do ln -s "/home/andrei/Data/HeteroR/resistance_genes/"$D"/rgi_table.json" "/home/andrei/Data/HeteroR/resistance_genes/linked/"$D"_rgi_table.json"; done && 
   rgi heatmap -i linked -o heatmap -cat gene_family -clus samples
   ```

## DAG example

the most recent version

![dag](figures/dag_full.png)

## Heatmap

the most recent version

The most recent version of RGI heatmap
![resistance genes heatmap](figures/heatmap54.png)

AMR genes categorised by AMR Gene Family and samples have been clustered hierarchically (see SciPy documentation). 
Yellow represents a perfect hit, teal represents a strict hit, purple represents no hit. 
Genes with asterisks (*) appear multiple times because they belong to more than one AMR Gene Family category in the antibiotic resistance ontology (ARO).