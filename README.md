# HeteroR
code for heteroresistance project

project structure:

`data_raw/`

`data_filtered/`

`qualcheck_reads/`

`assemblies/`

`qualcheck_assemblies/`

`workflow/`

`strain_list.txt`

`config.yaml`

workflow:

1. mount ARGOS
2. copy files
3. prepare files
4. get coverage
5. (create config)
6. run the pipeline on these files (+ back-mapping, spades)