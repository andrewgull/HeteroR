# Machine learning detection of unstable antibiotic heteroresistance in *E. coli*

This repository contains code and certain types data for the project published in %journalname%.

The pipelines are created using [Snakemake](https://snakemake.readthedocs.io/en/stable) v7.32.4

Data analysis was performed using *R* v4.4.1 and machine learning was performed using *tidymodels* v1.2.0.

To reproduce the main analysis (genomes assembly and anotation), you have to run the main pipeline

with this command to use conda environments:

```
snakemake --snakefile workflow/snakefile.smk --configfile workflow/config.yaml --use-conda 
```

or with this command to use Apptainer containers instead:

```
snakemake --snakefile workflow/snakefile.smk --configfile workflow/config.yaml --use-singularity
```

for more details on posiible options and installation of snakemake (number of threads, running on computer clusters, installing snakemake environment etc, see the official [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/))

## Snakemake pipelines

### Main pipeline

File: `workflow/snakefile.smk`

Purpose: assembling and annotating *E. coli* genomes (resistance genes, IS elements, direct repeats) from both short and long sequencing reads.

Configuration file: `workflow/config.yaml`

To run the pipeline short and long reads should be in `resources/data_raw/{strain}/short/` and `resources/data_raw/{strain}/long/` directories.

DAG:

![main dag](images/dag.png)

### Phylogeny pipeline

File: `workflow/phylogeny.smk`

Purpose: phylogenetic analysis of the samples including 27 reference strains.

Configuration file: `workflow/config_phylogeny.yaml`

### Analysis of mutants

File: `mutants.smk`

Purpose: analysis of the HR mutants.

Configuration file: `workflow/config_mutants.yaml`

DAG: 

![mut dag](images/dag_mutants.png)

## Data analysis and machine learning

R package *renv* was used to achieve reproducibility of the analysis in the `notebooks`

To get the same versions of packages, use `renv::restore()` and then you can run code in the notebooks:

 - For feature generation - `notebooks/modelling/features.qmd`.

 - For exploratory data analysis of the features - `notebooks/modelling/EDA.qmd`,

 - For training and validation procedures - `notebooks/modelling/training_and_validation.Rmd`,

 - For analysis of the models - `notebooks/modelling/models_analysis.Rmd`

The final features table used for training/validating and testing is here: `notebooks/modelling/data/features_strain.csv`

## Raw reads availability

The raw sequencing reads used in this project will be available from NCBI's SRA under BioProject PRJNA1165464
