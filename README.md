# Sequence based detection of unstable antibiotic heteroresistance

The pipeline is created using [Snakemake](https://snakemake.readthedocs.io/en/stable) - a Python-based workflow management system for reproducible and scalable data analysis. [The "rolling" paper reference](https://f1000research.com/articles/10-33/v2)

Data analysis and modelling were performed using R and *tidyverse*.

## Main workflow

![main dag](images/dag.png)

## Data analysis / ML workflow

see 

`notebooks/modelling/training_and_validation.Rmd`, 

`notebooks/modelling/EDA.qmd`,

`notebooks/modelling/features.qmd`

## Mutants analysis workflow

![mut dag](images/dag_mutants.png)