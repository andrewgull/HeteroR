# Machine learning detection of unstable antibiotic heteroresistance in *E. coli*

The pipelines are created using [Snakemake](https://snakemake.readthedocs.io/en/stable)

Data analysis and modelling are performed using R and *tidyverse*.

## Snakemake pipelines

### Main pipeline

File: `workflow/snakefile.smk`

Purpose: assembling and annotating *E. coli* genomes (resistance genes, IS elements, direct repeats)

Configuration file: `workflow/config.yaml`

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

For feature generation see `notebooks/modelling/features.qmd`.

For exploratory data analysis of the features, see file `notebooks/modelling/EDA.qmd`,

For training and validation procedures, see `notebooks/modelling/training_and_validation.Rmd`,

For analysis of the models, see `notebooks/modelling/models_analysis.Rmd`

Features table: `notebooks/modelling/data/features_strain.csv`

