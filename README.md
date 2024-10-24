# Machine learning detection of unstable antibiotic heteroresistance in *E. coli*

This repository contains code and certain types data for the project published in %journalname%.

The pipelines were created using [Snakemake](https://snakemake.readthedocs.io/en/stable) v7.32.4

Data analysis was performed using *R* v4.4.1 and machine learning was performed using [tidymodels](https://www.tidymodels.org/) v1.2.0.

## How to run the main analysis

First, ensure that the raw (short and long) reads are placed in `resources/data_raw/{strain}/short/` and `resources/data_raw/{strain}/long/` directories inside the project's directory.

After this is done, you can run the main pipeline (genome assembly, anotation of resistance genes, repeats and insertion sequences):

a) with this command to use conda environments:

```
# navigate to the project's directory
snakemake --snakefile workflow/snakefile.smk --configfile workflow/config.yaml --use-conda 
```

or  

b) with this command to use Apptainer containers instead:

```
# navigate to the project's directory
snakemake --snakefile workflow/snakefile.smk --configfile workflow/config.yaml --use-singularity
```

for more details on installation of snakemake and available options (number of threads, running on computer clusters etc), see the official [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/)

Then you can run the three R notebooks:

1. generate features table: `notebooks/modelling/features.qmd` (it is already available in `notebooks/modelling/data/features_strain.csv`)
2. (optional) exploratory data anlysis: `notebooks/modelling/EDA.qmd`
3. run training and validation:`notebooks/modelling/training_and_validation.Rmd`
4. comparison and analysis of the models: `notebooks/modelling/models_analysis.Rmd`

To install the same versions of R packages as were used in these notebooks, install *renv* package first and then run `renv::restore()` ([here](https://rstudio.github.io/renv/index.html) you can find *renv* documentation).

## How to run additional analyses

### Phylogenetic analysis

```
snakemake --snakefile workflow/phylogeny.smk --configfile workflow/config_phylogeny.yaml --use-conda 
```

**NB**: 27 reference strain are required.

### Analysis of the HR mutants

```
# analysis of the HR mutants
snakemake --snakefile workflow/mutants.smk --configfile workflow/config_mutants.yaml --use-conda 
```

## Raw data availability

The raw sequencing reads used in this project are available from NCBI's SRA under BioProjects PRJNA1165464 and PRJNA1083935.

## Rule graphs

1. The main analysis

![main dag](images/dag.png)

2. HR mutants analysis

![mut dag](images/dag_mutants.png)
