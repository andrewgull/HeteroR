# Machine learning detection of unstable antibiotic heteroresistance in *E. coli*

This repository contains code and certain types data for the project published in %journalname%.

The pipelines were created using [Snakemake](https://snakemake.readthedocs.io/en/stable) v8.23.1

Data analysis was performed using *R* v4.4.1 and machine learning was performed using [tidymodels](https://www.tidymodels.org/) v1.2.0.

This project contains 3 pipelines:

- (1) to run hybrid assembly of genomes and annotation of resistance genes, direct repeats and IS elements.

- (2) to run analysis of HR mutants

- (3) to run core-genome based phylogeny

## How to run the assembly-annotation pipeline (1)

1. create a directory for your project (in all following steps, we will assume that you are inside of that directory).

2. navigate to this directory and download the repository using: `git clone https://github.com/andrewgull/HeteroR`

3. create `resources/raw/` and place the raw data (short and long reads) there. Naming convention: long reads have `.fastq.gz` extension and short reads have `.fq.gz`extension.

4. install [conda/mamba](https://github.com/conda-forge/miniforge#mambaforge) and [snakemake](https://snakemake.readthedocs.io/en/stable)

5. activate snakemake environment: `conda activate snakemake`

6. run the pipeline using this command:

```bash
# use both conda and Apptainer
# substitute $N with a number you want to use
snakemake --snakefile workflow/assembly-annotation.smk --use-conda --use-singularity --cores $N
```

or

```bash
# use only conda
# substitute $N with a number you want to use
snakemake --snakefile workflow/assembly-annotation.smk --use-conda --cores $N
```

After the main pipeline has finished, you can run the three R notebooks (but not necessarily all of them):

1. to generate features table: `notebooks/modelling/features.qmd` (it is already available in `notebooks/modelling/data/features_strain.csv`)
2. for exploratory data anlysis: `notebooks/modelling/EDA.qmd`
3. to run training and validation:`notebooks/modelling/training_and_validation.Rmd`
4. for comparison and analysis of the models: `notebooks/modelling/models_analysis.Rmd`

To ensure that you use the same versions of R packages as were used in these notebooks, install *renv* package and run `renv::restore()` ([here](https://rstudio.github.io/renv/index.html) you can find *renv* documentation).

## How to run the additional analyses

### Analysis of the HR mutants (2)

```bash
# analysis of the HR mutants
snakemake --snakefile workflow/mutants.smk --use-conda --use-singularity --cores $N
```

### Phylogenetic analysis (3)

```bash
snakemake --snakefile workflow/phylogeny.smk --use-conda --use-singularity --cores $N
```

**NB**: 31 reference strain are required (see the publication for reference numbers).

## Configuration & Settings

Lists of strain names used in each of the pipelines can be found in `configs/strains_*.txt` files.

Settings of each software tool used can be found in `configs/config_*.yaml` files.

## Raw data availability

The raw sequencing reads used in this project are available from NCBI's SRA under BioProjects PRJNA1165464 and PRJNA1083935.

## Models and features table

The pre-compiled features table is available here: `notebooks/modelling/data/features_strain.csv`

The final models (trained LLR and GBT) are available here: `notebooks/modelling/models/llr_final.rds` & `notebooks/modelling/models/gbt_final.rds`

## Rule graphs

1. The main analysis
![main dag](images/dag.png)
2. HR mutants analysis
![mut dag](images/dag_mutants.png)
