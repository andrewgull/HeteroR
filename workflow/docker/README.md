# Containers

This directory has the Dockerfiles used to build the local containers for the analysis pipelines. It also makes use of external containers from [StaphB](https://hub.docker.com/u/staphb). These containers ensure reproducibility by locking down specific versions of software.

## Available Containers

| Image Name | Purpose | Key Software Included |
| :--- | :--- | :--- |
| **annotation** | Genome annotation and repeat detection | `trnascan-se`, `isescan`, `ismapper`, `genericrepeatfinder` |
| **assembly** | Hybrid assembly | `flye`, `medaka`, `polypolish`, `unicycler`, `bwa`, `seqkit` |
| **biopython** | Python-based sequence processing | `biopython`, `pandas`, `numpy`, `bcbio-gff` |
| **biostrings** | R-based sequence manipulation | `R`, `biostrings`, `dplyr` |
| **card_rgi** | Resistance gene identification | `rgi` |
| **default** | General purpose utilities | `bedtools`, `bcftools`, `bowtie2`, `bwa`, `fastp`, `filtlong`, `pigz`, `samtools`, `seqkit`, `spades` |
| **iqtree2** (StaphB) | Phylogenomic software | `iqtree2 v2.3.1` |
| **prokka** (StaphB) | Rapid prokaryotic genome annotation | `prokka v1.14.6` |
| **roary** (StaphB) | High speed pangenome pipeline | `roary v3.13.0` |
| **rscripts** | Data wrangling and plotting | `R`, `ggplot2`, `dplyr`, `data.table`, `ggpubr`, `readr` |

## Prerequisites

To rebuild these images, you need:
1.  **Docker** (requires sudo/root privileges).
2.  **Apptainer** (formerly Singularity) to convert the images to `.sif` format.

## How to Build the Containers

We provide a helper script to automate the build process.

**1. Run from the Project Root:**
```bash
bash workflow/docker/build_sif.sh
```

This script will use Dockerfiles located under this directory to build the containers and place them in `resources/apptainer/` directory where `snakemake` will look for them.