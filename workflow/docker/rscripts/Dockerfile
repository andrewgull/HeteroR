ARG VERSION=4.4.3
ARG CONDA_ENV=rscripts

# The build-stage image:
FROM mambaorg/micromamba:1.5.8 AS build

ARG CONDA_ENV=rscripts

# Install tools using micromamba
# R requires the environment to be at the same path as it was built.
# We use a consistent path /opt/conda/envs/${CONDA_ENV} for both build and runtime.
# We ensure the parent directory is writable so micromamba can create the environment.
USER root
RUN mkdir -p /opt/conda/envs && \
    chown mambauser:mambauser /opt/conda/envs
USER mambauser

RUN micromamba create -y -p /opt/conda/envs/${CONDA_ENV} -c bioconda -c conda-forge \
    r-base=4.3 \
    r-ggplot2=3.4.4 \
    r-stringr=1.5.1 \
    r-dplyr=1.1.4 \
    r-data.table=1.14.8 \
    r-ggpubr=0.6.0 \
    r-readr=2.1.5 \
    r-r.utils=2.12.3 \
    r-purrr=1.0.2 \
    bioconductor-treeio=1.26.0 \
    bioconductor-ggtree=3.10.0 \
    && micromamba clean --all --yes

# The runtime-stage image
FROM debian:bookworm-slim AS runtime

ARG VERSION
ARG CONDA_ENV
ENV VERSION=${VERSION}
ENV PATH="/opt/conda/envs/${CONDA_ENV}/bin:${PATH}"

# Create a non-root user
RUN groupadd -r appuser && useradd -r -g appuser appuser

# Copy the environment from the build stage
# We copy from the SAME path to the SAME path to ensure R works
COPY --from=build --chown=appuser:appuser /opt/conda/envs/${CONDA_ENV} /opt/conda/envs/${CONDA_ENV}

################## METADATA ######################
LABEL maintainer="andrei.guliaev@scilifelab.uu.se"
LABEL version=${VERSION}
LABEL env=${CONDA_ENV}

# Switch to non-root user
USER appuser
