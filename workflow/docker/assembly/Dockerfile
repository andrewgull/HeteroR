ARG VERSION=0.0.1
ARG ENV_NAME=hybrid_assembly

# The build-stage image:
FROM continuumio/miniconda3:23.5.2-0 AS build

ARG ENV_NAME

# Install conda-pack:
RUN conda install -c conda-forge conda-pack conda-libmamba-solver && \
    conda config --set solver libmamba && \
    conda create --name ${ENV_NAME} -c bioconda -c conda-forge \
    flye=2.9.1 \
    medaka=1.6.0 \
    polypolish=0.5.0 \
    unicycler=0.4.8 \
    bwa=0.7.17 \
    seqkit=2.2.0

# The runtime-stage image; we can use Debian as the
# base image since the Conda env also includes Python
# for us.
FROM debian:buster-slim AS runtime

ARG VERSION
ARG ENV_NAME
ENV VERSION=${VERSION}

################## METADATA ######################
LABEL maintainer="andrei.guliaev@scilifelab.uu.se"
LABEL version=${VERSION}
LABEL env=${ENV_NAME}
   
COPY --from=build /opt/conda/envs/${ENV_NAME}/ /opt/conda/envs/${ENV_NAME}/
ENV PATH=${PATH}:/opt/conda/envs/${ENV_NAME}/bin