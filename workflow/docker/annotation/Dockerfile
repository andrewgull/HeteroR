ARG VERSION=0.0.1
ARG ENV_NAME=annotation

# The build-stage image:
FROM continuumio/miniconda3:23.5.2-0 AS build

ARG ENV_NAME

# Install conda-pack:
RUN conda install -c conda-forge conda-pack conda-libmamba-solver && \
    conda config --set solver libmamba && \
    conda create --name ${ENV_NAME} -c bioconda -c conda-forge \
    genericrepeatfinder=1.0 \
    isescan=1.7.2.3 \
    prokka=1.14.6 \
    rgi=5.2.0 \
    trnascan-se=2.0.9 \
    && micromamba clean --all --yes

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