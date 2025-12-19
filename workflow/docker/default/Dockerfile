ARG VERSION=1.0
ARG ENV_NAME=default

# The build-stage image:
FROM mambaorg/micromamba:1.5.8 AS build

ARG ENV_NAME=default

# Install tools using micromamba
# We install to /tmp/env to avoid permission issues in the build container
RUN micromamba create -y -p /tmp/env -c bioconda -c conda-forge \
    bedtools=2.31.1 \
    bwa=0.7.17 \
    fastp=0.20.1 \
    filtlong=0.2.1 \
    pigz=2.6 \
    samtools=1.13 \
    seqkit=2.0.0 \
    spades=3.15.3 \
    && micromamba clean --all --yes

# The runtime-stage image
FROM debian:bookworm-slim AS runtime

ARG VERSION
ARG ENV_NAME
ENV VERSION=${VERSION}
ENV PATH="/opt/conda/envs/${ENV_NAME}/bin:${PATH}"

# Create a non-root user
RUN groupadd -r appuser && useradd -r -g appuser appuser

# Copy the environment from the build stage
# We copy from /tmp/env to the final location
COPY --from=build --chown=appuser:appuser /tmp/env /opt/conda/envs/${ENV_NAME}

################## METADATA ######################
LABEL maintainer="andrei.guliaev@scilifelab.uu.se"
LABEL version=${VERSION}
LABEL env=${ENV_NAME}

# Switch to non-root user
USER appuser

# No default command specified, but the environment is on the PATH