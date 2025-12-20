
# Build arguments
ARG VERSION=2.62.0
ARG CONDA_ENV=biostrings

########################
# Build Stage
########################
FROM mambaorg/micromamba:1.5.8 AS build

# Re-declare ARG for this stage
ARG CONDA_ENV

# Ensure environment directory exists and is writable
USER root
RUN mkdir -p /opt/conda/envs && chown mambauser:mambauser /opt/conda/envs
USER mambauser

# Create the micromamba environment at a fixed path using CONDA_ENV
RUN micromamba create -y -p /opt/conda/envs/${CONDA_ENV} -c bioconda -c conda-forge \
    r-base=4.0 \
    r-dplyr=1.0.8 \
    bioconductor-biostrings=${VERSION} \
    && micromamba clean --all --yes

########################
# Runtime Stage
########################
FROM debian:bookworm-slim AS runtime

# Re-declare ARG for this stage
ARG VERSION
ARG CONDA_ENV

# Set environment variables
ENV VERSION=${VERSION}
ENV CONDA_ENV=${CONDA_ENV}
ENV PATH="/opt/conda/envs/${CONDA_ENV}/bin:${PATH}"

# Create non-root user
RUN groupadd -r appuser && useradd -r -g appuser appuser

# Copy environment from build stage (same path, same name)
COPY --from=build --chown=appuser:appuser /opt/conda/envs/${CONDA_ENV} /opt/conda/envs/${CONDA_ENV}

########################
# Metadata
########################
LABEL maintainer="andrei.guliaev@scilifelab.uu.se"
LABEL version=${VERSION}
LABEL env=${CONDA_ENV}

# Switch to non-root user
USER appuser
