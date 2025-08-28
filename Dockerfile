# === STAGE 1: Builder with Ubuntu + build deps + CMake + build Crux ===
FROM ubuntu:22.04 AS builder

# Install build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    wget \
    ca-certificates \
    curl \
    subversion \
    libcurl4-openssl-dev \
    libssl-dev \
    uuid-dev \
    zlib1g-dev \
    libpulse-dev \
    unzip \
  && rm -rf /var/lib/apt/lists/*

# Install CMake 3.26.4
RUN wget https://github.com/Kitware/CMake/releases/download/v3.26.4/cmake-3.26.4-linux-x86_64.sh \
    && chmod +x cmake-3.26.4-linux-x86_64.sh \
    && ./cmake-3.26.4-linux-x86_64.sh --skip-license --prefix=/usr/local \
    && rm cmake-3.26.4-linux-x86_64.sh

# Build Crux
COPY crux-4.3.Source.tar.gz /crux/
WORKDIR /crux
RUN tar -zxvf crux-4.3.Source.tar.gz \
  && cd crux-4.3.Source \
  && cmake . \
  && make -j$(nproc) \
  && make install

# === STAGE 2: Runtime with Micromamba + Python + Crux ===
FROM mambaorg/micromamba:1.5.7

# Runtime dependencies
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
    libssl3 \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

# Copy Crux binary
COPY --from=builder /usr/local/bin/crux /usr/local/bin/crux

# Install hypedsearch conda environment
USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml \
    && micromamba clean --all --yes \
    && rm /tmp/env.yaml

# Make micromamba Python globally available
USER root
RUN ln -sf /opt/conda/bin/python /usr/local/bin/python

# Switch back to micromamba user
USER $MAMBA_USER

# Ensure Crux is in PATH
ENV PATH="/usr/local/bin:${PATH}"
