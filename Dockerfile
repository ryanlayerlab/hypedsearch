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
FROM ghcr.io/astral-sh/uv:python3.11-trixie

# Copy Crux binary from starge 1 to this stage
COPY --from=builder /usr/local/bin/crux /usr/local/bin/crux

WORKDIR /app

# Copy the project into the image
COPY pyproject.toml ./

# Disable development dependencies
ENV UV_NO_DEV=1

# Sync the project into a new environment, asserting the lockfile is up to date
WORKDIR /app
RUN uv sync