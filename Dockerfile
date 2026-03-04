# Varis — Structural Evidence for Every Variant
# Multi-stage build for reproducible deployment

FROM python:3.11-slim AS base

WORKDIR /app

# System dependencies for structural biology tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    dssp \
    blast2 \
    hmmer \
    clustalo \
    wget \
    g++ \
    git \
    && rm -rf /var/lib/apt/lists/*

# Optional: compile EvoEF2 from source (MIT licensed)
ARG INSTALL_EVOEF2=true
RUN if [ "$INSTALL_EVOEF2" = "true" ]; then \
    git clone --depth 1 https://github.com/tommyhuangthu/EvoEF2.git /tmp/EvoEF2 && \
    cd /tmp/EvoEF2 && g++ -O3 -ffast-math -o EvoEF2 src/*.cpp && \
    cp EvoEF2 /usr/local/bin/ && \
    cp -r library /usr/local/bin/library && \
    rm -rf /tmp/EvoEF2; \
    fi

# Python dependencies
COPY pyproject.toml .
RUN pip install --no-cache-dir ".[all]"

# Application code
COPY varis/ varis/
COPY tests/ tests/

# Data directories
RUN mkdir -p data/structures data/models data/logs data/validation

EXPOSE 8000

# Default: run the API
CMD ["uvicorn", "varis.m6_platform.api.main:create_app", "--host", "0.0.0.0", "--port", "8000"]
