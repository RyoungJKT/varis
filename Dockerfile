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
    && rm -rf /var/lib/apt/lists/*

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
