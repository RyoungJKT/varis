"""M6: VarisDB Platform — Database, API, frontend, visualization, reports.

Built in layers. Each layer is independently deployable:
Layer 1: Database + API (FastAPI + PostgreSQL)
Layer 2: Basic frontend (React + Tailwind, or Streamlit escape hatch)
Layer 3: Charts (Plotly)
Layer 4: 3D viewer (Mol*)
Layer 5: PDF reports (WeasyPrint)
Layer 6: Search engine (Elasticsearch, fallback: PostgreSQL full-text)
Layer 7: Batch processing (Celery + Redis)

Minimum viable VarisDB = Layer 1 + Layer 2.
"""
