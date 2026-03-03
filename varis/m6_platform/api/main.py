"""FastAPI Application — VarisDB REST API.

Endpoints:
  GET  /health                          - Health check
  GET  /api/v1/investigations/{id}      - Full investigation payload
  POST /api/v1/investigations           - Submit new investigation
  GET  /api/v1/variants?q={query}       - Search variants
  GET  /api/v1/variants/stats           - Database statistics
  GET  /api/v1/jobs/{job_id}            - Job status
"""

import logging
import os
import time
from collections import defaultdict
from contextlib import asynccontextmanager
from typing import Optional

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

from varis.m6_platform.api.models import (
    InvestigationResponse,
    InvestigationRequest,
    SearchResponse,
    StatsResponse,
    JobStatusResponse,
    VariantSummary,
)
from varis.m6_platform.api.database import (
    init_db,
    get_variant_record,
    search_variants,
    get_variant_count,
    get_classification_counts,
    get_gene_count,
    get_job,
    get_job_count,
    variant_exists,
    mark_stale_jobs_failed,
)
from varis.m6_platform.api.validation import validate_variant_input
from varis.m6_platform.api.investigation import build_investigation_response
from varis.m6_platform.api.worker import InvestigationWorker
from varis.models.variant_record import VariantRecord

logger = logging.getLogger(__name__)

# Module-level state set during app creation
_session_factory = None
_worker = None


def create_app(database_url: Optional[str] = None) -> FastAPI:
    """Create and configure the FastAPI application.

    Args:
        database_url: SQLAlchemy database URL. Defaults to env var or SQLite.

    Returns:
        Configured FastAPI application.
    """
    global _session_factory, _worker

    db_url = database_url or os.getenv("DATABASE_URL", "sqlite:///data/varis.db")
    cors_origins = os.getenv("CORS_ORIGINS", "http://localhost:5173").split(",")
    max_workers = int(os.getenv("MAX_INVESTIGATION_WORKERS", "2"))
    rate_limit = int(os.getenv("RATE_LIMIT_PER_MINUTE", "10"))

    engine, session_factory = init_db(db_url)
    _session_factory = session_factory

    # Initialize worker eagerly so it's available for TestClient
    session = session_factory()
    mark_stale_jobs_failed(session)
    session.close()
    _worker = InvestigationWorker(session_factory, max_workers=max_workers)

    @asynccontextmanager
    async def lifespan(app: FastAPI):
        logger.info("VarisDB API started")
        yield
        if _worker:
            _worker.shutdown()
        logger.info("VarisDB API stopped")

    app = FastAPI(
        title="VarisDB API",
        version="1.0.0",
        description="Variant structural investigation database",
        lifespan=lifespan,
    )

    app.add_middleware(
        CORSMiddleware,
        allow_origins=cors_origins,
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    # Rate limiter state
    _rate_counts: dict[str, list[float]] = defaultdict(list)

    def _check_rate_limit(client_ip: str) -> bool:
        """Check if client has exceeded rate limit."""
        now = time.time()
        window = now - 60
        _rate_counts[client_ip] = [t for t in _rate_counts[client_ip] if t > window]
        if len(_rate_counts[client_ip]) >= rate_limit:
            return False
        _rate_counts[client_ip].append(now)
        return True

    # === Health Check ===

    @app.get("/health")
    def health_check():
        return {"status": "ok", "service": "varisdb"}

    # === Investigation Endpoints ===

    @app.get("/api/v1/investigations/{variant_id}", response_model=InvestigationResponse)
    def get_investigation(variant_id: str):
        session = session_factory()
        try:
            record_dict = get_variant_record(session, variant_id)
            if record_dict is None:
                raise HTTPException(status_code=404, detail=f"Variant '{variant_id}' not found")
            record = VariantRecord.from_dict(record_dict)
            return build_investigation_response(record)
        finally:
            session.close()

    @app.post("/api/v1/investigations", status_code=202)
    def submit_investigation(request: InvestigationRequest, req: Request):
        client_ip = req.client.host if req.client else "unknown"
        if not _check_rate_limit(client_ip):
            raise HTTPException(status_code=429, detail="Rate limit exceeded. Try again later.")

        validation = validate_variant_input(request.gene, request.hgvs)
        if not validation["valid"]:
            raise HTTPException(status_code=422, detail=validation["error"])

        session = session_factory()
        try:
            vid = f"{request.gene}_{request.hgvs}"
            if variant_exists(session, vid):
                return {"variant_id": vid, "status": "exists"}

            job_id = _worker.submit(request.gene, request.hgvs, session)
            return {"job_id": job_id, "status": "queued"}
        finally:
            session.close()

    # === Variant Search ===

    @app.get("/api/v1/variants", response_model=SearchResponse)
    def search(q: str = "", page: int = 1, limit: int = 20):
        session = session_factory()
        try:
            if not q:
                return SearchResponse(results=[], total=0, page=page, limit=limit)
            results = search_variants(session, q, limit=limit)
            summaries = [
                VariantSummary(
                    variant_id=r["variant_id"],
                    gene=r["gene"],
                    hgvs_protein=r.get("hgvs_protein"),
                    classification=r.get("classification"),
                    score=r.get("score"),
                    clinvar_id=r.get("clinvar_id"),
                    investigation_timestamp=r.get("investigation_timestamp"),
                )
                for r in results
            ]
            return SearchResponse(results=summaries, total=len(summaries), page=page, limit=limit)
        finally:
            session.close()

    # === Stats ===

    @app.get("/api/v1/variants/stats", response_model=StatsResponse)
    def get_stats():
        session = session_factory()
        try:
            total = get_variant_count(session)
            counts = get_classification_counts(session)
            genes = get_gene_count(session)
            jobs = get_job_count(session)
            return StatsResponse(
                total_variants=total,
                total_pathogenic=counts.get("likely_pathogenic", 0),
                total_benign=counts.get("likely_benign", 0),
                total_uncertain=counts.get("uncertain", 0),
                total_jobs=jobs,
                genes_covered=genes,
            )
        finally:
            session.close()

    # === Job Status ===

    @app.get("/api/v1/jobs/{job_id}", response_model=JobStatusResponse)
    def get_job_status(job_id: str):
        session = session_factory()
        try:
            job = get_job(session, job_id)
            if job is None:
                raise HTTPException(status_code=404, detail=f"Job '{job_id}' not found")
            return JobStatusResponse(**job)
        finally:
            session.close()

    return app


# For uvicorn: uvicorn varis.m6_platform.api.main:app --reload --port 8000
app = create_app()
