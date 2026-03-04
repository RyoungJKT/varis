"""Database — SQLAlchemy ORM for VarisDB.

Stores variant investigation results and job tracking.
Dialect-agnostic: SQLite for development, PostgreSQL for production.
"""

import json
import logging
import os
import uuid
from datetime import datetime, timezone
from typing import Optional

from sqlalchemy import Column, String, Float, Text, DateTime, create_engine
from sqlalchemy.orm import DeclarativeBase, Session, sessionmaker

logger = logging.getLogger(__name__)

DEFAULT_DATABASE_URL = os.getenv("DATABASE_URL", "sqlite:///data/varis.db")


# =============================================================================
# ORM Models
# =============================================================================

class Base(DeclarativeBase):
    pass


class VariantRow(Base):
    """A completed variant investigation stored in the database."""

    __tablename__ = "variants"

    variant_id = Column(String(200), primary_key=True)
    gene_symbol = Column(String(50), index=True)
    hgvs_protein = Column(String(100))
    clinvar_id = Column(String(50), nullable=True)
    classification = Column(String(50), nullable=True)
    score_ensemble = Column(Float, nullable=True)
    ensemble_version = Column(String(20), nullable=True)
    pipeline_version = Column(String(20), nullable=True)
    record_json = Column(Text, nullable=False)
    created_at = Column(DateTime, default=lambda: datetime.now(timezone.utc))
    updated_at = Column(DateTime, default=lambda: datetime.now(timezone.utc),
                        onupdate=lambda: datetime.now(timezone.utc))


class JobRow(Base):
    """A background investigation job."""

    __tablename__ = "jobs"

    job_id = Column(String(36), primary_key=True)
    variant_id = Column(String(200), nullable=True)
    status = Column(String(20), default="queued", index=True)
    current_step = Column(String(100), nullable=True)
    error_message = Column(Text, nullable=True)
    created_at = Column(DateTime, default=lambda: datetime.now(timezone.utc))
    updated_at = Column(DateTime, default=lambda: datetime.now(timezone.utc),
                        onupdate=lambda: datetime.now(timezone.utc))


# =============================================================================
# Initialization
# =============================================================================

def init_db(database_url: Optional[str] = None) -> tuple:
    """Initialize database engine and create tables.

    Args:
        database_url: SQLAlchemy database URL. Defaults to DATABASE_URL env var.

    Returns:
        Tuple of (engine, SessionLocal factory).
    """
    url = database_url or DEFAULT_DATABASE_URL
    connect_args = {}
    if url.startswith("sqlite"):
        connect_args["check_same_thread"] = False
    engine = create_engine(url, connect_args=connect_args)
    Base.metadata.create_all(engine)
    session_factory = sessionmaker(bind=engine)
    logger.info(f"Database initialized: {url}")
    return engine, session_factory


def get_session(session_factory) -> Session:
    """Get a database session from the factory.

    Args:
        session_factory: SQLAlchemy sessionmaker instance.

    Returns:
        A new database session.
    """
    return session_factory()


# =============================================================================
# Variant Operations
# =============================================================================

def save_variant_record(session: Session, variant_record) -> str:
    """Save a completed investigation to the database.

    Args:
        session: SQLAlchemy session.
        variant_record: VariantRecord instance.

    Returns:
        The variant_id of the saved record.
    """
    try:
        record_dict = variant_record.to_dict()
        record_json = json.dumps(record_dict, default=str)

        existing = session.query(VariantRow).filter_by(
            variant_id=variant_record.variant_id
        ).first()

        if existing:
            existing.gene_symbol = variant_record.gene_symbol
            existing.hgvs_protein = variant_record.hgvs_protein
            existing.clinvar_id = variant_record.clinvar_id
            existing.classification = variant_record.classification
            existing.score_ensemble = variant_record.score_ensemble
            existing.ensemble_version = variant_record.ensemble_version
            existing.pipeline_version = variant_record.pipeline_version
            existing.record_json = record_json
            existing.updated_at = datetime.now(timezone.utc)
        else:
            row = VariantRow(
                variant_id=variant_record.variant_id,
                gene_symbol=variant_record.gene_symbol,
                hgvs_protein=variant_record.hgvs_protein,
                clinvar_id=variant_record.clinvar_id,
                classification=variant_record.classification,
                score_ensemble=variant_record.score_ensemble,
                ensemble_version=variant_record.ensemble_version,
                pipeline_version=variant_record.pipeline_version,
                record_json=record_json,
            )
            session.add(row)

        session.commit()
        return variant_record.variant_id
    except Exception as e:
        session.rollback()
        logger.error(f"Failed to save variant record: {e}")
        raise


def get_variant_record(session: Session, variant_id: str) -> Optional[dict]:
    """Retrieve a variant record from the database.

    Args:
        session: SQLAlchemy session.
        variant_id: The variant identifier.

    Returns:
        Dict of the full variant record, or None if not found.
    """
    row = session.query(VariantRow).filter_by(variant_id=variant_id).first()
    if row is None:
        return None
    return json.loads(row.record_json)


def search_variants(session: Session, query: str, limit: int = 20) -> list[dict]:
    """Search variants by gene symbol, HGVS, or ClinVar ID.

    Args:
        session: SQLAlchemy session.
        query: Search query string.
        limit: Maximum number of results.

    Returns:
        List of variant summary dicts.
    """
    pattern = f"%{query}%"
    rows = (
        session.query(VariantRow)
        .filter(
            (VariantRow.gene_symbol.ilike(pattern))
            | (VariantRow.hgvs_protein.ilike(pattern))
            | (VariantRow.variant_id.ilike(pattern))
            | (VariantRow.clinvar_id.ilike(pattern))
        )
        .limit(limit)
        .all()
    )
    return [
        {
            "variant_id": r.variant_id,
            "gene": r.gene_symbol,
            "hgvs_protein": r.hgvs_protein,
            "classification": r.classification,
            "score": r.score_ensemble,
            "clinvar_id": r.clinvar_id,
            "investigation_timestamp": r.created_at.isoformat() if r.created_at else None,
        }
        for r in rows
    ]


def list_variants(
    session: Session,
    page: int = 1,
    limit: int = 20,
    sort_by: str = "created_at",
) -> tuple[list[dict], int]:
    """List all variants with pagination.

    Args:
        session: SQLAlchemy session.
        page: Page number (1-indexed).
        limit: Maximum number of results per page.
        sort_by: Column to sort by (default: created_at descending).

    Returns:
        Tuple of (list of variant summary dicts, total count).
    """
    total = session.query(VariantRow).count()

    sort_column = getattr(VariantRow, sort_by, VariantRow.created_at)
    offset = (page - 1) * limit
    rows = (
        session.query(VariantRow)
        .order_by(sort_column.desc())
        .offset(offset)
        .limit(limit)
        .all()
    )
    results = [
        {
            "variant_id": r.variant_id,
            "gene": r.gene_symbol,
            "hgvs_protein": r.hgvs_protein,
            "classification": r.classification,
            "score": r.score_ensemble,
            "clinvar_id": r.clinvar_id,
            "investigation_timestamp": r.created_at.isoformat() if r.created_at else None,
        }
        for r in rows
    ]
    return results, total


def variant_exists(session: Session, variant_id: str) -> bool:
    """Check if a variant already exists in the database.

    Args:
        session: SQLAlchemy session.
        variant_id: The variant identifier.

    Returns:
        True if the variant exists, False otherwise.
    """
    return session.query(VariantRow).filter_by(variant_id=variant_id).first() is not None


def get_variant_count(session: Session) -> int:
    """Get total number of variants in the database.

    Args:
        session: SQLAlchemy session.

    Returns:
        Total variant count.
    """
    return session.query(VariantRow).count()


def get_classification_counts(session: Session) -> dict[str, int]:
    """Get counts by classification.

    Args:
        session: SQLAlchemy session.

    Returns:
        Dict mapping classification to count.
    """
    rows = session.query(VariantRow).all()
    counts: dict[str, int] = {}
    for r in rows:
        cls = r.classification or "unknown"
        counts[cls] = counts.get(cls, 0) + 1
    return counts


def get_gene_count(session: Session) -> int:
    """Get number of distinct genes.

    Args:
        session: SQLAlchemy session.

    Returns:
        Number of distinct genes.
    """
    from sqlalchemy import func
    result = session.query(func.count(VariantRow.gene_symbol.distinct())).scalar()
    return result or 0


# =============================================================================
# Job Operations
# =============================================================================

def create_job(session: Session, variant_id: str) -> str:
    """Create a new investigation job.

    Args:
        session: SQLAlchemy session.
        variant_id: The variant being investigated.

    Returns:
        The job_id (UUID string).
    """
    job_id = str(uuid.uuid4())
    row = JobRow(
        job_id=job_id,
        variant_id=variant_id,
        status="queued",
    )
    session.add(row)
    session.commit()
    return job_id


def get_job(session: Session, job_id: str) -> Optional[dict]:
    """Retrieve job status.

    Args:
        session: SQLAlchemy session.
        job_id: The job identifier.

    Returns:
        Dict with job details, or None if not found.
    """
    row = session.query(JobRow).filter_by(job_id=job_id).first()
    if row is None:
        return None
    return {
        "job_id": row.job_id,
        "variant_id": row.variant_id,
        "status": row.status,
        "current_step": row.current_step,
        "error_message": row.error_message,
        "created_at": row.created_at.isoformat() if row.created_at else None,
        "updated_at": row.updated_at.isoformat() if row.updated_at else None,
    }


def update_job_status(session: Session, job_id: str, status: str,
                      current_step: Optional[str] = None,
                      error_message: Optional[str] = None) -> None:
    """Update job status and progress.

    Args:
        session: SQLAlchemy session.
        job_id: The job identifier.
        status: New status string.
        current_step: Current pipeline step description.
        error_message: Error message if failed.
    """
    row = session.query(JobRow).filter_by(job_id=job_id).first()
    if row:
        row.status = status
        if current_step is not None:
            row.current_step = current_step
        if error_message is not None:
            row.error_message = error_message
        row.updated_at = datetime.now(timezone.utc)
        session.commit()


def get_job_count(session: Session) -> int:
    """Get total number of jobs.

    Args:
        session: SQLAlchemy session.

    Returns:
        Total job count.
    """
    return session.query(JobRow).count()


def mark_stale_jobs_failed(session: Session) -> int:
    """On startup, mark any running jobs as failed.

    Args:
        session: SQLAlchemy session.

    Returns:
        Number of jobs marked as failed.
    """
    stale = (
        session.query(JobRow)
        .filter(JobRow.status.in_(["queued", "running"]))
        .all()
    )
    count = 0
    for job in stale:
        job.status = "failed"
        job.error_message = "Server restarted while job was in progress"
        job.updated_at = datetime.now(timezone.utc)
        count += 1
    if count:
        session.commit()
        logger.info(f"Marked {count} stale jobs as failed")
    return count
