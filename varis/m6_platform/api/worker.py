"""Background worker for async variant investigations.

Uses ThreadPoolExecutor for concurrent pipeline execution.
Updates job status in the database as the pipeline progresses.
"""

import logging
from concurrent.futures import ThreadPoolExecutor

from varis.m6_platform.api.database import (
    create_job,
    update_job_status,
    save_variant_record,
    mark_stale_jobs_failed,
)

logger = logging.getLogger(__name__)

# Pipeline steps shown in job status
PIPELINE_STEPS = [
    "M1: Data Ingestion",
    "M2: Structure Engine",
    "M3: Structural Analysis",
    "M4: Conservation",
    "M5: ML Scoring",
]


def run_pipeline(gene: str, hgvs: str):
    """Run the M1-M5 investigation pipeline.

    Args:
        gene: Gene symbol.
        hgvs: HGVS protein notation.

    Returns:
        Completed VariantRecord.
    """
    from varis.pipeline import investigate
    return investigate(gene, hgvs)


class InvestigationWorker:
    """Manages background investigation jobs with a thread pool.

    Args:
        session_factory: SQLAlchemy sessionmaker for creating DB sessions.
        max_workers: Maximum concurrent investigations.
    """

    def __init__(self, session_factory, max_workers: int = 2):
        self._session_factory = session_factory
        self._executor = ThreadPoolExecutor(max_workers=max_workers)
        logger.info(f"InvestigationWorker started with {max_workers} workers")

    def submit(self, gene: str, hgvs: str, session) -> str:
        """Submit a new variant for investigation.

        Args:
            gene: Gene symbol.
            hgvs: HGVS protein notation.
            session: Active DB session for creating the job.

        Returns:
            The job_id string.
        """
        variant_id = f"{gene}_{hgvs}"
        job_id = create_job(session, variant_id)
        self._executor.submit(self._run_investigation, job_id, gene, hgvs)
        logger.info(f"Submitted job {job_id} for {variant_id}")
        return job_id

    def _run_investigation(self, job_id: str, gene: str, hgvs: str) -> None:
        """Worker thread: run pipeline and update job status.

        Args:
            job_id: The job identifier.
            gene: Gene symbol.
            hgvs: HGVS protein notation.
        """
        session = self._session_factory()
        try:
            update_job_status(session, job_id, "running", current_step=PIPELINE_STEPS[0])

            record = run_pipeline(gene, hgvs)

            save_variant_record(session, record)
            update_job_status(session, job_id, "succeeded",
                              current_step="Complete")
            logger.info(f"Job {job_id} succeeded: {record.variant_id}")
        except Exception as e:
            logger.error(f"Job {job_id} failed: {e}")
            update_job_status(session, job_id, "failed",
                              error_message=str(e))
        finally:
            session.close()

    def mark_stale_jobs(self, session) -> int:
        """Mark stale jobs as failed on startup.

        Args:
            session: Active DB session.

        Returns:
            Number of jobs marked as failed.
        """
        return mark_stale_jobs_failed(session)

    def shutdown(self) -> None:
        """Gracefully shut down the thread pool."""
        self._executor.shutdown(wait=False)
        logger.info("InvestigationWorker shut down")
