"""M7: Self-Evolution — Auto-retrain, tool discovery, and evolution logging.

Three autonomous loops at different independence levels:
Loop 1: Auto-Retraining (fully autonomous) — weekly ClinVar scan → retrain → deploy
Loop 2: Tool Discovery (semi-autonomous) — scout arXiv/GitHub → LLM evaluates → propose
Loop 3: Auto-Integration (partially autonomous) — install → wrapper → benchmark → integrate

Priority order:
  1. Evolution Log (must-have — records everything)
  2. Auto-retrain Loop 1 (must-have)
  3. ClinVar submission formatter (must-have)
  4. Open-source packaging (must-have)
  5. Tool scout Loop 2 (stretch goal)
  6. Auto-integration Loop 3 (stretch goal)
"""
import logging
from pathlib import Path
from typing import Optional

from varis.config import MODELS_DIR

logger = logging.getLogger(__name__)


def run(
    archive_root: Optional[Path] = None,
    log_db: Optional[object] = None,
) -> dict:
    """Run M7 auto-retrain loop — thin wrapper for pipeline integration.

    This is the entry point that the main pipeline orchestrator calls to
    trigger a retraining cycle. It delegates to auto_retrain.run_retrain_loop().

    Args:
        archive_root: Root directory for the model archive. Defaults to MODELS_DIR.
        log_db: SQLAlchemy session factory for evolution log. If None, initializes
            from DATABASE_URL config.

    Returns:
        Dict with keys: completed, decision, reason, metrics, duration.
        If an unexpected error occurs, returns a dict with completed=False.
    """
    try:
        from varis.m7_evolution.auto_retrain import run_retrain_loop

        if archive_root is None:
            archive_root = MODELS_DIR

        result = run_retrain_loop(archive_root=archive_root, log_db=log_db)
        logger.info(
            "M7 run complete: completed=%s, decision=%s",
            result.get("completed"),
            result.get("decision"),
        )
        return result
    except Exception as e:
        logger.error("M7 run failed: %s", e)
        return {
            "completed": False,
            "decision": None,
            "reason": f"M7 run failed: {e}",
            "metrics": None,
            "duration": 0.0,
        }
