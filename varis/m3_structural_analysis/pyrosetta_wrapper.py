"""PyRosetta Wrapper — Fallback ΔΔG stability prediction (requires license).

This is a stub. Full implementation will be added when a PyRosetta academic
license is obtained.

Populates: ddg_pyrosetta (when available).
"""
import logging
from varis.models.variant_record import NullReason, VariantRecord

logger = logging.getLogger(__name__)


def run_pyrosetta(variant_record: VariantRecord) -> VariantRecord:
    """Run PyRosetta ΔΔG prediction (stub — requires license).

    Checks whether the pyrosetta package is importable. If not,
    records tool_missing and returns. If found, records not_attempted
    (full implementation pending license).

    Args:
        variant_record: The shared VariantRecord with pdb_path set.

    Returns:
        The variant_record with ddg_pyrosetta updated.
    """
    try:
        import pyrosetta  # noqa: F401
    except ImportError:
        logger.info("PyRosetta not installed — requires academic license")
        variant_record.set_with_reason("ddg_pyrosetta", None, NullReason.TOOL_MISSING)
        return variant_record

    # TODO: Full PyRosetta implementation when license is obtained
    logger.info("PyRosetta found but full implementation pending")
    variant_record.set_with_reason("ddg_pyrosetta", None, NullReason.NOT_ATTEMPTED)
    return variant_record
