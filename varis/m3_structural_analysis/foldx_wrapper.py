"""FoldX Wrapper — ΔΔG stability prediction (requires license).

This is a stub that checks for the FoldX binary. If not found, sets
ddg_available=False with reason tool_missing. Full implementation will be
added when a FoldX academic license is obtained.

Populates: ddg_foldx (when available).
Sets: ddg_available, ddg_missing_reason.
"""
import logging
import shutil
from varis.models.variant_record import NullReason, VariantRecord

logger = logging.getLogger(__name__)


def run_foldx(variant_record: VariantRecord) -> VariantRecord:
    """Run FoldX ΔΔG prediction (stub — requires license).

    Checks for the FoldX binary in PATH. If not found, records
    tool_missing and returns. If found, records not_attempted
    (full implementation pending license).

    Args:
        variant_record: The shared VariantRecord with pdb_path set.

    Returns:
        The variant_record with ddg_foldx and ddg feature status updated.
    """
    if not shutil.which("foldx") and not shutil.which("FoldX"):
        logger.info("FoldX binary not found — requires academic license")
        variant_record.set_with_reason("ddg_foldx", None, NullReason.TOOL_MISSING)
        variant_record.set_feature_status("ddg", False, NullReason.TOOL_MISSING)
        return variant_record

    # TODO: Full FoldX implementation when license is obtained
    logger.info("FoldX binary found but full implementation pending")
    variant_record.set_with_reason("ddg_foldx", None, NullReason.NOT_ATTEMPTED)
    variant_record.set_feature_status("ddg", False, NullReason.NOT_ATTEMPTED)
    return variant_record
