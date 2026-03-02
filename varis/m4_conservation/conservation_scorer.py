"""Conservation Scorer — Calculates per-position conservation from alignment.

Uses Shannon entropy and percent identity at the mutation position.
Priority: 3 (easy — custom Python, no external tool)
Populates: conservation_score, position_entropy, conserved_across_mammals
"""
import logging
from varis.models.variant_record import VariantRecord
logger = logging.getLogger(__name__)

def score_conservation(variant_record: VariantRecord, alignment: dict) -> VariantRecord:
    """Calculate conservation score at the mutation position from MSA."""
    pass

def _shannon_entropy(column: list[str]) -> float:
    """Calculate Shannon entropy for a column in the alignment."""
    pass
