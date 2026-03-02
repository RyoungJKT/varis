"""Multiple Sequence Alignment — Aligns orthologous sequences.

Uses Clustal Omega (primary) or MAFFT (fallback) to create MSA.
Priority: 2 (easy-medium)
"""
import logging
from varis.models.variant_record import VariantRecord
logger = logging.getLogger(__name__)

def run_alignment(variant_record: VariantRecord, homologs: list[str]) -> tuple[VariantRecord, dict | None]:
    """Align homologous sequences. Returns (record, alignment_data)."""
    pass
