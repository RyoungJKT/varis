"""BLAST Client — Finds homologous protein sequences across species.

Uses NCBI BLAST web API (no local installation required).
Priority: 1 (build first — easy via web API)
Fallback: PSI-BLAST, or skip to ConSurf pre-computed scores.
"""
import logging
from varis.models.variant_record import VariantRecord
logger = logging.getLogger(__name__)

def run_blast(variant_record: VariantRecord) -> tuple[VariantRecord, list[str] | None]:
    """Run BLAST to find orthologs. Returns (record, list_of_homolog_sequences)."""
    pass
