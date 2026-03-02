"""ClinVar Submitter — Formats and submits evidence to ClinVar.

Generates ClinVar-format submission from Varis investigation results.
High-confidence predictions are submitted with structural evidence.
"""
import logging
from varis.models.variant_record import VariantRecord
logger = logging.getLogger(__name__)

def format_clinvar_submission(variant_record: VariantRecord) -> dict:
    """Format variant investigation into ClinVar submission format."""
    pass

def submit_to_clinvar(submission: dict) -> str | None:
    """Submit evidence to ClinVar API. Returns submission ID."""
    pass
