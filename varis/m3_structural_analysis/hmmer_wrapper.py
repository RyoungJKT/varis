"""HMMER — Functional domain identification via Pfam.

Structural question: Which structural/functional domain contains the mutation?
Mutations in critical domains (kinase, BRCT, SH3) are more likely pathogenic.

Priority: 8 (medium — large database download)
Fallback: InterProScan API, then UniProt domain annotations
Populates: domain_name, domain_id, domain_criticality
"""
import logging
from varis.models.variant_record import VariantRecord
logger = logging.getLogger(__name__)

def run_hmmer(variant_record: VariantRecord) -> VariantRecord:
    """Identify functional domain at mutation position using HMMER + Pfam."""
    pass
