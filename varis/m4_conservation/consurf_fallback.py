"""ConSurf Fallback — Pre-computed conservation scores when BLAST/alignment fails.

This is the safety net. ConSurf has pre-computed conservation for many proteins.
Always available as a fallback since it's a web API lookup.
"""
import logging
from varis.models.variant_record import VariantRecord
logger = logging.getLogger(__name__)

def fetch_consurf(variant_record: VariantRecord) -> VariantRecord:
    """Fetch pre-computed conservation score from ConSurf database."""
    pass
