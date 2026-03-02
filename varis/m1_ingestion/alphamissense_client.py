"""AlphaMissense Client — Retrieves pre-computed pathogenicity scores.

AlphaMissense scores are used as ONE of ~15 features in Varis's ML ensemble.
The relationship is complementary: AlphaMissense pre-screens, Varis investigates why.

Populates: alphamissense_score, alphamissense_class.
"""

import logging
from varis.models.variant_record import VariantRecord

logger = logging.getLogger(__name__)


def fetch_alphamissense(variant_record: VariantRecord) -> VariantRecord:
    """Look up AlphaMissense pre-computed pathogenicity score.

    Args:
        variant_record: Must have gene_symbol and hgvs_protein set.

    Returns:
        VariantRecord with AlphaMissense score and classification.
    """
    pass


def _lookup_score(gene: str, position: int, ref_aa: str, alt_aa: str) -> tuple[float, str] | None:
    """Look up AlphaMissense score from local database or API.

    Args:
        gene: Gene symbol.
        position: Residue position.
        ref_aa: Reference amino acid (single letter).
        alt_aa: Alternate amino acid (single letter).

    Returns:
        Tuple of (score, classification) or None if not found.
    """
    pass
