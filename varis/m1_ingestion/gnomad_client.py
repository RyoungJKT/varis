"""gnomAD Client — Retrieves population allele frequencies from gnomAD (Broad Institute).

Populates: gnomad_frequency, gnomad_popmax, gnomad_homozygotes.
"""

import logging
from varis.models.variant_record import VariantRecord
from varis.config import GNOMAD_API_URL

logger = logging.getLogger(__name__)


def fetch_gnomad(variant_record: VariantRecord) -> VariantRecord:
    """Query gnomAD for allele frequency across populations.

    Args:
        variant_record: Must have gene_symbol and residue_position set.

    Returns:
        VariantRecord with gnomAD frequency fields populated.
    """
    pass


def _query_gnomad_graphql(gene: str, position: int, ref_aa: str, alt_aa: str) -> dict | None:
    """Execute GraphQL query against gnomAD API.

    Args:
        gene: Gene symbol.
        position: Residue position.
        ref_aa: Reference amino acid (single letter).
        alt_aa: Alternate amino acid (single letter).

    Returns:
        gnomAD response dict, or None on failure.
    """
    pass
