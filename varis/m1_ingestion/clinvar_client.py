"""ClinVar Client — Retrieves variant classification and clinical data from ClinVar (NIH).

Uses NCBI E-utilities API to search and fetch variant records.
Populates: clinvar_id, clinvar_classification, clinvar_review_status, clinvar_conditions.
"""

import logging
from varis.models.variant_record import VariantRecord
from varis.config import CLINVAR_BASE_URL, NCBI_API_KEY

logger = logging.getLogger(__name__)


def fetch_clinvar(variant_record: VariantRecord) -> VariantRecord:
    """Query ClinVar for variant classification and clinical significance.

    Args:
        variant_record: Must have gene_symbol and hgvs_protein set.

    Returns:
        VariantRecord with ClinVar fields populated, or None fields on failure.
    """
    pass


def _search_clinvar(gene: str, variant: str) -> str | None:
    """Search ClinVar for a variant and return the ClinVar variation ID.

    Args:
        gene: Gene symbol, e.g., "BRCA1".
        variant: HGVS protein notation, e.g., "p.Arg1699Trp".

    Returns:
        ClinVar variation ID string, or None if not found.
    """
    pass


def _fetch_clinvar_record(variation_id: str) -> dict | None:
    """Fetch full ClinVar record for a given variation ID.

    Args:
        variation_id: ClinVar variation ID.

    Returns:
        Parsed ClinVar record as dict, or None on failure.
    """
    pass
