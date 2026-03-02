"""UniProt Client — Retrieves protein sequence and functional annotations.

Populates: uniprot_id, protein_name, protein_sequence, protein_length, protein_function.
"""

import logging
from varis.models.variant_record import VariantRecord
from varis.config import UNIPROT_BASE_URL

logger = logging.getLogger(__name__)


def fetch_uniprot(variant_record: VariantRecord) -> VariantRecord:
    """Query UniProt for protein sequence and metadata.

    Args:
        variant_record: Must have gene_symbol set.

    Returns:
        VariantRecord with protein fields populated.
    """
    pass


def _search_uniprot(gene: str) -> str | None:
    """Search UniProt for a human gene and return the accession ID.

    Args:
        gene: Gene symbol, e.g., "BRCA1".

    Returns:
        UniProt accession ID (e.g., "P38398"), or None if not found.
    """
    pass


def _fetch_protein_data(accession: str) -> dict | None:
    """Fetch protein entry from UniProt REST API.

    Args:
        accession: UniProt accession ID.

    Returns:
        Parsed protein data dict, or None on failure.
    """
    pass
