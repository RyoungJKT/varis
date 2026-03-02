"""UniProt Client — Retrieves protein sequence and functional annotations.

Queries the UniProt REST API to resolve a gene symbol to a reviewed (Swiss-Prot)
accession, then fetches the full protein entry for sequence, name, and function.

Populates: uniprot_id, protein_name, protein_sequence, protein_length, protein_function.
"""

from __future__ import annotations

import logging
from typing import Optional

import httpx

from varis.models.variant_record import VariantRecord, NullReason
from varis.config import UNIPROT_BASE_URL

logger = logging.getLogger(__name__)

_TIMEOUT = 15.0


def fetch_uniprot(
    variant_record: VariantRecord,
    client: Optional[httpx.Client] = None,
) -> VariantRecord:
    """Query UniProt for protein sequence and metadata.

    Searches UniProt for the gene symbol, preferring reviewed (Swiss-Prot) entries,
    then fetches the full protein entry to extract name, sequence, and function.

    Args:
        variant_record: Must have gene_symbol set.
        client: Optional httpx.Client for dependency injection (e.g., testing).
            If None, a new client is created internally.

    Returns:
        VariantRecord with protein fields populated, or None fields with
        reasons on failure.
    """
    gene = variant_record.gene_symbol
    if not gene:
        logger.warning("No gene_symbol set on variant record; skipping UniProt lookup")
        for field_name in ("uniprot_id", "protein_name", "protein_sequence",
                           "protein_length", "protein_function"):
            variant_record.set_with_reason(
                field_name, None, NullReason.UPSTREAM_DEPENDENCY_FAILED
            )
        return variant_record

    try:
        accession = _search_uniprot(gene, client=client)
        if accession is None:
            logger.warning("UniProt search returned no results for gene %s", gene)
            for field_name in ("uniprot_id", "protein_name", "protein_sequence",
                               "protein_length", "protein_function"):
                variant_record.set_with_reason(
                    field_name, None, NullReason.NO_DATA_AVAILABLE
                )
            return variant_record

        variant_record.set_with_reason("uniprot_id", accession)

        protein_data = _fetch_protein_data(accession, client=client)
        if protein_data is None:
            logger.warning(
                "UniProt fetch returned no data for accession %s", accession
            )
            for field_name in ("protein_name", "protein_sequence",
                               "protein_length", "protein_function"):
                variant_record.set_with_reason(
                    field_name, None, NullReason.NO_DATA_AVAILABLE
                )
            return variant_record

        # Extract protein name
        protein_name = _extract_protein_name(protein_data)
        variant_record.set_with_reason(
            "protein_name",
            protein_name,
            NullReason.NO_DATA_AVAILABLE if protein_name is None else None,
        )

        # Extract sequence
        sequence = _extract_sequence(protein_data)
        variant_record.set_with_reason(
            "protein_sequence",
            sequence,
            NullReason.NO_DATA_AVAILABLE if sequence is None else None,
        )

        # Derive protein length from sequence
        protein_length = len(sequence) if sequence else None
        variant_record.set_with_reason(
            "protein_length",
            protein_length,
            NullReason.NO_DATA_AVAILABLE if protein_length is None else None,
        )

        # Extract function
        protein_function = _extract_function(protein_data)
        variant_record.set_with_reason(
            "protein_function",
            protein_function,
            NullReason.NO_DATA_AVAILABLE if protein_function is None else None,
        )

        logger.info(
            "UniProt data retrieved for %s: accession=%s, name=%s, length=%s",
            gene,
            accession,
            protein_name,
            protein_length,
        )

    except Exception as e:
        logger.warning(
            "UniProt lookup failed for %s: %s", gene, e
        )
        for field_name in ("uniprot_id", "protein_name", "protein_sequence",
                           "protein_length", "protein_function"):
            variant_record.set_with_reason(
                field_name, None, NullReason.TOOL_CRASHED
            )

    return variant_record


def _search_uniprot(
    gene: str,
    client: Optional[httpx.Client] = None,
) -> Optional[str]:
    """Search UniProt for a human gene and return the best accession ID.

    Prefers reviewed (Swiss-Prot) entries over unreviewed (TrEMBL).

    Args:
        gene: Gene symbol, e.g., "BRCA1".
        client: Optional httpx.Client for dependency injection.

    Returns:
        UniProt accession ID (e.g., "P38398"), or None if not found.
    """
    url = f"{UNIPROT_BASE_URL}/uniprotkb/search"
    params = {
        "query": f"gene_exact:{gene} AND organism_id:9606",
        "format": "json",
        "fields": "accession,reviewed,protein_name,gene_names",
        "size": "5",
    }

    should_close = False
    if client is None:
        client = httpx.Client(timeout=_TIMEOUT)
        should_close = True

    try:
        response = client.get(url, params=params)
        response.raise_for_status()
        data = response.json()
    finally:
        if should_close:
            client.close()

    results = data.get("results", [])
    if not results:
        return None

    # Prefer reviewed (Swiss-Prot) entries
    reviewed = [
        r for r in results
        if r.get("entryType") == "UniProtKB reviewed (Swiss-Prot)"
    ]

    if reviewed:
        return reviewed[0].get("primaryAccession")

    # Fall back to the first result (unreviewed TrEMBL)
    return results[0].get("primaryAccession")


def _fetch_protein_data(
    accession: str,
    client: Optional[httpx.Client] = None,
) -> Optional[dict]:
    """Fetch full protein entry from UniProt REST API.

    Args:
        accession: UniProt accession ID, e.g., "P38398".
        client: Optional httpx.Client for dependency injection.

    Returns:
        Parsed protein data dict, or None on failure.
    """
    url = f"{UNIPROT_BASE_URL}/uniprotkb/{accession}.json"

    should_close = False
    if client is None:
        client = httpx.Client(timeout=_TIMEOUT)
        should_close = True

    try:
        response = client.get(url)
        response.raise_for_status()
        return response.json()
    finally:
        if should_close:
            client.close()


def _extract_protein_name(protein_data: dict) -> Optional[str]:
    """Extract the recommended protein name from UniProt entry data.

    Falls back to submittedName if recommendedName is absent.

    Args:
        protein_data: Full UniProt protein entry as parsed JSON dict.

    Returns:
        Protein name string, or None if not found.
    """
    description = protein_data.get("proteinDescription", {})

    # Primary: recommendedName
    rec_name = description.get("recommendedName", {})
    full_name = rec_name.get("fullName", {})
    name = full_name.get("value")
    if name:
        return name

    # Fallback: submittedName (common for unreviewed entries)
    submitted_names = description.get("submissionNames", [])
    if submitted_names:
        submitted_full = submitted_names[0].get("fullName", {})
        return submitted_full.get("value")

    return None


def _extract_sequence(protein_data: dict) -> Optional[str]:
    """Extract the protein amino acid sequence from UniProt entry data.

    Args:
        protein_data: Full UniProt protein entry as parsed JSON dict.

    Returns:
        Amino acid sequence string, or None if not found.
    """
    sequence_block = protein_data.get("sequence", {})
    return sequence_block.get("value")


def _extract_function(protein_data: dict) -> Optional[str]:
    """Extract functional description from UniProt entry comments.

    Looks for comments with commentType == "FUNCTION" and concatenates
    all function text values.

    Args:
        protein_data: Full UniProt protein entry as parsed JSON dict.

    Returns:
        Function description string, or None if no function annotation found.
    """
    comments = protein_data.get("comments", [])
    function_texts: list[str] = []

    for comment in comments:
        if comment.get("commentType") == "FUNCTION":
            texts = comment.get("texts", [])
            for text_entry in texts:
                value = text_entry.get("value")
                if value:
                    function_texts.append(value)

    if function_texts:
        return " ".join(function_texts)

    return None
