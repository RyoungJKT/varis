"""UniProt Orthologs Client — Fetch ortholog sequences from UniProt REST API.

Replaces the BLAST approach with pre-computed homologs from UniProt,
which is significantly faster than submitting a BLAST job.

Queries UniRef clusters at 50% identity to retrieve homologous sequences,
parses FASTA response to extract sequence IDs, sequences, and taxonomy
(OX= field from headers).

Priority: 1 (fast — uses pre-computed UniProt data)
Fallback: BLAST client, then ConSurf pre-computed scores.
"""

import logging
import re
from typing import Optional

import httpx

from varis.models.variant_record import NullReason, VariantRecord

logger = logging.getLogger(__name__)

# =============================================================================
# CONSTANTS
# =============================================================================

_UNIREF_STREAM_URL = (
    "https://rest.uniprot.org/uniref/stream"
    "?query=member:{uniprot_id}+AND+identity:0.5"
    "&format=fasta&size=100"
)
_TIMEOUT_SECONDS = 60
_MAX_ORTHOLOGS = 100
_MIN_SEQUENCES = 10  # Including query

# Regex to extract OX=NNNNN (taxonomy ID) from UniProt FASTA headers
_OX_PATTERN = re.compile(r"OX=(\d+)")


# =============================================================================
# PUBLIC API
# =============================================================================

def fetch_orthologs(
    variant_record: VariantRecord,
    client: Optional[httpx.Client] = None,
) -> tuple[VariantRecord, Optional[dict]]:
    """Fetch ortholog sequences from UniProt for conservation analysis.

    Queries the UniProt UniRef API for homologous sequences at 50% identity,
    parses the FASTA response, and returns a dict suitable for downstream
    alignment and conservation scoring.

    Args:
        variant_record: The shared VariantRecord (needs uniprot_id from M1).
        client: Optional httpx.Client for dependency injection in tests.
            If None, a new client is created with default timeout.

    Returns:
        A tuple of (variant_record, orthologs_dict) where orthologs_dict
        has keys:
          - "sequences": dict of seq_id -> sequence string
          - "query_id": str identifying the query protein
          - "taxonomy": dict of seq_id -> NCBI taxon ID (int)
        Returns (variant_record, None) if uniprot_id is missing, too few
        sequences found, or any error occurs.
    """
    uniprot_id = variant_record.uniprot_id
    if uniprot_id is None:
        logger.warning(
            "No uniprot_id on variant %s — skipping ortholog fetch",
            variant_record.variant_id,
        )
        return variant_record, None

    try:
        fasta_text = _fetch_fasta(uniprot_id, client)
        if fasta_text is None:
            return variant_record, None

        entries = _parse_fasta(fasta_text)
        if not entries:
            logger.warning(
                "No sequences parsed from UniProt response for %s",
                uniprot_id,
            )
            return variant_record, None

        # Cap at MAX_ORTHOLOGS
        entries = entries[:_MAX_ORTHOLOGS]

        # Build result dict
        sequences: dict[str, str] = {}
        taxonomy: dict[str, int] = {}

        for entry in entries:
            seq_id = entry["id"]
            sequences[seq_id] = entry["sequence"]
            if entry["taxon_id"] is not None:
                taxonomy[seq_id] = entry["taxon_id"]

        # Add query protein sequence if available and not already present
        query_id = f"query_{uniprot_id}"
        if variant_record.protein_sequence:
            sequences[query_id] = variant_record.protein_sequence
        else:
            # Use uniprot_id as query reference without sequence
            query_id = next(iter(sequences)) if sequences else query_id

        total_sequences = len(sequences)
        if total_sequences < _MIN_SEQUENCES:
            logger.warning(
                "Only %d sequences found for %s (minimum %d required)",
                total_sequences, uniprot_id, _MIN_SEQUENCES,
            )
            return variant_record, None

        orthologs = {
            "sequences": sequences,
            "query_id": query_id,
            "taxonomy": taxonomy,
        }

        logger.info(
            "Fetched %d orthologs for %s (%d with taxonomy)",
            len(sequences) - 1,  # exclude query
            uniprot_id,
            len(taxonomy),
        )
        return variant_record, orthologs

    except Exception as e:
        logger.warning(
            "Ortholog fetch failed for %s: %s",
            variant_record.variant_id,
            e,
        )
        return variant_record, None


# =============================================================================
# PRIVATE HELPERS
# =============================================================================

def _fetch_fasta(
    uniprot_id: str,
    client: Optional[httpx.Client] = None,
) -> Optional[str]:
    """Fetch FASTA text from UniProt UniRef API.

    Args:
        uniprot_id: UniProt accession (e.g., "P38398").
        client: Optional httpx.Client for DI. If None, creates one.

    Returns:
        Raw FASTA text, or None on failure.
    """
    url = _UNIREF_STREAM_URL.format(uniprot_id=uniprot_id)

    owns_client = client is None
    if owns_client:
        client = httpx.Client(
            timeout=_TIMEOUT_SECONDS,
            follow_redirects=True,
        )

    try:
        response = client.get(url)

        if response.status_code == 404:
            logger.warning("UniProt returned 404 for %s", uniprot_id)
            return None

        if response.status_code != 200:
            logger.warning(
                "UniProt returned HTTP %d for %s",
                response.status_code, uniprot_id,
            )
            return None

        return response.text

    except httpx.TimeoutException:
        logger.warning("UniProt request timed out for %s", uniprot_id)
        return None
    except httpx.HTTPError as e:
        logger.warning("HTTP error fetching orthologs for %s: %s", uniprot_id, e)
        return None
    finally:
        if owns_client:
            client.close()


def _parse_fasta(fasta_text: str) -> list[dict]:
    """Parse FASTA text into a list of sequence entries.

    Each entry is a dict with keys:
      - "id": sequence identifier (accession from sp|ACCESSION|NAME)
      - "sequence": amino acid sequence (all lines joined)
      - "taxon_id": NCBI taxonomy ID from OX= field, or None

    Args:
        fasta_text: Raw FASTA-formatted text from UniProt.

    Returns:
        List of parsed entries. Empty list if no valid entries found.
    """
    entries: list[dict] = []
    current_header: Optional[str] = None
    current_lines: list[str] = []

    for line in fasta_text.splitlines():
        line = line.strip()
        if not line:
            continue

        if line.startswith(">"):
            # Save previous entry
            if current_header is not None:
                entry = _header_to_entry(current_header, "".join(current_lines))
                if entry is not None:
                    entries.append(entry)

            current_header = line[1:]  # Strip leading '>'
            current_lines = []
        else:
            current_lines.append(line)

    # Save last entry
    if current_header is not None:
        entry = _header_to_entry(current_header, "".join(current_lines))
        if entry is not None:
            entries.append(entry)

    return entries


def _header_to_entry(header: str, sequence: str) -> Optional[dict]:
    """Parse a single FASTA header + sequence into an entry dict.

    Extracts the sequence ID (from sp|ACCESSION|NAME or first word)
    and the taxonomy ID (from OX=NNNNN).

    Args:
        header: The FASTA header line (without leading '>').
        sequence: The concatenated sequence lines.

    Returns:
        Dict with "id", "sequence", "taxon_id", or None if sequence is empty.
    """
    if not sequence:
        return None

    # Extract sequence ID: try sp|ACCESSION|NAME format first
    parts = header.split("|")
    if len(parts) >= 2:
        seq_id = parts[1]
    else:
        # Fallback: first whitespace-delimited token
        seq_id = header.split()[0] if header.split() else header

    # Extract taxonomy ID from OX=NNNNN
    taxon_id: Optional[int] = None
    ox_match = _OX_PATTERN.search(header)
    if ox_match:
        taxon_id = int(ox_match.group(1))

    return {
        "id": seq_id,
        "sequence": sequence,
        "taxon_id": taxon_id,
    }
