"""AlphaMissense Client — Retrieves pre-computed pathogenicity scores.

AlphaMissense scores are used as ONE of ~15 features in Varis's ML ensemble.
The relationship is complementary: AlphaMissense pre-screens, Varis investigates why.

Phase 1 uses the hegelab.org community web API for individual variant lookups.

Populates: alphamissense_score, alphamissense_class.
"""

from __future__ import annotations

import logging
from typing import Optional

import httpx

from varis.models.variant_record import VariantRecord, NullReason

logger = logging.getLogger(__name__)

_TIMEOUT = 15.0
_ALPHAMISSENSE_API_URL = "https://alphamissense.hegelab.org/api/variant"


def fetch_alphamissense(
    variant_record: VariantRecord,
    client: Optional[httpx.Client] = None,
) -> VariantRecord:
    """Look up AlphaMissense pre-computed pathogenicity score.

    Requires uniprot_id, residue_position, ref_aa_single, and alt_aa_single
    to be set on the variant record. If any are missing, both output fields
    are set to None with UPSTREAM_DEPENDENCY_FAILED.

    Args:
        variant_record: Must have uniprot_id, residue_position, ref_aa_single,
            and alt_aa_single set.
        client: Optional httpx.Client for dependency injection (e.g., testing).
            If None, a new client is created internally.

    Returns:
        VariantRecord with alphamissense_score and alphamissense_class populated,
        or None fields with reasons on failure.
    """
    uniprot_id = variant_record.uniprot_id
    position = variant_record.residue_position
    ref_aa = variant_record.ref_aa_single
    alt_aa = variant_record.alt_aa_single

    # Check that all required upstream fields are available
    if not all([uniprot_id, position, ref_aa, alt_aa]):
        missing = []
        if not uniprot_id:
            missing.append("uniprot_id")
        if not position:
            missing.append("residue_position")
        if not ref_aa:
            missing.append("ref_aa_single")
        if not alt_aa:
            missing.append("alt_aa_single")
        logger.warning(
            "AlphaMissense lookup skipped: missing upstream fields %s",
            ", ".join(missing),
        )
        for field_name in ("alphamissense_score", "alphamissense_class"):
            variant_record.set_with_reason(
                field_name, None, NullReason.UPSTREAM_DEPENDENCY_FAILED
            )
        return variant_record

    try:
        result = _lookup_score(uniprot_id, position, ref_aa, alt_aa, client=client)

        if result is not None:
            score, classification = result
            variant_record.set_with_reason("alphamissense_score", score)
            variant_record.set_with_reason("alphamissense_class", classification)
            logger.info(
                "AlphaMissense score for %s %s%d%s: %.4f (%s)",
                uniprot_id,
                ref_aa,
                position,
                alt_aa,
                score,
                classification,
            )
        else:
            logger.warning(
                "AlphaMissense returned no data for %s %s%d%s",
                uniprot_id,
                ref_aa,
                position,
                alt_aa,
            )
            for field_name in ("alphamissense_score", "alphamissense_class"):
                variant_record.set_with_reason(
                    field_name, None, NullReason.NO_DATA_AVAILABLE
                )

    except Exception as e:
        logger.warning(
            "AlphaMissense lookup failed for %s %s%d%s: %s",
            uniprot_id,
            ref_aa,
            position,
            alt_aa,
            e,
        )
        for field_name in ("alphamissense_score", "alphamissense_class"):
            variant_record.set_with_reason(
                field_name, None, NullReason.TOOL_CRASHED
            )

    return variant_record


def _lookup_score(
    uniprot_id: str,
    position: int,
    ref_aa: str,
    alt_aa: str,
    client: Optional[httpx.Client] = None,
) -> tuple[float, str] | None:
    """Look up AlphaMissense score from the hegelab.org community API.

    Constructs a variant string (e.g., "R1699W") and queries the API
    endpoint for that specific UniProt ID and variant.

    Args:
        uniprot_id: UniProt accession, e.g., "P38398".
        position: Residue position in the protein.
        ref_aa: Reference amino acid (single letter).
        alt_aa: Alternate amino acid (single letter).
        client: Optional httpx.Client for dependency injection.

    Returns:
        Tuple of (score, classification) where score is a float 0-1 and
        classification is a normalized string (e.g., "likely_pathogenic"),
        or None if the variant is not found or the request fails.
    """
    variant_str = f"{ref_aa}{position}{alt_aa}"
    url = f"{_ALPHAMISSENSE_API_URL}/{uniprot_id}/{variant_str}"

    should_close = False
    if client is None:
        client = httpx.Client(timeout=_TIMEOUT)
        should_close = True

    try:
        response = client.get(url)

        # 404 means variant not in AlphaMissense database
        if response.status_code == 404:
            logger.debug(
                "AlphaMissense: variant %s/%s not found (404)",
                uniprot_id,
                variant_str,
            )
            return None

        response.raise_for_status()
        data = response.json()

        score = float(data["am_pathogenicity"])
        classification = data["am_class"]

        # Normalize classification: lowercase, replace spaces with underscores
        classification = classification.lower().replace(" ", "_")

        return (score, classification)

    except httpx.HTTPStatusError as e:
        logger.warning(
            "AlphaMissense HTTP error for %s/%s: %s",
            uniprot_id,
            variant_str,
            e,
        )
        return None
    except (ValueError, KeyError) as e:
        logger.warning(
            "AlphaMissense parse error for %s/%s: %s",
            uniprot_id,
            variant_str,
            e,
        )
        return None
    finally:
        if should_close:
            client.close()
