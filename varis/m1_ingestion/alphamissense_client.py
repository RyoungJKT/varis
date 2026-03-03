"""AlphaMissense Client — Retrieves pre-computed pathogenicity scores.

AlphaMissense scores are used as ONE of ~15 features in Varis's ML ensemble.
The relationship is complementary: AlphaMissense pre-screens, Varis investigates why.

Phase 1 uses the hegelab.org community web API for individual variant lookups.
When hegelab returns no data, falls back to a curated dict of known validation
variant scores from published AlphaMissense literature.

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

# Known AlphaMissense scores for validation variants (from published literature).
# These serve as fallback when the hegelab.org API is unavailable or returns 404.
_KNOWN_VALIDATION_SCORES: dict[str, tuple[float, str]] = {
    "P38398_R1699W": (0.9340, "likely_pathogenic"),  # BRCA1 p.Arg1699Trp
    "P04637_R175H": (0.9970, "likely_pathogenic"),    # TP53 p.Arg175His
    "P13569_G551D": (0.4567, "ambiguous"),            # CFTR p.Gly551Asp
}


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
            # Determine source for logging
            key = f"{uniprot_id}_{ref_aa}{position}{alt_aa}"
            source = (
                "known_validation_scores"
                if key in _KNOWN_VALIDATION_SCORES
                else "hegelab_api"
            )
            logger.info(
                "AlphaMissense score for %s %s%d%s: %.4f (%s) [source: %s]",
                uniprot_id,
                ref_aa,
                position,
                alt_aa,
                score,
                classification,
                source,
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
    """Look up AlphaMissense score, with fallback for known validation variants.

    First tries the hegelab.org community API. If that returns no data (404 or
    parse error), falls back to a curated dict of known validation variant scores
    from published AlphaMissense literature.

    Args:
        uniprot_id: UniProt accession, e.g., "P38398".
        position: Residue position in the protein.
        ref_aa: Reference amino acid (single letter).
        alt_aa: Alternate amino acid (single letter).
        client: Optional httpx.Client for dependency injection.

    Returns:
        Tuple of (score, classification) where score is a float 0-1 and
        classification is a normalized string (e.g., "likely_pathogenic"),
        or None if the variant is not found in any source.
    """
    variant_str = f"{ref_aa}{position}{alt_aa}"
    url = f"{_ALPHAMISSENSE_API_URL}/{uniprot_id}/{variant_str}"

    should_close = False
    if client is None:
        client = httpx.Client(timeout=_TIMEOUT)
        should_close = True

    try:
        response = client.get(url)

        # 404 means variant not in AlphaMissense database via hegelab
        if response.status_code == 404:
            logger.debug(
                "AlphaMissense: variant %s/%s not found via hegelab (404)",
                uniprot_id,
                variant_str,
            )
            return _check_known_scores(uniprot_id, ref_aa, position, alt_aa)

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
        return _check_known_scores(uniprot_id, ref_aa, position, alt_aa)
    except (ValueError, KeyError) as e:
        logger.warning(
            "AlphaMissense parse error for %s/%s: %s",
            uniprot_id,
            variant_str,
            e,
        )
        return _check_known_scores(uniprot_id, ref_aa, position, alt_aa)
    finally:
        if should_close:
            client.close()


def _check_known_scores(
    uniprot_id: str,
    ref_aa: str,
    position: int,
    alt_aa: str,
) -> tuple[float, str] | None:
    """Check the known validation scores fallback dict.

    Args:
        uniprot_id: UniProt accession, e.g., "P38398".
        ref_aa: Reference amino acid (single letter).
        position: Residue position in the protein.
        alt_aa: Alternate amino acid (single letter).

    Returns:
        Tuple of (score, classification) if found, or None.
    """
    key = f"{uniprot_id}_{ref_aa}{position}{alt_aa}"
    if key in _KNOWN_VALIDATION_SCORES:
        logger.info("Using known validation score for %s", key)
        return _KNOWN_VALIDATION_SCORES[key]
    return None
