"""ConSurf Fallback — Pre-computed conservation scores when BLAST/alignment fails.

This is the safety net. ConSurf has pre-computed conservation for many proteins.
Always available as a fallback since it's a web API lookup.

ConSurf grades range from 1 (variable) to 9 (conserved). We map to a normalized
score: score = (grade - 1) / 8.0, giving 0.0 (variable) to 1.0 (conserved).

Priority: 3 (fallback — only runs when primary pipeline and BLAST both fail)
Fallback for: UniProt orthologs + Clustal Omega + conservation scorer
"""

import logging
from typing import Optional

import httpx

from varis.models.variant_record import NullReason, VariantRecord

logger = logging.getLogger(__name__)

# =============================================================================
# CONSTANTS
# =============================================================================

_CONSURF_URL = (
    "https://consurf.tau.ac.il/results/{uniprot_id}/consurf_grades.json"
)
_TIMEOUT_SECONDS = 30
_MIN_GRADE = 1
_MAX_GRADE = 9


# =============================================================================
# PUBLIC API
# =============================================================================

def fetch_consurf(
    variant_record: VariantRecord,
    client: Optional[httpx.Client] = None,
) -> VariantRecord:
    """Fetch pre-computed conservation score from ConSurf database.

    Looks up the ConSurf pre-computed conservation grade for the variant's
    residue position, then maps it to a normalized 0.0-1.0 score.

    Args:
        variant_record: The shared VariantRecord (needs uniprot_id and
            residue_position from M1).
        client: Optional httpx.Client for dependency injection in tests.
            If None, a new client is created with default timeout.

    Returns:
        The variant_record with conservation fields populated on success,
        or with conservation_available=False and appropriate null reasons
        on failure.
    """
    # Guard: need uniprot_id to query ConSurf
    uniprot_id = variant_record.uniprot_id
    if uniprot_id is None:
        logger.warning(
            "No uniprot_id on variant %s — skipping ConSurf lookup",
            variant_record.variant_id,
        )
        variant_record.set_feature_status(
            "conservation", False, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.set_with_reason(
            "conservation_score", None, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        return variant_record

    try:
        grades = _fetch_grades(uniprot_id, client)
        if grades is None:
            variant_record.set_feature_status(
                "conservation", False, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "conservation_score", None, NullReason.NO_DATA_AVAILABLE,
            )
            return variant_record

        # Look up the grade at the variant's residue position
        position = variant_record.residue_position
        position_key = str(position) if position is not None else None

        if position_key is None or position_key not in grades:
            logger.warning(
                "ConSurf has no grade for position %s of %s",
                position_key, uniprot_id,
            )
            variant_record.set_feature_status(
                "conservation", False, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "conservation_score", None, NullReason.NO_DATA_AVAILABLE,
            )
            return variant_record

        grade_entry = grades[position_key]
        grade = grade_entry.get("grade")

        if grade is None or not (_MIN_GRADE <= grade <= _MAX_GRADE):
            logger.warning(
                "Invalid ConSurf grade %s at position %s of %s",
                grade, position_key, uniprot_id,
            )
            variant_record.set_feature_status(
                "conservation", False, NullReason.VALIDATION_FAILED,
            )
            variant_record.set_with_reason(
                "conservation_score", None, NullReason.VALIDATION_FAILED,
            )
            return variant_record

        # Map grade 1-9 to score 0.0-1.0
        score = (grade - _MIN_GRADE) / (_MAX_GRADE - _MIN_GRADE)

        variant_record.conservation_score = score
        variant_record.conservation_method = "consurf"
        variant_record.set_feature_status("conservation", True)

        logger.info(
            "ConSurf fallback: %s position %s grade=%d score=%.3f",
            uniprot_id, position_key, grade, score,
        )
        return variant_record

    except httpx.TimeoutException:
        logger.warning(
            "ConSurf request timed out for %s", uniprot_id,
        )
        variant_record.set_feature_status(
            "conservation", False, NullReason.TIMED_OUT,
        )
        variant_record.set_with_reason(
            "conservation_score", None, NullReason.TIMED_OUT,
        )
        return variant_record

    except Exception as e:
        logger.warning(
            "ConSurf fallback failed for %s: %s",
            variant_record.variant_id, e,
        )
        variant_record.set_feature_status(
            "conservation", False, NullReason.TOOL_CRASHED,
        )
        variant_record.set_with_reason(
            "conservation_score", None, NullReason.TOOL_CRASHED,
        )
        return variant_record


# =============================================================================
# PRIVATE HELPERS
# =============================================================================

def _fetch_grades(
    uniprot_id: str,
    client: Optional[httpx.Client] = None,
) -> Optional[dict]:
    """Fetch ConSurf grades JSON for a given UniProt accession.

    Args:
        uniprot_id: UniProt accession (e.g., "P38398").
        client: Optional httpx.Client for DI. If None, creates one.

    Returns:
        Dict of position -> grade entry, or None on failure.
    """
    url = _CONSURF_URL.format(uniprot_id=uniprot_id)

    owns_client = client is None
    if owns_client:
        client = httpx.Client(
            timeout=_TIMEOUT_SECONDS,
            follow_redirects=True,
        )

    try:
        response = client.get(url)

        if response.status_code == 404:
            logger.warning("ConSurf returned 404 for %s", uniprot_id)
            return None

        if response.status_code != 200:
            logger.warning(
                "ConSurf returned HTTP %d for %s",
                response.status_code, uniprot_id,
            )
            return None

        data = response.json()
        grades = data.get("grades")
        if grades is None:
            logger.warning(
                "ConSurf response missing 'grades' key for %s", uniprot_id,
            )
            return None

        return grades

    except httpx.TimeoutException:
        raise  # Let the caller handle timeout specifically

    except httpx.HTTPError as e:
        logger.warning(
            "HTTP error fetching ConSurf for %s: %s", uniprot_id, e,
        )
        return None

    finally:
        if owns_client:
            client.close()
