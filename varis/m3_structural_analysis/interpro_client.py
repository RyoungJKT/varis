"""InterPro Client — Functional domain identification via InterPro REST API.

Structural question: Which structural/functional domain contains the mutation?
Mutations in critical domains (kinase, BRCT, SH3) are more likely pathogenic.

This replaces the HMMER wrapper (which required a 1.5GB Pfam database download).
The InterPro API returns Pfam domain annotations for a UniProt protein directly.

Populates: domain_name, domain_id, domain_criticality, domain_start, domain_end,
           domain_available
"""

import logging
from typing import Optional

import httpx

from varis.models.variant_record import NullReason, VariantRecord

logger = logging.getLogger(__name__)

# =============================================================================
# INTERPRO API CONFIGURATION
# =============================================================================
INTERPRO_BASE_URL = "https://www.ebi.ac.uk/interpro/api"
INTERPRO_TIMEOUT = 30.0

# =============================================================================
# DOMAIN CRITICALITY CLASSIFICATION
# =============================================================================
# Curated lists of Pfam domain IDs by functional importance.
# Mutations in "critical" domains are strongly associated with pathogenicity.

_CRITICAL_DOMAINS: set[str] = {
    "PF00533",  # BRCT — BRCA1 C-terminal domain
    "PF00069",  # Pkinase — Protein kinase domain
    "PF07714",  # PK_Tyr_Ser-Thr — Protein tyrosine/serine-threonine kinase
    "PF00870",  # P53 — P53 DNA-binding domain
    "PF00076",  # RRM_1 — RNA recognition motif
    "PF00097",  # zf-C3HC4 — Zinc finger, RING-type
    "PF00096",  # zf-C2H2 — Zinc finger, C2H2 type
    "PF00400",  # WD40 — WD domain, G-beta repeat
    "PF00271",  # Helicase_C — Helicase conserved C-terminal domain
    "PF00250",  # Fork_head — Forkhead domain
    "PF00104",  # Hormone_recep — Ligand-binding domain of nuclear hormone receptor
    "PF00046",  # Homeobox — Homeobox domain
    "PF00105",  # zf-C4 — Zinc finger, C4 type
    "PF00010",  # HLH — Helix-loop-helix DNA-binding domain
    "PF00125",  # Histone — Core histone H2A/H2B/H3/H4
    "PF00626",  # Gelsolin — Gelsolin repeat
    "PF00412",  # LIM — LIM domain
}

_IMPORTANT_DOMAINS: set[str] = {
    "PF00651",  # BTB — BTB/POZ domain
    "PF00397",  # WW — WW domain
    "PF00595",  # PDZ — PDZ domain
    "PF00018",  # SH3_1 — SH3 domain
    "PF00017",  # SH2 — SH2 domain
    "PF00169",  # PH — PH domain
    "PF00307",  # CH — Calponin homology domain
    "PF00023",  # Ank — Ankyrin repeat
    "PF00560",  # LRR_1 — Leucine Rich Repeat
    "PF00855",  # PWWP — PWWP domain
    "PF01344",  # Kelch_1 — Kelch motif
    "PF00856",  # SET — SET domain
    "PF00628",  # PHD — PHD-finger
    "PF00439",  # Bromodomain — Bromodomain
}


def _classify_criticality(pfam_id: str) -> str:
    """Classify domain criticality based on Pfam ID.

    Args:
        pfam_id: Pfam accession, e.g. "PF00533".

    Returns:
        One of "critical", "important", or "peripheral".
    """
    if pfam_id in _CRITICAL_DOMAINS:
        return "critical"
    if pfam_id in _IMPORTANT_DOMAINS:
        return "important"
    return "peripheral"


def _find_domain_at_position(
    results: list[dict],
    position: int,
) -> Optional[dict]:
    """Find the Pfam domain entry that contains the given residue position.

    Searches through InterPro API results to find a domain fragment where
    start <= position <= end.

    Args:
        results: The "results" list from the InterPro API JSON response.
        position: The residue position to look up.

    Returns:
        A dict with keys: accession, name, start, end — or None if not found.
    """
    for entry in results:
        metadata = entry.get("metadata", {})
        entry_type = metadata.get("type", "")

        # Only consider domain entries (not family, homologous_superfamily, etc.)
        if entry_type != "domain":
            continue

        accession = metadata.get("accession", "")
        name = metadata.get("name", "")

        proteins = entry.get("proteins", [])
        for protein in proteins:
            locations = protein.get("entry_protein_locations", [])
            for location in locations:
                fragments = location.get("fragments", [])
                for fragment in fragments:
                    start = fragment.get("start")
                    end = fragment.get("end")
                    if start is not None and end is not None:
                        if start <= position <= end:
                            return {
                                "accession": accession,
                                "name": name,
                                "start": start,
                                "end": end,
                            }

    return None


def run_interpro(
    variant_record: VariantRecord,
    client: Optional[httpx.Client] = None,
) -> VariantRecord:
    """Identify functional domain at mutation position via InterPro REST API.

    Queries the InterPro API for Pfam domain annotations on the protein
    identified by uniprot_id, then checks whether the residue_position falls
    within any annotated domain.

    Args:
        variant_record: The shared VariantRecord (needs uniprot_id and
            residue_position).
        client: Optional httpx.Client for dependency injection (testing).
            If None, a new client is created.

    Returns:
        The variant_record with domain fields populated (or None on failure).
    """
    uniprot_id = variant_record.uniprot_id
    position = variant_record.residue_position

    # --- Guard: no UniProt ID ---
    if uniprot_id is None:
        logger.warning("InterPro skipped: no uniprot_id in variant record.")
        variant_record.set_feature_status(
            "domain", False, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.set_with_reason(
            "domain_name", None, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.set_with_reason(
            "domain_id", None, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.set_with_reason(
            "domain_start", None, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.set_with_reason(
            "domain_end", None, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.set_with_reason(
            "domain_criticality", None, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        return variant_record

    # --- Guard: no residue position ---
    if position is None:
        logger.warning("InterPro skipped: no residue_position in variant record.")
        variant_record.set_feature_status(
            "domain", False, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.set_with_reason(
            "domain_name", None, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.set_with_reason(
            "domain_id", None, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.set_with_reason(
            "domain_start", None, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.set_with_reason(
            "domain_end", None, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.set_with_reason(
            "domain_criticality", None, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        return variant_record

    # --- Query InterPro API ---
    url = f"{INTERPRO_BASE_URL}/entry/pfam/protein/uniprot/{uniprot_id}"
    own_client = False

    try:
        if client is None:
            client = httpx.Client(timeout=INTERPRO_TIMEOUT)
            own_client = True

        response = client.get(url)

        if response.status_code == 404:
            logger.info(
                "InterPro: no Pfam entries found for UniProt %s.", uniprot_id,
            )
            variant_record.set_feature_status(
                "domain", False, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "domain_name", None, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "domain_id", None, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "domain_start", None, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "domain_end", None, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "domain_criticality", None, NullReason.NO_DATA_AVAILABLE,
            )
            return variant_record

        response.raise_for_status()
        data = response.json()
        results = data.get("results", [])

        # --- Find domain containing the mutation position ---
        hit = _find_domain_at_position(results, position)

        if hit is None:
            logger.info(
                "InterPro: position %d not in any Pfam domain for %s.",
                position, uniprot_id,
            )
            variant_record.set_feature_status(
                "domain", False, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "domain_name", None, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "domain_id", None, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "domain_start", None, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "domain_end", None, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "domain_criticality", None, NullReason.NO_DATA_AVAILABLE,
            )
            return variant_record

        # --- Populate domain fields ---
        criticality = _classify_criticality(hit["accession"])

        variant_record.domain_name = hit["name"]
        variant_record.domain_id = hit["accession"]
        variant_record.domain_start = hit["start"]
        variant_record.domain_end = hit["end"]
        variant_record.domain_criticality = criticality
        variant_record.set_feature_status("domain", True)

        logger.info(
            "InterPro: position %d in %s (%s), boundaries %d-%d, "
            "criticality=%s",
            position, hit["name"], hit["accession"],
            hit["start"], hit["end"], criticality,
        )

    except httpx.TimeoutException as e:
        logger.warning(
            "InterPro timed out for %s: %s",
            variant_record.variant_id or "unknown", e,
        )
        variant_record.set_feature_status(
            "domain", False, NullReason.TIMED_OUT,
        )
        variant_record.set_with_reason(
            "domain_name", None, NullReason.TIMED_OUT,
        )
        variant_record.set_with_reason(
            "domain_id", None, NullReason.TIMED_OUT,
        )
        variant_record.set_with_reason(
            "domain_start", None, NullReason.TIMED_OUT,
        )
        variant_record.set_with_reason(
            "domain_end", None, NullReason.TIMED_OUT,
        )
        variant_record.set_with_reason(
            "domain_criticality", None, NullReason.TIMED_OUT,
        )

    except httpx.HTTPStatusError as e:
        logger.warning(
            "InterPro HTTP error for %s: %s",
            variant_record.variant_id or "unknown", e,
        )
        reason = NullReason.RATE_LIMITED if e.response.status_code == 429 else NullReason.TOOL_CRASHED
        variant_record.set_feature_status("domain", False, reason)
        variant_record.set_with_reason("domain_name", None, reason)
        variant_record.set_with_reason("domain_id", None, reason)
        variant_record.set_with_reason("domain_start", None, reason)
        variant_record.set_with_reason("domain_end", None, reason)
        variant_record.set_with_reason("domain_criticality", None, reason)

    except Exception as e:
        logger.warning(
            "InterPro failed for %s: %s",
            variant_record.variant_id or "unknown", e,
        )
        variant_record.set_feature_status(
            "domain", False, NullReason.TOOL_CRASHED,
        )
        variant_record.set_with_reason(
            "domain_name", None, NullReason.TOOL_CRASHED,
        )
        variant_record.set_with_reason(
            "domain_id", None, NullReason.TOOL_CRASHED,
        )
        variant_record.set_with_reason(
            "domain_start", None, NullReason.TOOL_CRASHED,
        )
        variant_record.set_with_reason(
            "domain_end", None, NullReason.TOOL_CRASHED,
        )
        variant_record.set_with_reason(
            "domain_criticality", None, NullReason.TOOL_CRASHED,
        )

    finally:
        if own_client and client is not None:
            client.close()

    return variant_record
