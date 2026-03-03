"""ESMFold Predictor — Fallback structure prediction via Meta's ESMFold API.

Uses the ESM Atlas API to predict protein structures for sequences up to 400
amino acids. This is the fallback path when AlphaFold DB does not have a
pre-computed structure for the protein.

Most clinically relevant proteins (e.g., BRCA1 at 1863aa) exceed ESMFold's
400aa limit, so this fallback will often be skipped with a clear reason code.
"""

import logging
from pathlib import Path

import httpx

from varis.config import ESMFOLD_API_URL, STRUCTURES_DIR
from varis.models.variant_record import NullReason, VariantRecord

logger = logging.getLogger(__name__)

# ESMFold API has a practical limit around 400 residues
ESMFOLD_MAX_SEQUENCE_LENGTH = 400


def predict_esmfold(
    variant_record: VariantRecord,
    client: httpx.Client | None = None,
) -> VariantRecord:
    """Predict structure using ESMFold API. Needs protein_sequence from M1.

    Skips prediction if a structure already exists (pdb_path is set),
    if no protein sequence is available, or if the sequence exceeds
    ESMFold's 400-residue limit.

    Args:
        variant_record: The shared VariantRecord with protein_sequence populated.
        client: Optional httpx.Client for dependency injection (testing).
                If not provided, a new client is created and closed internally.

    Returns:
        The variant_record with pdb_path, structure_source, and
        structure_source_url populated on success, or null fields with
        reason codes on failure.
    """
    # Skip if structure already exists
    if variant_record.pdb_path is not None:
        logger.info("Structure already exists at %s, skipping ESMFold", variant_record.pdb_path)
        return variant_record

    # Skip if no protein sequence available
    if variant_record.protein_sequence is None:
        logger.warning("No protein_sequence available for ESMFold prediction")
        variant_record.set_with_reason(
            "pdb_path", None, NullReason.UPSTREAM_DEPENDENCY_FAILED
        )
        return variant_record

    sequence = variant_record.protein_sequence

    # Skip if sequence is too long for ESMFold
    if len(sequence) > ESMFOLD_MAX_SEQUENCE_LENGTH:
        logger.info(
            "Sequence length %d exceeds ESMFold limit of %d, skipping",
            len(sequence),
            ESMFOLD_MAX_SEQUENCE_LENGTH,
        )
        variant_record.set_with_reason(
            "pdb_path", None, NullReason.INTENTIONALLY_SKIPPED
        )
        return variant_record

    # Determine the UniProt ID for the output filename
    uniprot_id = variant_record.uniprot_id or "unknown"

    # Create client if not provided via DI
    client_created_internally = client is None
    if client_created_internally:
        client = httpx.Client(timeout=120.0)

    try:
        logger.info(
            "Predicting structure with ESMFold for %s (%d residues)",
            uniprot_id,
            len(sequence),
        )

        response = client.post(
            ESMFOLD_API_URL,
            content=sequence,
        )

        if response.status_code == 200 and "ATOM" in response.text:
            # Ensure output directory exists
            output_dir = Path(STRUCTURES_DIR)
            output_dir.mkdir(parents=True, exist_ok=True)

            pdb_filename = f"ESMFold-{uniprot_id}.pdb"
            pdb_path = output_dir / pdb_filename
            pdb_path.write_text(response.text)

            variant_record.pdb_path = str(pdb_path)
            variant_record.structure_source = "esmfold"
            variant_record.structure_source_url = ESMFOLD_API_URL

            logger.info("ESMFold structure saved to %s", pdb_path)
        else:
            logger.warning(
                "ESMFold API returned status %d or no ATOM records",
                response.status_code,
            )
            variant_record.set_with_reason(
                "pdb_path", None, NullReason.TOOL_CRASHED
            )

    except Exception as e:
        logger.warning(
            "ESMFold prediction failed for %s: %s",
            variant_record.variant_id or "unknown",
            e,
        )
        variant_record.set_with_reason(
            "pdb_path", None, NullReason.TOOL_CRASHED
        )

    finally:
        if client_created_internally:
            client.close()

    return variant_record
