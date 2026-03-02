"""AlphaFold DB Client — Downloads predicted 3D protein structures.

Downloads PDB files from the AlphaFold Protein Structure Database.
Queries the AlphaFold API to get the correct PDB URL (model version may change).
Caches downloaded files to avoid re-downloading for the same protein.

This is the primary structure source. ESMFold (in M2) is the fallback.

Populates: structure_source, pdb_path (stored in data/structures/).
"""

import logging
from pathlib import Path
import httpx
from varis.models.variant_record import VariantRecord, NullReason
from varis.config import STRUCTURES_DIR

logger = logging.getLogger(__name__)

_TIMEOUT = 30.0
_ALPHAFOLD_API_URL = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"


def fetch_alphafold_structure(variant_record: VariantRecord,
                              client: httpx.Client | None = None) -> VariantRecord:
    """Download AlphaFold predicted structure for the protein.

    Args:
        variant_record: Must have uniprot_id set (from uniprot_client).
        client: Optional httpx.Client for dependency injection.

    Returns:
        VariantRecord with pdb_path set to downloaded file, or None if unavailable.
    """
    uniprot_id = variant_record.uniprot_id
    if not uniprot_id:
        logger.info("No UniProt ID available — skipping AlphaFold download")
        for f in ("structure_source", "pdb_path"):
            variant_record.set_with_reason(f, None, NullReason.UPSTREAM_DEPENDENCY_FAILED)
        return variant_record

    try:
        pdb_path = _download_pdb(uniprot_id, STRUCTURES_DIR, client)
        if pdb_path:
            variant_record.structure_source = "alphafold"
            variant_record.pdb_path = str(pdb_path)
            logger.info(f"AlphaFold structure: {pdb_path}")
        else:
            logger.warning(f"AlphaFold structure not available for {uniprot_id}")
            for f in ("structure_source", "pdb_path"):
                variant_record.set_with_reason(f, None, NullReason.NO_DATA_AVAILABLE)

    except Exception as e:
        logger.warning(f"AlphaFold download failed for {uniprot_id}: {e}")
        for f in ("structure_source", "pdb_path"):
            variant_record.set_with_reason(f, None, NullReason.TOOL_CRASHED)

    return variant_record


def _download_pdb(uniprot_id: str, output_dir: Path,
                  client: httpx.Client | None = None) -> Path | None:
    """Download PDB file from AlphaFold DB, with caching.

    Queries the AlphaFold API to discover the correct PDB URL (model version
    may change over time), then downloads and caches the file.

    Args:
        uniprot_id: UniProt accession ID.
        output_dir: Directory to save the PDB file.
        client: Optional httpx.Client.

    Returns:
        Path to downloaded PDB file, or None on failure.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Check for any cached AlphaFold PDB for this protein
    cached = list(output_dir.glob(f"AF-{uniprot_id}-F1-model_v*.pdb"))
    if cached and cached[0].stat().st_size > 0:
        logger.info(f"Using cached structure: {cached[0]}")
        return cached[0]

    should_close = False
    if client is None:
        client = httpx.Client(timeout=_TIMEOUT, follow_redirects=True)
        should_close = True

    try:
        # Query the AlphaFold API to get the correct PDB URL
        api_url = _ALPHAFOLD_API_URL.format(uniprot_id=uniprot_id)
        api_resp = client.get(api_url)
        if api_resp.status_code == 404:
            logger.info(f"No AlphaFold entry for {uniprot_id}")
            return None
        api_resp.raise_for_status()

        data = api_resp.json()
        if not data:
            return None

        entry = data[0] if isinstance(data, list) else data
        pdb_url = entry.get("pdbUrl")
        if not pdb_url:
            logger.warning(f"AlphaFold API returned no pdbUrl for {uniprot_id}")
            return None

        # Extract filename from URL
        filename = pdb_url.rsplit("/", 1)[-1]
        output_path = output_dir / filename

        # Download the PDB file
        resp = client.get(pdb_url)
        if resp.status_code == 404:
            return None
        resp.raise_for_status()

        output_path.write_bytes(resp.content)
        logger.info(f"Downloaded AlphaFold structure: {output_path} ({len(resp.content)} bytes)")
        return output_path
    finally:
        if should_close:
            client.close()
