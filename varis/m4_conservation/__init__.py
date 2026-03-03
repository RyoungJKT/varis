"""M4: Conservation Engine — Calculates evolutionary conservation.

CRITICAL INDEPENDENCE: M4 depends on M1 (needs protein sequence), NOT on M2 or M3.
If all structural analysis fails, conservation still works.
If conservation fails, structural analysis still works.

Depends on: M1 (needs protein_sequence and/or uniprot_id)
Populates: conservation_score, conservation_method, num_orthologs, position_entropy,
           conserved_across_mammals, msa_num_sequences, msa_gap_fraction_at_site,
           msa_column_index, conservation_available

Pipeline: cache check → UniProt orthologs → Clustal Omega → entropy scoring → ConSurf fallback
Cache: protein-level by uniprot_id — every variant in the same protein reuses the alignment.
"""

import json
import logging
from pathlib import Path
from typing import Optional

from varis.config import DATA_DIR
from varis.models.variant_record import NullReason, VariantRecord

logger = logging.getLogger(__name__)

# =============================================================================
# CONSTANTS
# =============================================================================

_CACHE_DIR = DATA_DIR / "conservation"


# =============================================================================
# PUBLIC API
# =============================================================================

def run(variant_record: VariantRecord) -> VariantRecord:
    """Execute M4: evolutionary conservation analysis. Independent of M2/M3.

    Orchestrates the full conservation pipeline:
      1. Check cache for pre-computed scores by uniprot_id
      2. Fetch orthologs from UniProt
      3. Align with Clustal Omega
      4. Score conservation (Shannon entropy)
      5. Fallback to ConSurf if primary pipeline fails
      6. Save results to cache

    Args:
        variant_record: The shared VariantRecord (needs protein_sequence
            and/or uniprot_id from M1).

    Returns:
        The variant_record with conservation fields populated, or with
        appropriate null reasons if all approaches fail.
    """
    # Guard: need protein_sequence or uniprot_id to do anything
    if variant_record.protein_sequence is None and variant_record.uniprot_id is None:
        logger.warning(
            "No protein_sequence or uniprot_id for %s — M4 cannot run",
            variant_record.variant_id,
        )
        variant_record.set_feature_status(
            "conservation", False, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.set_with_reason(
            "conservation_score", None, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.mark_module_failed("M4")
        return variant_record

    # Step 0: Check cache
    cached = _load_cache(variant_record)
    if cached is not None:
        logger.info(
            "Loaded conservation from cache for %s position %s",
            variant_record.uniprot_id,
            variant_record.residue_position,
        )
        variant_record.mark_module_completed("M4")
        return variant_record

    # Lazy imports — only load sub-modules when actually needed
    from varis.m4_conservation.uniprot_orthologs import fetch_orthologs
    from varis.m4_conservation.clustal_client import run_alignment
    from varis.m4_conservation.conservation_scorer import score_conservation
    from varis.m4_conservation.consurf_fallback import fetch_consurf

    # Step 1: Find orthologs via UniProt
    orthologs = None
    try:
        variant_record, orthologs = fetch_orthologs(variant_record)
    except Exception as e:
        logger.warning("UniProt orthologs failed: %s", e)
        variant_record.mark_module_failed("M4.orthologs")

    # Step 2: Multiple sequence alignment
    alignment = None
    if orthologs is not None:
        try:
            variant_record, alignment = run_alignment(variant_record, orthologs)
        except Exception as e:
            logger.warning("Alignment failed: %s", e)
            variant_record.mark_module_failed("M4.alignment")

    # Step 3: Score conservation from alignment
    if alignment is not None:
        try:
            variant_record = score_conservation(variant_record, alignment)
        except Exception as e:
            logger.warning("Conservation scoring failed: %s", e)
            variant_record.mark_module_failed("M4.scoring")

    # Fallback: ConSurf pre-computed scores if primary pipeline failed
    if variant_record.conservation_score is None:
        try:
            variant_record = fetch_consurf(variant_record)
        except Exception as e:
            logger.warning("ConSurf fallback failed: %s", e)
            variant_record.mark_module_failed("M4.consurf")

    # Save cache and mark module status
    if variant_record.conservation_score is not None:
        _save_cache(variant_record)
        variant_record.mark_module_completed("M4")
    else:
        variant_record.mark_module_failed("M4")

    return variant_record


# =============================================================================
# CACHE HELPERS
# =============================================================================

def _cache_path(uniprot_id: str) -> Path:
    """Return the cache file path for a given UniProt accession.

    Args:
        uniprot_id: UniProt accession (e.g., "P38398").

    Returns:
        Path to the cache JSON file.
    """
    return _CACHE_DIR / f"{uniprot_id}_scores.json"


def _load_cache(variant_record: VariantRecord) -> Optional[dict]:
    """Load cached conservation scores for the variant's position.

    Checks for a cache file at data/conservation/{uniprot_id}_scores.json.
    If the file exists and contains a score for the variant's residue
    position, populates the variant_record fields from the cache and
    returns the position data.

    Args:
        variant_record: The shared VariantRecord (needs uniprot_id and
            residue_position).

    Returns:
        The position cache entry dict if found, or None if cache miss.
    """
    uniprot_id = variant_record.uniprot_id
    position = variant_record.residue_position

    if uniprot_id is None or position is None:
        return None

    cache_file = _cache_path(uniprot_id)
    if not cache_file.exists():
        return None

    try:
        with open(cache_file, "r") as f:
            cache_data = json.load(f)

        positions = cache_data.get("positions", {})
        position_key = str(position)

        if position_key not in positions:
            return None

        pos_data = positions[position_key]

        # Populate the variant record from cache
        variant_record.conservation_score = pos_data.get("conservation_score")
        variant_record.position_entropy = pos_data.get("position_entropy")
        variant_record.msa_column_index = pos_data.get("msa_column_index")
        variant_record.msa_gap_fraction_at_site = pos_data.get("msa_gap_fraction_at_site")
        variant_record.conserved_across_mammals = pos_data.get("conserved_across_mammals")

        # Populate protein-level fields from cache root
        variant_record.conservation_method = cache_data.get("method")
        variant_record.num_orthologs = cache_data.get("num_orthologs")
        variant_record.msa_num_sequences = cache_data.get("msa_num_sequences")

        variant_record.set_feature_status("conservation", True)

        return pos_data

    except (json.JSONDecodeError, OSError, KeyError) as e:
        logger.warning(
            "Cache read failed for %s: %s", uniprot_id, e,
        )
        return None


def _save_cache(variant_record: VariantRecord) -> None:
    """Save conservation scores to cache for the variant's protein.

    Writes (or updates) the cache file at
    data/conservation/{uniprot_id}_scores.json. If the file already exists,
    the position entry is merged; protein-level fields are updated.

    Args:
        variant_record: The shared VariantRecord with conservation fields
            populated.
    """
    uniprot_id = variant_record.uniprot_id
    position = variant_record.residue_position

    if uniprot_id is None or position is None:
        return

    cache_file = _cache_path(uniprot_id)

    try:
        # Ensure cache directory exists
        _CACHE_DIR.mkdir(parents=True, exist_ok=True)

        # Load existing cache or start fresh
        cache_data: dict = {}
        if cache_file.exists():
            try:
                with open(cache_file, "r") as f:
                    cache_data = json.load(f)
            except (json.JSONDecodeError, OSError):
                cache_data = {}

        # Update protein-level fields
        cache_data["method"] = variant_record.conservation_method
        cache_data["num_orthologs"] = variant_record.num_orthologs
        cache_data["msa_num_sequences"] = variant_record.msa_num_sequences

        # Ensure positions dict exists
        if "positions" not in cache_data:
            cache_data["positions"] = {}

        # Add/update position entry
        position_key = str(position)
        cache_data["positions"][position_key] = {
            "conservation_score": variant_record.conservation_score,
            "position_entropy": variant_record.position_entropy,
            "msa_column_index": variant_record.msa_column_index,
            "msa_gap_fraction_at_site": variant_record.msa_gap_fraction_at_site,
            "conserved_across_mammals": variant_record.conserved_across_mammals,
        }

        with open(cache_file, "w") as f:
            json.dump(cache_data, f, indent=2)

        logger.info(
            "Saved conservation cache for %s position %s",
            uniprot_id, position_key,
        )

    except OSError as e:
        logger.warning(
            "Failed to save conservation cache for %s: %s",
            uniprot_id, e,
        )
