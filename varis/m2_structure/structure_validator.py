"""Structure Validator — Validates PDB and extracts quality metrics.

Checks that the mutation residue exists in the structure, extracts pLDDT
scores from B-factor column (AlphaFold/ESMFold only), and computes a
structure quality summary.

Populates: mutation_site_present, mutation_site_plddt, plddt_mean,
           plddt_available, mutation_site_confidence_bucket, numbering_scheme,
           structure_quality_summary, preparation_steps, pdb_hash.
"""

from __future__ import annotations

import hashlib
import logging
from pathlib import Path
from typing import Optional

from Bio.PDB import PDBParser

from varis.config import PLDDT_CONFIDENCE_THRESHOLD
from varis.models.variant_record import NullReason, VariantRecord

logger = logging.getLogger(__name__)

# Sources whose B-factor column stores pLDDT scores
_PREDICTED_SOURCES = ("alphafold", "esmfold")

# pLDDT confidence buckets (thresholds from AlphaFold documentation)
_PLDDT_HIGH = 90.0
_PLDDT_MEDIUM = 70.0  # matches PLDDT_CONFIDENCE_THRESHOLD default


def validate_structure(record: VariantRecord) -> VariantRecord:
    """Validate a PDB structure and extract quality metrics.

    Parses the PDB file referenced by ``record.pdb_path``, checks that the
    mutation residue is present, extracts pLDDT scores when applicable, and
    computes a quality summary.

    Args:
        record: VariantRecord with ``pdb_path``, ``residue_position``, and
            ``structure_source`` already populated (by M1 / M2 retriever).

    Returns:
        The same VariantRecord with structure validation fields populated.
        On failure, fields are set to None with appropriate reason codes.
    """
    try:
        return _validate_structure_inner(record)
    except Exception as exc:
        logger.warning(
            "Structure validation crashed for %s: %s",
            record.variant_id or "unknown",
            exc,
        )
        _set_all_none(record, NullReason.TOOL_CRASHED)
        return record


# ---------------------------------------------------------------------------
# Internal implementation
# ---------------------------------------------------------------------------

def _validate_structure_inner(record: VariantRecord) -> VariantRecord:
    """Core validation logic — separated so the outer function can catch."""
    pdb_path = Path(record.pdb_path)
    if not pdb_path.exists():
        logger.warning("PDB file not found: %s", pdb_path)
        _set_all_none(record, NullReason.NO_DATA_AVAILABLE)
        return record

    position: Optional[int] = record.residue_position
    if position is None:
        logger.warning("No residue_position set; cannot validate mutation site")
        _set_all_none(record, NullReason.UPSTREAM_DEPENDENCY_FAILED)
        return record

    # --- Parse PDB ---
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(pdb_path))
    model = structure[0]
    chain = next(iter(model.get_chains()))

    # Standard residues only (filter heteroatoms / water)
    standard_residues = [r for r in chain.get_residues() if r.id[0] == " "]
    residue_ids = {r.id[1] for r in standard_residues}

    # --- Check mutation site ---
    if position not in residue_ids:
        logger.info(
            "Residue %d not found in structure (%d standard residues, range %d–%d)",
            position,
            len(standard_residues),
            min(residue_ids) if residue_ids else 0,
            max(residue_ids) if residue_ids else 0,
        )
        record.mutation_site_present = False
        # Record why: the residue was not found in the structure.
        # set_with_reason only records for None values, but False is meaningful
        # and downstream modules need to know the reason for the failure.
        if record.null_reasons is None:
            record.null_reasons = {}
        record.null_reasons["mutation_site_present"] = NullReason.NO_DATA_AVAILABLE
        record.coordinate_mapping_confidence = "failed"
        # Still compute hash even if site is missing
        record.pdb_hash = _compute_hash(pdb_path)
        # pLDDT fields not meaningful without mutation site
        _set_plddt_none(record, NullReason.NO_DATA_AVAILABLE)
        return record

    # Site is present
    record.mutation_site_present = True
    record.numbering_scheme = "uniprot_canonical"
    record.preparation_steps = ["validated"]

    # --- pLDDT extraction (predicted structures only) ---
    source = (record.structure_source or "").lower()
    if source in _PREDICTED_SOURCES:
        record.plddt_available = True
        _extract_plddt(record, chain, standard_residues, position)
    else:
        record.plddt_available = False
        record.set_with_reason(
            "mutation_site_plddt", None, NullReason.INTENTIONALLY_SKIPPED
        )
        record.set_with_reason(
            "plddt_mean", None, NullReason.INTENTIONALLY_SKIPPED
        )
        record.set_with_reason(
            "mutation_site_confidence_bucket", None, NullReason.INTENTIONALLY_SKIPPED
        )

    # --- Structure quality summary ---
    record.structure_quality_summary = _build_quality_summary(record, len(standard_residues))

    # --- PDB hash ---
    record.pdb_hash = _compute_hash(pdb_path)

    return record


def _extract_plddt(
    record: VariantRecord,
    chain,
    standard_residues: list,
    position: int,
) -> None:
    """Extract pLDDT scores from B-factor column of a predicted structure.

    Args:
        record: VariantRecord to populate.
        chain: BioPython Chain object.
        standard_residues: Pre-filtered list of standard residues.
        position: Residue position for the mutation site.
    """
    # Site pLDDT from CA atom B-factor
    try:
        site_residue = chain[(" ", position, " ")]
        ca_atom = site_residue["CA"]
        site_plddt = ca_atom.get_bfactor()
        record.mutation_site_plddt = round(site_plddt, 2)
    except (KeyError, Exception) as exc:
        logger.warning("Could not extract site pLDDT at position %d: %s", position, exc)
        record.set_with_reason("mutation_site_plddt", None, NullReason.TOOL_CRASHED)
        site_plddt = None

    # Mean pLDDT from all CA atoms
    ca_bfactors = []
    for residue in standard_residues:
        try:
            ca_bfactors.append(residue["CA"].get_bfactor())
        except KeyError:
            continue  # Some residues may lack CA (e.g. non-standard)

    if ca_bfactors:
        mean_plddt = sum(ca_bfactors) / len(ca_bfactors)
        record.plddt_mean = round(mean_plddt, 2)
    else:
        record.set_with_reason("plddt_mean", None, NullReason.TOOL_CRASHED)

    # Confidence bucket for mutation site
    if site_plddt is not None:
        if site_plddt >= _PLDDT_HIGH:
            record.mutation_site_confidence_bucket = "high"
        elif site_plddt >= _PLDDT_MEDIUM:
            record.mutation_site_confidence_bucket = "medium"
        else:
            record.mutation_site_confidence_bucket = "low"
    else:
        record.set_with_reason(
            "mutation_site_confidence_bucket", None, NullReason.TOOL_CRASHED
        )


def _build_quality_summary(record: VariantRecord, num_residues: int) -> str:
    """Build a JSON-serialisable quality summary dict, returned as string.

    Args:
        record: VariantRecord with pLDDT fields possibly populated.
        num_residues: Number of standard residues in the structure.

    Returns:
        JSON string with quality metrics.
    """
    import json

    summary: dict = {
        "plddt_mean": record.plddt_mean,
        "plddt_site": record.mutation_site_plddt,
        "percent_low_confidence": None,
        "num_residues": num_residues,
    }

    # Compute percent_low_confidence only when pLDDT is available
    if record.plddt_available and record.pdb_path:
        try:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("protein", str(record.pdb_path))
            model = structure[0]
            chain = next(iter(model.get_chains()))
            standard_residues = [r for r in chain.get_residues() if r.id[0] == " "]
            low_count = 0
            total_count = 0
            for residue in standard_residues:
                try:
                    bfactor = residue["CA"].get_bfactor()
                    total_count += 1
                    if bfactor < PLDDT_CONFIDENCE_THRESHOLD:
                        low_count += 1
                except KeyError:
                    continue
            if total_count > 0:
                summary["percent_low_confidence"] = round(
                    low_count / total_count * 100, 2
                )
        except Exception as exc:
            logger.warning("Could not compute percent_low_confidence: %s", exc)

    return json.dumps(summary)


def _compute_hash(file_path: Path) -> str:
    """Compute SHA-256 hash of a file.

    Args:
        file_path: Path to the file.

    Returns:
        Hex-encoded SHA-256 digest (64 characters).
    """
    sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


def _set_plddt_none(record: VariantRecord, reason: str) -> None:
    """Set all pLDDT-related fields to None with the given reason."""
    record.set_with_reason("mutation_site_plddt", None, reason)
    record.set_with_reason("plddt_mean", None, reason)
    record.set_with_reason("mutation_site_confidence_bucket", None, reason)
    record.plddt_available = False


def _set_all_none(record: VariantRecord, reason: str) -> None:
    """Set all validator output fields to None with the given reason."""
    record.set_with_reason("mutation_site_present", None, reason)
    record.set_with_reason("numbering_scheme", None, reason)
    record.set_with_reason("structure_quality_summary", None, reason)
    record.set_with_reason("preparation_steps", None, reason)
    record.set_with_reason("pdb_hash", None, reason)
    _set_plddt_none(record, reason)
