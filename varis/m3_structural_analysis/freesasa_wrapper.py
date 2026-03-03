"""FreeSASA — Solvent Accessible Surface Area calculation.

Structural question: How buried is the mutated residue?
Buried residues (low SASA) in the hydrophobic core are structurally critical.
Mutations there are devastating. Surface residues are more tolerant.

Priority: 1 (build first — easy to install, always works)
Populates: solvent_accessibility_relative, burial_category, sasa_available
"""

import logging
from typing import Optional

from varis.models.variant_record import NullReason, VariantRecord

logger = logging.getLogger(__name__)

# Max SASA values per residue type in extended conformation (Tien et al. 2013).
# Used to normalize absolute SASA to relative SASA (0.0 – 1.0).
_MAX_SASA: dict[str, float] = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}

# Burial threshold: relative SASA below this is "core", above is "surface".
_BURIAL_THRESHOLD: float = 0.25


def _get_pdb_path(record: VariantRecord) -> Optional[str]:
    """Return the best available PDB path (prefer fixed structure).

    Args:
        record: The VariantRecord with structure paths.

    Returns:
        Path string to the PDB file, or None if no structure available.
    """
    return record.pdb_fixed_path or record.pdb_path


def _find_residue_area(
    residue_areas: dict,
    chain_id: str,
    residue_number: int,
) -> Optional[object]:
    """Find the ResidueArea for a specific chain and residue number.

    FreeSASA's residueAreas() returns a nested dict:
      { chain_id_str: { residue_number_str: ResidueArea } }

    The key format may vary by version, so we try several approaches.

    Args:
        residue_areas: The dict returned by freesasa Result.residueAreas().
        chain_id: Chain identifier, e.g. "A".
        residue_number: Residue number to find, e.g. 1699.

    Returns:
        The matching ResidueArea object, or None if not found.
    """
    # Approach 1: nested dict keyed by chain, then residue number as string
    chain_data = residue_areas.get(chain_id)
    if chain_data is not None and isinstance(chain_data, dict):
        # Try string key
        res_key = str(residue_number)
        if res_key in chain_data:
            return chain_data[res_key]
        # Try int key
        if residue_number in chain_data:
            return chain_data[residue_number]
        # Iterate and match by residueNumber attribute
        for _key, area in chain_data.items():
            if hasattr(area, "residueNumber") and area.residueNumber == residue_number:
                return area

    # Approach 2: flat dict with composite keys (older versions)
    # Keys might be "A,1699" or "A, 1699"
    for fmt in [f"{chain_id},{residue_number}", f"{chain_id}, {residue_number}"]:
        if fmt in residue_areas:
            return residue_areas[fmt]

    # Approach 3: iterate all entries looking for matching residue number
    for key, value in residue_areas.items():
        if isinstance(value, dict):
            # Nested dict format — already tried above for the right chain
            continue
        if hasattr(value, "residueNumber") and value.residueNumber == residue_number:
            return value

    return None


def run_freesasa(variant_record: VariantRecord) -> VariantRecord:
    """Calculate solvent accessibility for the mutated residue.

    Computes absolute SASA using FreeSASA, normalizes by residue-type maximum
    (Tien et al. 2013), and classifies burial: core (<0.25) / surface (>=0.25).

    Args:
        variant_record: The shared VariantRecord (needs pdb_path and residue_position).

    Returns:
        The variant_record with SASA fields populated (or None on failure).
    """
    pdb_path = _get_pdb_path(variant_record)
    if pdb_path is None:
        logger.warning("FreeSASA skipped: no PDB structure available.")
        variant_record.set_feature_status(
            "sasa", False, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.set_with_reason(
            "solvent_accessibility_relative", None,
            NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.set_with_reason(
            "burial_category", None, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        return variant_record

    residue_pos = variant_record.residue_position
    if residue_pos is None:
        logger.warning("FreeSASA skipped: no residue_position in variant record.")
        variant_record.set_feature_status(
            "sasa", False, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.set_with_reason(
            "solvent_accessibility_relative", None,
            NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.set_with_reason(
            "burial_category", None, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        return variant_record

    try:
        import freesasa
    except ImportError:
        logger.warning("FreeSASA skipped: freesasa package not installed.")
        variant_record.set_feature_status("sasa", False, NullReason.TOOL_MISSING)
        variant_record.set_with_reason(
            "solvent_accessibility_relative", None, NullReason.TOOL_MISSING,
        )
        variant_record.set_with_reason(
            "burial_category", None, NullReason.TOOL_MISSING,
        )
        return variant_record

    try:
        structure = freesasa.Structure(pdb_path)
        result = freesasa.calc(structure)
        residue_areas = result.residueAreas()

        # Find the residue area for the mutation site.
        # AlphaFold structures typically use chain "A".
        chain_id = "A"
        area = _find_residue_area(residue_areas, chain_id, residue_pos)

        if area is None:
            logger.warning(
                "FreeSASA: residue %d not found in chain %s of %s.",
                residue_pos, chain_id, pdb_path,
            )
            variant_record.set_feature_status(
                "sasa", False, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "solvent_accessibility_relative", None,
                NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "burial_category", None, NullReason.NO_DATA_AVAILABLE,
            )
            return variant_record

        # Get absolute SASA for this residue
        absolute_sasa: float = area.total

        # Get residue type for normalization
        residue_type: str = area.residueType  # e.g. "ARG"
        max_sasa = _MAX_SASA.get(residue_type)

        if max_sasa is None or max_sasa <= 0:
            logger.warning(
                "FreeSASA: unknown residue type %r at position %d, "
                "cannot normalize SASA.",
                residue_type, residue_pos,
            )
            variant_record.set_feature_status(
                "sasa", False, NullReason.VALIDATION_FAILED,
            )
            variant_record.set_with_reason(
                "solvent_accessibility_relative", None,
                NullReason.VALIDATION_FAILED,
            )
            variant_record.set_with_reason(
                "burial_category", None, NullReason.VALIDATION_FAILED,
            )
            return variant_record

        # Compute relative SASA, clamped to [0.0, 1.0]
        relative_sasa: float = min(max(absolute_sasa / max_sasa, 0.0), 1.0)

        # Classify burial
        burial: str = "core" if relative_sasa < _BURIAL_THRESHOLD else "surface"

        # Write results
        variant_record.solvent_accessibility_relative = relative_sasa
        variant_record.burial_category = burial
        variant_record.set_feature_status("sasa", True)

        logger.info(
            "FreeSASA: residue %d (%s) — absolute=%.2f, relative=%.3f, "
            "category=%s",
            residue_pos, residue_type, absolute_sasa, relative_sasa, burial,
        )

    except Exception as e:
        logger.warning(
            "FreeSASA failed for %s: %s",
            variant_record.variant_id or "unknown", e,
        )
        variant_record.set_feature_status("sasa", False, NullReason.TOOL_CRASHED)
        variant_record.set_with_reason(
            "solvent_accessibility_relative", None, NullReason.TOOL_CRASHED,
        )
        variant_record.set_with_reason(
            "burial_category", None, NullReason.TOOL_CRASHED,
        )

    return variant_record
