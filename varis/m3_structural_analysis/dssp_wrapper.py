"""DSSP Wrapper — Assigns secondary structure at mutation site.

Uses BioPython's DSSP interface with the mkdssp binary to assign
secondary structure codes (H=helix, E=sheet, C/other=coil).

Populates: secondary_structure, secondary_structure_name.
Sets: dssp_available, dssp_missing_reason.
"""

import logging
import shutil
from typing import Optional

from varis.models.variant_record import NullReason, VariantRecord

logger = logging.getLogger(__name__)

# Map DSSP 8-state codes to 3-state names
_DSSP_TO_NAME: dict[str, str] = {
    "H": "helix", "G": "helix", "I": "helix",  # All helix types
    "E": "sheet", "B": "sheet",                   # Sheet types
    "T": "coil", "S": "coil", "-": "coil", "C": "coil",  # Coil/other
}


def _get_pdb_path(record: VariantRecord) -> Optional[str]:
    """Return the best available PDB path (prefer fixed structure).

    Args:
        record: The VariantRecord with structure paths.

    Returns:
        Path string to the PDB file, or None if no structure available.
    """
    return record.pdb_fixed_path or record.pdb_path


def run_dssp(variant_record: VariantRecord) -> VariantRecord:
    """Assign secondary structure at mutation site using DSSP.

    Uses BioPython's DSSP interface with the external mkdssp binary to
    compute 8-state secondary structure, then maps to 3-state (helix/sheet/coil).

    Args:
        variant_record: The shared VariantRecord (needs pdb_path and residue_position).

    Returns:
        The variant_record with DSSP fields populated (or None on failure).
    """
    pdb_path = _get_pdb_path(variant_record)
    if pdb_path is None:
        logger.warning("DSSP skipped: no PDB structure available.")
        variant_record.set_feature_status(
            "dssp", False, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.set_with_reason(
            "secondary_structure", None,
            NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.set_with_reason(
            "secondary_structure_name", None,
            NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        return variant_record

    residue_pos = variant_record.residue_position
    if residue_pos is None:
        logger.warning("DSSP skipped: no residue_position in variant record.")
        variant_record.set_feature_status(
            "dssp", False, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.set_with_reason(
            "secondary_structure", None,
            NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.set_with_reason(
            "secondary_structure_name", None,
            NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        return variant_record

    # Check mkdssp binary is installed
    if shutil.which("mkdssp") is None:
        logger.warning("DSSP skipped: mkdssp binary not found on PATH.")
        variant_record.set_feature_status(
            "dssp", False, NullReason.TOOL_MISSING,
        )
        variant_record.set_with_reason(
            "secondary_structure", None, NullReason.TOOL_MISSING,
        )
        variant_record.set_with_reason(
            "secondary_structure_name", None, NullReason.TOOL_MISSING,
        )
        return variant_record

    try:
        from Bio.PDB import PDBParser
        from Bio.PDB.DSSP import DSSP
    except ImportError:
        logger.warning("DSSP skipped: BioPython not installed.")
        variant_record.set_feature_status(
            "dssp", False, NullReason.TOOL_MISSING,
        )
        variant_record.set_with_reason(
            "secondary_structure", None, NullReason.TOOL_MISSING,
        )
        variant_record.set_with_reason(
            "secondary_structure_name", None, NullReason.TOOL_MISSING,
        )
        return variant_record

    try:
        # Parse the PDB structure
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("variant", pdb_path)
        model = structure[0]

        # Run DSSP on the first model
        dssp = DSSP(model, pdb_path, dssp="mkdssp")

        # Get the first chain ID from the model
        chains = list(model.get_chains())
        if not chains:
            logger.warning("DSSP: no chains found in structure %s.", pdb_path)
            variant_record.set_feature_status(
                "dssp", False, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "secondary_structure", None,
                NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "secondary_structure_name", None,
                NullReason.NO_DATA_AVAILABLE,
            )
            return variant_record

        chain_id = chains[0].id

        # DSSP keys are (chain_id, (' ', residue_number, ' '))
        dssp_key = (chain_id, (' ', residue_pos, ' '))

        if dssp_key not in dssp:
            logger.warning(
                "DSSP: residue %d not found in chain %s of %s.",
                residue_pos, chain_id, pdb_path,
            )
            variant_record.set_feature_status(
                "dssp", False, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "secondary_structure", None,
                NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "secondary_structure_name", None,
                NullReason.NO_DATA_AVAILABLE,
            )
            return variant_record

        # DSSP tuple: (dssp index, amino acid, secondary structure, ...)
        # Index 2 is the secondary structure code
        ss_code = dssp[dssp_key][2]

        # Map to 3-state name; default to coil if code is unexpected
        ss_name = _DSSP_TO_NAME.get(ss_code, "coil")

        # Normalize unexpected codes to "-"
        if ss_code not in _DSSP_TO_NAME:
            logger.info(
                "DSSP: unexpected code %r at residue %d, defaulting to '-' (coil).",
                ss_code, residue_pos,
            )
            ss_code = "-"

        # Write results
        variant_record.secondary_structure = ss_code
        variant_record.secondary_structure_name = ss_name
        variant_record.set_feature_status("dssp", True)

        logger.info(
            "DSSP: residue %d — code=%s, name=%s",
            residue_pos, ss_code, ss_name,
        )

    except Exception as e:
        logger.warning(
            "DSSP failed for %s: %s",
            variant_record.variant_id or "unknown", e,
        )
        variant_record.set_feature_status(
            "dssp", False, NullReason.TOOL_CRASHED,
        )
        variant_record.set_with_reason(
            "secondary_structure", None, NullReason.TOOL_CRASHED,
        )
        variant_record.set_with_reason(
            "secondary_structure_name", None, NullReason.TOOL_CRASHED,
        )

    return variant_record
