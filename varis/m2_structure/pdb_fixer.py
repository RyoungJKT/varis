"""PDB Fixer — Adds missing atoms/hydrogens using OpenMM PDBFixer.

Conditional execution: only runs when missing heavy atoms are detected.
AlphaFold PDBs are typically clean and will skip without modification.

CRITICAL DESIGN DECISION: This module NEVER reconstructs missing residues
or loops. That is modeling, not fixing, and creates false structural certainty.
We only add missing heavy atoms within existing residues and add hydrogens.

Requires: pdbfixer, openmm (optional dependencies — gracefully handles ImportError)
"""

import logging
from pathlib import Path

from varis.models.variant_record import NullReason, VariantRecord

logger = logging.getLogger(__name__)


def fix_structure(variant_record: VariantRecord) -> VariantRecord:
    """Conditionally repair a PDB structure by adding missing atoms/hydrogens.

    Only runs when missing heavy atoms are detected in existing residues.
    Never reconstructs missing residues or loops (that's modeling, not fixing).

    Args:
        variant_record: The shared VariantRecord with pdb_path set.

    Returns:
        The variant_record with pdb_fixed_path populated if fixing was applied,
        or unchanged if the structure was already clean or fixing was skipped.
    """
    # -------------------------------------------------------------------------
    # Step 1: Check pdb_path exists
    # -------------------------------------------------------------------------
    if variant_record.pdb_path is None:
        logger.info("No pdb_path set — skipping PDB fixer")
        variant_record.set_with_reason(
            "pdb_fixed_path", None, NullReason.UPSTREAM_DEPENDENCY_FAILED
        )
        return variant_record

    pdb_path = Path(variant_record.pdb_path)
    if not pdb_path.exists():
        logger.warning("PDB file does not exist: %s", pdb_path)
        variant_record.set_with_reason(
            "pdb_fixed_path", None, NullReason.UPSTREAM_DEPENDENCY_FAILED
        )
        return variant_record

    # -------------------------------------------------------------------------
    # Step 2: Try importing pdbfixer (lazy import — may not be installed)
    # -------------------------------------------------------------------------
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
    except ImportError:
        logger.info(
            "pdbfixer/openmm not installed — skipping PDB repair. "
            "Install with: conda install -c conda-forge pdbfixer openmm"
        )
        variant_record.set_with_reason(
            "pdb_fixed_path", None, NullReason.TOOL_MISSING
        )
        return variant_record

    # -------------------------------------------------------------------------
    # Steps 3-9: Run PDBFixer conditionally
    # -------------------------------------------------------------------------
    try:
        result = _run_fixer(variant_record, pdb_path, PDBFixer, PDBFile)
        return result
    except Exception as e:
        logger.warning(
            "PDBFixer crashed for %s: %s",
            variant_record.variant_id or "unknown",
            e,
        )
        variant_record.set_with_reason(
            "pdb_fixed_path", None, NullReason.TOOL_CRASHED
        )
        return variant_record


def _run_fixer(
    variant_record: VariantRecord,
    pdb_path: Path,
    PDBFixer: type,
    PDBFile: type,
) -> VariantRecord:
    """Internal: run the actual PDBFixer logic.

    Separated from fix_structure() so that the outer function handles
    all exception catching cleanly.

    Args:
        variant_record: The shared VariantRecord.
        pdb_path: Path to the input PDB file.
        PDBFixer: The PDBFixer class (passed to avoid re-importing).
        PDBFile: The openmm.app.PDBFile class (passed to avoid re-importing).

    Returns:
        The variant_record with pdb_fixed_path set if atoms were added.
    """
    # Step 3: Create PDBFixer instance
    fixer = PDBFixer(filename=str(pdb_path))

    # Step 4: Find missing residues, then CLEAR them — we do NOT fill gaps
    fixer.findMissingResidues()
    fixer.missingResidues = {}

    # Step 5: Find missing atoms in existing residues
    fixer.findMissingAtoms()
    missing_atoms = fixer.missingAtoms
    missing_terminals = fixer.missingTerminals

    n_missing_atoms = sum(len(atoms) for atoms in missing_atoms.values())
    n_missing_terminals = sum(len(atoms) for atoms in missing_terminals.values())
    total_missing = n_missing_atoms + n_missing_terminals

    # Step 6: If no missing atoms, return without fixing
    if total_missing == 0:
        logger.info(
            "No missing atoms detected in %s — structure is clean",
            pdb_path.name,
        )
        # Structure is clean: pdb_fixed_path stays None (or same as pdb_path)
        # but we do NOT set a null reason — this is the expected good case
        variant_record.pdb_fixed_path = None
        return variant_record

    # Step 7: Add missing atoms and hydrogens
    logger.info(
        "Adding %d missing atoms (%d heavy, %d terminal) to %s",
        total_missing,
        n_missing_atoms,
        n_missing_terminals,
        pdb_path.name,
    )
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)

    # Step 8: Save to .fixed.pdb
    fixed_path = pdb_path.with_suffix(".fixed.pdb")
    with open(fixed_path, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    logger.info("Fixed structure saved to %s", fixed_path)

    # Step 9: Update VariantRecord
    variant_record.pdb_fixed_path = str(fixed_path)

    if variant_record.preparation_steps is None:
        variant_record.preparation_steps = []
    variant_record.preparation_steps.append("added_missing_atoms")
    variant_record.preparation_steps.append("added_hydrogens")

    return variant_record
