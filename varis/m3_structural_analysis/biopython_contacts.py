"""BioPython Contacts — Counts WT local environment at mutation site.

Measures the wild-type local environment: heavy-atom contacts within 4.5A,
hydrogen bonds, and packing density. Does NOT model the mutant -- that would
require FoldX/PyRosetta to be scientifically valid.

Definition: Heavy atoms only (element != H). Contact threshold: 4.5A.
H-bond: donor-acceptor distance <=3.5A with N/O/S atoms.

Populates: contacts_wt, hbonds_wt, packing_density.
Sets: contacts_available, contacts_missing_reason.
"""

import logging
from typing import Optional

from varis.models.variant_record import NullReason, VariantRecord

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
_CONTACT_THRESHOLD = 4.5   # Angstroms
_HBOND_THRESHOLD = 3.5     # Angstroms
_HBOND_ELEMENTS = {"N", "O", "S"}


def _get_pdb_path(record: VariantRecord) -> Optional[str]:
    """Return the best available PDB path (prefer fixed structure).

    Args:
        record: The VariantRecord with structure paths.

    Returns:
        Path string to the PDB file, or None if no structure available.
    """
    return record.pdb_fixed_path or record.pdb_path


def run_contacts(variant_record: VariantRecord) -> VariantRecord:
    """Count heavy-atom contacts, H-bonds, and packing density at mutation site.

    Uses BioPython's NeighborSearch for spatial queries on the wild-type
    structure. Only heavy atoms (element != H) are considered.

    Args:
        variant_record: The shared VariantRecord (needs pdb_path and
            residue_position).

    Returns:
        The variant_record with contacts fields populated (or None on failure).
    """
    pdb_path = _get_pdb_path(variant_record)
    if pdb_path is None:
        logger.warning("Contacts skipped: no PDB structure available.")
        variant_record.set_feature_status(
            "contacts", False, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.set_with_reason(
            "contacts_wt", None, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.set_with_reason(
            "hbonds_wt", None, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.set_with_reason(
            "packing_density", None, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        return variant_record

    residue_pos = variant_record.residue_position
    if residue_pos is None:
        logger.warning("Contacts skipped: no residue_position in variant record.")
        variant_record.set_feature_status(
            "contacts", False, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.set_with_reason(
            "contacts_wt", None, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.set_with_reason(
            "hbonds_wt", None, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        variant_record.set_with_reason(
            "packing_density", None, NullReason.UPSTREAM_DEPENDENCY_FAILED,
        )
        return variant_record

    try:
        from Bio.PDB import PDBParser, NeighborSearch
    except ImportError:
        logger.warning("Contacts skipped: BioPython not installed.")
        variant_record.set_feature_status(
            "contacts", False, NullReason.TOOL_MISSING,
        )
        variant_record.set_with_reason(
            "contacts_wt", None, NullReason.TOOL_MISSING,
        )
        variant_record.set_with_reason(
            "hbonds_wt", None, NullReason.TOOL_MISSING,
        )
        variant_record.set_with_reason(
            "packing_density", None, NullReason.TOOL_MISSING,
        )
        return variant_record

    try:
        # Parse the PDB structure
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("variant", pdb_path)
        model = structure[0]

        # Get first chain
        chains = list(model.get_chains())
        if not chains:
            logger.warning("Contacts: no chains found in structure %s.", pdb_path)
            variant_record.set_feature_status(
                "contacts", False, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "contacts_wt", None, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "hbonds_wt", None, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "packing_density", None, NullReason.NO_DATA_AVAILABLE,
            )
            return variant_record

        chain = chains[0]

        # Find target residue at residue_position
        target_residue = None
        for residue in chain.get_residues():
            if residue.get_id()[1] == residue_pos:
                target_residue = residue
                break

        if target_residue is None:
            logger.warning(
                "Contacts: residue %d not found in chain %s of %s.",
                residue_pos, chain.id, pdb_path,
            )
            variant_record.set_feature_status(
                "contacts", False, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "contacts_wt", None, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "hbonds_wt", None, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "packing_density", None, NullReason.NO_DATA_AVAILABLE,
            )
            return variant_record

        # Collect all heavy atoms in the entire structure (element != "H")
        all_heavy_atoms = [
            atom for atom in model.get_atoms()
            if atom.element.strip() != "H"
        ]

        if not all_heavy_atoms:
            logger.warning("Contacts: no heavy atoms found in %s.", pdb_path)
            variant_record.set_feature_status(
                "contacts", False, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "contacts_wt", None, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "hbonds_wt", None, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "packing_density", None, NullReason.NO_DATA_AVAILABLE,
            )
            return variant_record

        # Get heavy atoms in the target residue
        target_heavy_atoms = [
            atom for atom in target_residue.get_atoms()
            if atom.element.strip() != "H"
        ]

        if not target_heavy_atoms:
            logger.warning(
                "Contacts: no heavy atoms in target residue %d.", residue_pos,
            )
            variant_record.set_feature_status(
                "contacts", False, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "contacts_wt", None, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "hbonds_wt", None, NullReason.NO_DATA_AVAILABLE,
            )
            variant_record.set_with_reason(
                "packing_density", None, NullReason.NO_DATA_AVAILABLE,
            )
            return variant_record

        # Build set of target atom ids for exclusion
        target_atom_ids = {id(atom) for atom in target_heavy_atoms}

        # Build NeighborSearch with all heavy atoms
        ns = NeighborSearch(all_heavy_atoms)

        # Count contacts: neighboring atoms NOT in the target residue
        contact_atom_ids: set[int] = set()
        for target_atom in target_heavy_atoms:
            neighbors = ns.search(target_atom.get_vector().get_array(),
                                  _CONTACT_THRESHOLD, level="A")
            for neighbor in neighbors:
                neighbor_id = id(neighbor)
                if neighbor_id not in target_atom_ids:
                    contact_atom_ids.add(neighbor_id)

        contacts_wt = len(contact_atom_ids)

        # Count H-bonds: contacts where both atoms are N/O/S and distance <= 3.5A
        # Use frozenset pairs to avoid double counting
        hbond_pairs: set[frozenset[int]] = set()
        for target_atom in target_heavy_atoms:
            if target_atom.element.strip() not in _HBOND_ELEMENTS:
                continue
            neighbors = ns.search(target_atom.get_vector().get_array(),
                                  _HBOND_THRESHOLD, level="A")
            for neighbor in neighbors:
                neighbor_id = id(neighbor)
                if neighbor_id in target_atom_ids:
                    continue
                if neighbor.element.strip() not in _HBOND_ELEMENTS:
                    continue
                pair = frozenset((id(target_atom), neighbor_id))
                hbond_pairs.add(pair)

        hbonds_wt = len(hbond_pairs)

        # Packing density = contacts / number of target heavy atoms
        packing_density = contacts_wt / len(target_heavy_atoms)

        # Write results
        variant_record.contacts_wt = contacts_wt
        variant_record.hbonds_wt = hbonds_wt
        variant_record.packing_density = packing_density
        variant_record.set_feature_status("contacts", True)

        logger.info(
            "Contacts: residue %d — contacts=%d, hbonds=%d, "
            "packing_density=%.2f",
            residue_pos, contacts_wt, hbonds_wt, packing_density,
        )

    except Exception as e:
        logger.warning(
            "Contacts failed for %s: %s",
            variant_record.variant_id or "unknown", e,
        )
        variant_record.set_feature_status(
            "contacts", False, NullReason.TOOL_CRASHED,
        )
        variant_record.set_with_reason(
            "contacts_wt", None, NullReason.TOOL_CRASHED,
        )
        variant_record.set_with_reason(
            "hbonds_wt", None, NullReason.TOOL_CRASHED,
        )
        variant_record.set_with_reason(
            "packing_density", None, NullReason.TOOL_CRASHED,
        )

    return variant_record
