"""PDB Fixer — Adds missing atoms/hydrogens using OpenMM PDBFixer."""
import logging
from varis.models.variant_record import VariantRecord
logger = logging.getLogger(__name__)

def fix_structure(variant_record: VariantRecord) -> VariantRecord:
    """Repair PDB with PDBFixer. Sets pdb_fixed_path."""
    pass
