"""DSSP — Secondary structure assignment from 3D coordinates.

Structural question: What secondary structure element does the mutation sit in?
A mutation that breaks an alpha-helix is structurally catastrophic.

Priority: 2 (easy, via mkdssp or mdtraj)
Populates: secondary_structure, secondary_structure_name, helix_disruption
"""
import logging
from varis.models.variant_record import VariantRecord
logger = logging.getLogger(__name__)

def run_dssp(variant_record: VariantRecord) -> VariantRecord:
    """Assign secondary structure at mutation site. H=helix, E=sheet, C=coil."""
    pass
