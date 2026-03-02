"""FreeSASA — Solvent Accessible Surface Area calculation.

Structural question: How buried is the mutated residue?
Buried residues (low SASA) in the hydrophobic core are structurally critical.
Mutations there are devastating. Surface residues are more tolerant.

Priority: 1 (build first — easy to install, always works)
Populates: solvent_accessibility, burial_category
"""
import logging
from varis.models.variant_record import VariantRecord
from varis.config import DDG_DAMAGING_THRESHOLD
logger = logging.getLogger(__name__)

def run_freesasa(variant_record: VariantRecord) -> VariantRecord:
    """Calculate solvent accessibility for the mutated residue.
    Classifies as core (<10%), interface (10-40%), or surface (>40%)."""
    pass
