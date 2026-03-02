"""PyRosetta — Independent second ΔΔG calculation with backbone relaxation.

Structural question: Same as FoldX, but with conformational sampling.
Provides a second opinion on stability. When FoldX and PyRosetta agree, confidence is high.

Priority: 5 (medium-hard — needs license)
Populates: ddg_pyrosetta
"""
import logging
from varis.models.variant_record import VariantRecord
logger = logging.getLogger(__name__)

def run_pyrosetta(variant_record: VariantRecord) -> VariantRecord:
    """Calculate ΔΔG using PyRosetta with backbone relaxation."""
    pass
