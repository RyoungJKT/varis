"""FoldX — Primary protein stability calculation (ΔΔG).

Structural question: How much does this mutation destabilize the folded protein?
ΔΔG > 2 kcal/mol = damaging. This is the single strongest structural feature.

Priority: 4 (medium difficulty — needs license acceptance)
Fallback: PyRosetta, then DDGun, then skip ΔΔG entirely
Populates: ddg_foldx
"""
import logging
from varis.models.variant_record import VariantRecord
from varis.config import DDG_DAMAGING_THRESHOLD
logger = logging.getLogger(__name__)

def run_foldx(variant_record: VariantRecord) -> VariantRecord:
    """Calculate ΔΔG using FoldX BuildModel. Needs pdb_fixed_path from M2."""
    pass

def _prepare_foldx_config(pdb_path: str, position: int, ref_aa: str, alt_aa: str) -> str:
    """Generate FoldX configuration file for mutation."""
    pass

def _parse_foldx_output(output_path: str) -> float | None:
    """Parse FoldX Dif_ output file to extract ΔΔG value."""
    pass
