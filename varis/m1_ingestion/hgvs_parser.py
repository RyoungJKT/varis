"""HGVS Parser — Extracts gene, position, and amino acid change from HGVS notation.

Parses protein-level HGVS like "p.Arg1699Trp" into structured components:
residue position (1699), reference AA (Arg), alternate AA (Trp), charge change.
"""

import re
import logging
from varis.models.variant_record import VariantRecord
from varis.config import AA_CHARGE, AA_THREE_TO_ONE

logger = logging.getLogger(__name__)

# Pattern: p.Arg1699Trp or p.R1699W
HGVS_PROTEIN_PATTERN = re.compile(
    r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})"
)
HGVS_PROTEIN_SINGLE = re.compile(
    r"p\.([A-Z])(\d+)([A-Z])"
)


def parse_hgvs(variant_record: VariantRecord) -> VariantRecord:
    """Parse HGVS protein notation and populate variant detail fields.

    Args:
        variant_record: Must have hgvs_protein set (e.g., "p.Arg1699Trp").

    Returns:
        VariantRecord with residue_position, ref/alt amino acids, and charge change.
    """
    pass


def _calculate_charge_change(ref_aa: str, alt_aa: str) -> str:
    """Calculate the electrostatic charge change between two amino acids.

    Args:
        ref_aa: Reference amino acid (three-letter code).
        alt_aa: Alternate amino acid (three-letter code).

    Returns:
        Human-readable charge change, e.g., "+ve → neutral".
    """
    pass
