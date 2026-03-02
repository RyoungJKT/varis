"""HGVS Parser — Extracts gene, position, and amino acid change from HGVS notation.

Parses protein-level HGVS like "p.Arg1699Trp" into structured components:
residue position (1699), reference AA (Arg), alternate AA (Trp), charge change.
"""

import re
import logging
from varis.models.variant_record import VariantRecord, NullReason
from varis.config import AA_CHARGE, AA_THREE_TO_ONE, AA_ONE_TO_THREE

logger = logging.getLogger(__name__)

# Pattern: p.Arg1699Trp (three-letter codes)
HGVS_PROTEIN_PATTERN = re.compile(
    r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})"
)
# Pattern: p.R1699W (single-letter codes)
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
    hgvs = variant_record.hgvs_protein
    if not hgvs:
        for f in ("residue_position", "ref_amino_acid", "alt_amino_acid",
                  "ref_aa_single", "alt_aa_single", "charge_change"):
            variant_record.set_with_reason(f, None, NullReason.VALIDATION_FAILED)
        return variant_record

    # Try three-letter pattern first, then single-letter
    match = HGVS_PROTEIN_PATTERN.match(hgvs)
    if match:
        ref_three, position_str, alt_three = match.groups()
        ref_single = AA_THREE_TO_ONE.get(ref_three)
        alt_single = AA_THREE_TO_ONE.get(alt_three)
        if not ref_single or not alt_single:
            logger.warning(f"Unknown amino acid code in {hgvs}")
            for f in ("residue_position", "ref_amino_acid", "alt_amino_acid",
                      "ref_aa_single", "alt_aa_single", "charge_change"):
                variant_record.set_with_reason(f, None, NullReason.VALIDATION_FAILED)
            return variant_record
    else:
        match = HGVS_PROTEIN_SINGLE.match(hgvs)
        if match:
            ref_single, position_str, alt_single = match.groups()
            ref_three = AA_ONE_TO_THREE.get(ref_single)
            alt_three = AA_ONE_TO_THREE.get(alt_single)
            if not ref_three or not alt_three:
                logger.warning(f"Unknown amino acid code in {hgvs}")
                for f in ("residue_position", "ref_amino_acid", "alt_amino_acid",
                          "ref_aa_single", "alt_aa_single", "charge_change"):
                    variant_record.set_with_reason(f, None, NullReason.VALIDATION_FAILED)
                return variant_record
        else:
            logger.warning(f"Could not parse HGVS notation: {hgvs}")
            for f in ("residue_position", "ref_amino_acid", "alt_amino_acid",
                      "ref_aa_single", "alt_aa_single", "charge_change"):
                variant_record.set_with_reason(f, None, NullReason.VALIDATION_FAILED)
            return variant_record

    variant_record.residue_position = int(position_str)
    variant_record.ref_amino_acid = ref_three
    variant_record.alt_amino_acid = alt_three
    variant_record.ref_aa_single = ref_single
    variant_record.alt_aa_single = alt_single
    variant_record.charge_change = _calculate_charge_change(ref_three, alt_three)
    variant_record.input_notation_normalized = f"p.{ref_three}{position_str}{alt_three}"

    logger.info(
        f"Parsed HGVS: {hgvs} → position={position_str}, "
        f"{ref_three}({ref_single}) → {alt_three}({alt_single})"
    )
    return variant_record


def _calculate_charge_change(ref_aa: str, alt_aa: str) -> str:
    """Calculate the electrostatic charge change between two amino acids.

    Args:
        ref_aa: Reference amino acid (three-letter code).
        alt_aa: Alternate amino acid (three-letter code).

    Returns:
        Human-readable charge change, e.g., "positive → neutral".
    """
    charge_labels = {"+": "positive", "-": "negative", "0": "neutral"}
    ref_charge = AA_CHARGE.get(ref_aa, "0")
    alt_charge = AA_CHARGE.get(alt_aa, "0")
    ref_label = charge_labels.get(ref_charge, "neutral")
    alt_label = charge_labels.get(alt_charge, "neutral")
    if ref_charge == alt_charge:
        return f"no change ({ref_label})"
    return f"{ref_label} → {alt_label}"
