"""ACMG Mapper — Maps structural evidence to suggested ACMG classification codes.

Translates Varis's structural findings into the clinical evidence framework
that genetic counselors use. AlphaMissense only provides PP3.
Varis suggests: PM1, PM5, PP3, PP2, PS3-proxy.

CRITICAL: Varis SUGGESTS evidence codes to support professional review. It does
NOT replace ACMG adjudication by qualified variant curation teams. All codes —
particularly PS3-proxy (computational structural evidence, not wet-lab functional
data) — must be clearly labelled as computational suggestions in every report,
API response, and user-facing output.

ACMG codes supported:
  PP3: Computational prediction supports pathogenic (ensemble score)
  PM1: Mutation in critical functional domain (HMMER)
  PM5: Known pathogenic variant at same position (ClinVar cross-ref)
  PS3-proxy: Structural damage evidence (ΔΔG destabilization) — NOT full PS3
  PP2: Low benign rate in this gene (gnomAD)
"""
import logging
from varis.models.variant_record import VariantRecord
from varis.config import DDG_DAMAGING_THRESHOLD, GNOMAD_RARE_THRESHOLD
logger = logging.getLogger(__name__)

def map_acmg_codes(variant_record: VariantRecord) -> VariantRecord:
    """Assign ACMG evidence codes based on structural investigation results.

    Each code has explicit criteria. Only assigned when evidence supports it.
    """
    pass

def _check_pp3(variant_record: VariantRecord) -> bool:
    """PP3: Multiple computational predictions support pathogenic."""
    pass

def _check_pm1(variant_record: VariantRecord) -> bool:
    """PM1: Located in critical functional domain without benign variation."""
    pass

def _check_ps3_proxy(variant_record: VariantRecord) -> bool:
    """PS3-proxy: ΔΔG > threshold indicates functional damage."""
    pass

def _check_pp2(variant_record: VariantRecord) -> bool:
    """PP2: Gene has low rate of benign missense variation."""
    pass
