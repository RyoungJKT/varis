"""Evidence Mapper — Maps computational evidence to suggested tags.

Evidence tags are inspired by ACMG criteria but are NOT clinical
adjudications. Each tag requires multiple independent signals to
prevent circular reasoning ("the model proves itself").

Tags:
  computational_support (PP3-like): conservation > 0.9 AND (buried core OR high score)
  rarity_evidence (PM2-like): gnomAD < 0.0001 or absent
  energetics_support: DDG > 2.0 kcal/mol
  domain_context (PM1-like): critical domain
"""
import logging

from varis.config import DDG_DAMAGING_THRESHOLD, GNOMAD_RARE_THRESHOLD
from varis.models.variant_record import VariantRecord

logger = logging.getLogger(__name__)


def map_evidence_tags(variant_record: VariantRecord) -> VariantRecord:
    """Assign evidence tags based on structural investigation results.

    Each tag has explicit criteria requiring multiple independent signals.

    Args:
        variant_record: Record with features and scores populated.

    Returns:
        VariantRecord with evidence_* fields and evidence_tags populated.
    """
    tags = []

    # Computational support: conservation > 0.9 AND (buried core OR high score)
    comp_support = _check_computational_support(variant_record)
    variant_record.evidence_computational_support = comp_support
    if comp_support:
        tags.append("computational_support")

    # Rarity: gnomAD < 0.0001 or absent
    rarity = _check_rarity(variant_record)
    variant_record.evidence_rarity = rarity
    if rarity:
        tags.append("rarity_evidence")

    # Energetics: DDG > 2.0
    energetics = _check_energetics(variant_record)
    variant_record.evidence_energetics = energetics
    if energetics:
        tags.append("energetics_support")

    # Domain context: critical domain
    domain = _check_domain_context(variant_record)
    variant_record.evidence_domain_context = domain
    if domain:
        tags.append("domain_context")

    variant_record.evidence_tags = tags
    return variant_record


def _check_computational_support(variant_record: VariantRecord) -> bool:
    """PP3-like: Multiple computational signals support pathogenicity.

    Requires conservation > 0.9 AND at least one of:
    - Buried in core (burial_category == "core")
    - High ensemble score (> 0.8)

    Args:
        variant_record: Record with conservation and structural data.

    Returns:
        True if computational support criteria met.
    """
    conservation = variant_record.conservation_score
    if conservation is None or conservation <= 0.9:
        return False

    buried = variant_record.burial_category == "core"
    high_score = (
        variant_record.score_ensemble is not None
        and variant_record.score_ensemble > 0.8
    )

    return buried or high_score


def _check_rarity(variant_record: VariantRecord) -> bool:
    """PM2-like: Variant is rare or absent from population databases.

    Args:
        variant_record: Record with gnomAD data.

    Returns:
        True if gnomAD frequency < threshold or absent.
    """
    freq = variant_record.gnomad_frequency
    if freq is None:
        return True  # Absent from gnomAD = rare
    return freq < GNOMAD_RARE_THRESHOLD


def _check_energetics(variant_record: VariantRecord) -> bool:
    """Energetics support: DDG indicates structural destabilization.

    Args:
        variant_record: Record with DDG data.

    Returns:
        True if DDG > threshold.
    """
    ddg = variant_record.ddg_mean
    if ddg is None:
        return False
    return ddg > DDG_DAMAGING_THRESHOLD


def _check_domain_context(variant_record: VariantRecord) -> bool:
    """PM1-like: Variant is in a critical functional domain.

    Args:
        variant_record: Record with domain data.

    Returns:
        True if domain_criticality is "critical".
    """
    return variant_record.domain_criticality == "critical"
