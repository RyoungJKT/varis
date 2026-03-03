"""Investigation response builder — maps VariantRecord to API response.

Converts the internal VariantRecord data structure into the 5-section
InvestigationResponse consumed by the React frontend.
"""

import logging
from typing import Optional

from varis.m6_platform.api.models import (
    InvestigationResponse,
    StructureSection,
    FeatureItem,
    PredictionSection,
    ExplanationItem,
    ProvenanceSection,
)

logger = logging.getLogger(__name__)


def build_investigation_response(record) -> InvestigationResponse:
    """Build an InvestigationResponse from a VariantRecord.

    Args:
        record: A VariantRecord instance (from varis.models.variant_record).

    Returns:
        InvestigationResponse with all 5 sections populated from available data.
    """
    return InvestigationResponse(
        variant_id=record.variant_id or f"{record.gene_symbol}_{record.hgvs_protein}",
        gene=record.gene_symbol or "",
        hgvs=record.hgvs_protein or "",
        clinvar_id=record.clinvar_id,
        structure=_build_structure(record),
        features=_build_features(record),
        prediction=_build_prediction(record),
        explanation=_build_explanation(record),
        evidence_tags=record.evidence_tags or [],
        provenance=_build_provenance(record),
    )


def _build_structure(record) -> StructureSection:
    """Build structure section from M2 data."""
    variant_id = record.variant_id or f"{record.gene_symbol}_{record.hgvs_protein}"
    pdb_url = f"/api/v1/structures/{variant_id}" if record.pdb_path else None

    return StructureSection(
        source=record.structure_source,
        chain="A",
        residue_index=record.structure_residue_position or record.residue_position,
        ref_aa=record.ref_aa_single,
        alt_aa=record.alt_aa_single,
        plddt_at_residue=record.mutation_site_plddt,
        plddt_mean=record.plddt_mean,
        confidence_bucket=record.mutation_site_confidence_bucket,
        coordinate_mapping_confidence=record.coordinate_mapping_confidence,
        uniprot_id=record.uniprot_id,
        pdb_url=pdb_url,
        normalization_warnings=record.normalization_warnings or [],
    )


def _build_features(record) -> list[FeatureItem]:
    """Build features list from structural and conservation data."""
    features = []
    null_reasons = record.null_reasons or {}

    # Structural features
    _add_feature(features, "ddg_mean", record.ddg_mean, "kcal/mol",
                 record.ddg_available, null_reasons.get("ddg_mean"),
                 _evidence_tag_for("ddg", record))
    _add_feature(features, "solvent_accessibility", record.solvent_accessibility_relative,
                 "relative", record.sasa_available, null_reasons.get("solvent_accessibility_relative"))
    _add_feature(features, "secondary_structure", None, None,
                 record.dssp_available, null_reasons.get("secondary_structure"),
                 value_override=record.secondary_structure_name)
    _add_feature(features, "contacts_wt", _to_float(record.contacts_wt), "count",
                 record.contacts_available, null_reasons.get("contacts_wt"))
    _add_feature(features, "hbonds_wt", _to_float(record.hbonds_wt), "count",
                 record.contacts_available, null_reasons.get("hbonds_wt"))
    _add_feature(features, "packing_density", record.packing_density, None,
                 record.contacts_available, null_reasons.get("packing_density"))

    # Domain
    _add_feature(features, "domain", None, None,
                 record.domain_available, null_reasons.get("domain_name"),
                 value_override=record.domain_name)

    # Conservation
    _add_feature(features, "conservation_score", record.conservation_score, None,
                 record.conservation_available, null_reasons.get("conservation_score"),
                 _evidence_tag_for("conservation", record))

    # Population
    _add_feature(features, "gnomad_frequency", record.gnomad_frequency, None,
                 record.gnomad_frequency is not None, null_reasons.get("gnomad_frequency"))

    # External scores
    _add_feature(features, "alphamissense_score", record.alphamissense_score, None,
                 record.alphamissense_score is not None,
                 null_reasons.get("alphamissense_score"))

    # pLDDT
    _add_feature(features, "mutation_site_plddt", record.mutation_site_plddt, None,
                 record.plddt_available, null_reasons.get("mutation_site_plddt"))

    return features


def _add_feature(features: list, name: str, value: Optional[float], units: Optional[str],
                 available: Optional[bool], null_reason: Optional[str] = None,
                 evidence_tag: Optional[str] = None,
                 value_override: Optional[str] = None) -> None:
    """Add a feature item to the list."""
    is_available = available if available is not None else (value is not None or value_override is not None)
    features.append(FeatureItem(
        name=name,
        value=value,
        units=units,
        evidence_tag=evidence_tag,
        available=is_available,
        null_reason=null_reason,
    ))


def _to_float(value) -> Optional[float]:
    """Convert int to float for FeatureItem compatibility."""
    return float(value) if value is not None else None


def _evidence_tag_for(feature_group: str, record) -> Optional[str]:
    """Get evidence tag for a feature if evidence criteria are met."""
    tags = record.evidence_tags or []
    if feature_group == "ddg" and "energetics_support" in tags:
        return "energetics_support"
    if feature_group == "conservation" and "computational_support" in tags:
        return "computational_support"
    return None


def _build_prediction(record) -> PredictionSection:
    """Build prediction section from M5 scoring data."""
    individual_scores = None
    if record.score_catboost is not None:
        individual_scores = {}
        if record.score_catboost is not None:
            individual_scores["catboost"] = record.score_catboost
        if record.score_xgboost is not None:
            individual_scores["xgboost"] = record.score_xgboost
        if record.score_lightgbm is not None:
            individual_scores["lightgbm"] = record.score_lightgbm

    return PredictionSection(
        score=record.score_ensemble,
        classification=record.classification,
        confidence_lower=record.confidence_lower,
        confidence_upper=record.confidence_upper,
        model_agreement=record.model_agreement,
        individual_scores=individual_scores,
        features_used=record.features_used,
        ensemble_version=record.ensemble_version,
    )


def _build_explanation(record) -> list[ExplanationItem]:
    """Build explanation section from SHAP top features.

    Returns items sorted by |shap| descending.
    """
    shap_features = record.shap_top_features
    if not shap_features:
        return []

    items = [
        ExplanationItem(
            feature=f["feature"],
            value=f.get("value"),
            shap=f["shap"],
        )
        for f in shap_features
    ]
    # Sort by |shap| descending
    items.sort(key=lambda x: abs(x.shap), reverse=True)
    return items


def _build_provenance(record) -> ProvenanceSection:
    """Build provenance section from metadata."""
    data_sources = []
    if record.clinvar_id:
        data_sources.append("ClinVar")
    if record.uniprot_id:
        data_sources.append("UniProt")
    if record.structure_source == "alphafold":
        data_sources.append("AlphaFold DB")
    elif record.structure_source == "esmfold":
        data_sources.append("ESMFold")
    if record.gnomad_frequency is not None:
        data_sources.append("gnomAD")
    if record.alphamissense_score is not None:
        data_sources.append("AlphaMissense")
    if record.conservation_score is not None:
        data_sources.append("Clustal Omega")

    return ProvenanceSection(
        data_sources=data_sources,
        modules_completed=record.modules_completed or [],
        modules_failed=record.modules_failed or [],
        pipeline_version=record.pipeline_version,
        investigation_timestamp=record.investigation_timestamp,
        processing_time_seconds=record.processing_time_seconds,
        tool_versions=record.tool_versions,
    )
