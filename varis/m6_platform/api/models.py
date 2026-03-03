"""Pydantic response models for the VarisDB API.

Defines all request/response schemas for FastAPI endpoints.
These models drive automatic validation and OpenAPI documentation.
"""

from typing import Optional
from pydantic import BaseModel, Field


# =============================================================================
# INVESTIGATION RESPONSE — 5-section payload for React frontend
# =============================================================================

class StructureSection(BaseModel):
    """3D structure metadata for the Mol* viewer."""

    source: Optional[str] = None
    chain: Optional[str] = None
    residue_index: Optional[int] = None
    ref_aa: Optional[str] = None
    alt_aa: Optional[str] = None
    plddt_at_residue: Optional[float] = None
    plddt_mean: Optional[float] = None
    confidence_bucket: Optional[str] = None
    coordinate_mapping_confidence: Optional[str] = None
    uniprot_id: Optional[str] = None
    normalization_warnings: list[str] = Field(default_factory=list)


class FeatureItem(BaseModel):
    """A single structural/conservation feature with availability info."""

    name: str
    value: Optional[float] = None
    units: Optional[str] = None
    evidence_tag: Optional[str] = None
    available: bool = True
    null_reason: Optional[str] = None


class PredictionSection(BaseModel):
    """ML ensemble prediction with confidence bounds."""

    score: Optional[float] = None
    classification: Optional[str] = None
    confidence_lower: Optional[float] = None
    confidence_upper: Optional[float] = None
    model_agreement: Optional[str] = None
    individual_scores: Optional[dict[str, float]] = None
    features_used: Optional[int] = None
    ensemble_version: Optional[str] = None


class ExplanationItem(BaseModel):
    """A single SHAP explanation entry."""

    feature: str
    value: Optional[float] = None
    shap: float


class ProvenanceSection(BaseModel):
    """Data source tracking and audit trail."""

    data_sources: list[str] = Field(default_factory=list)
    modules_completed: list[str] = Field(default_factory=list)
    modules_failed: list[str] = Field(default_factory=list)
    pipeline_version: Optional[str] = None
    investigation_timestamp: Optional[str] = None
    processing_time_seconds: Optional[float] = None
    tool_versions: Optional[dict[str, str]] = None


class InvestigationResponse(BaseModel):
    """Full investigation payload — consumed by the React frontend.

    Contains 5 sections: structure, features, prediction, explanation, provenance.
    """

    variant_id: str
    gene: str
    hgvs: str
    clinvar_id: Optional[str] = None
    structure: StructureSection = Field(default_factory=StructureSection)
    features: list[FeatureItem] = Field(default_factory=list)
    prediction: PredictionSection = Field(default_factory=PredictionSection)
    explanation: list[ExplanationItem] = Field(default_factory=list)
    evidence_tags: list[str] = Field(default_factory=list)
    provenance: ProvenanceSection = Field(default_factory=ProvenanceSection)


# =============================================================================
# SEARCH RESPONSE
# =============================================================================

class VariantSummary(BaseModel):
    """Lightweight variant record for search results."""

    variant_id: str
    gene: str
    hgvs_protein: Optional[str] = None
    classification: Optional[str] = None
    score: Optional[float] = None
    clinvar_id: Optional[str] = None
    investigation_timestamp: Optional[str] = None


class SearchResponse(BaseModel):
    """Paginated search results."""

    results: list[VariantSummary] = Field(default_factory=list)
    total: int = 0
    page: int = 1
    limit: int = 20


# =============================================================================
# JOB STATUS
# =============================================================================

class JobStatusResponse(BaseModel):
    """Status of an async investigation job."""

    job_id: str
    status: str  # "queued" | "running" | "succeeded" | "failed"
    variant_id: Optional[str] = None
    current_step: Optional[str] = None
    error_message: Optional[str] = None
    created_at: Optional[str] = None
    updated_at: Optional[str] = None


# =============================================================================
# REQUEST MODELS
# =============================================================================

class InvestigationRequest(BaseModel):
    """Request to submit a new variant investigation."""

    gene: str = Field(..., min_length=1, max_length=50)
    hgvs: str = Field(..., min_length=1, max_length=100)


# =============================================================================
# STATS RESPONSE
# =============================================================================

class StatsResponse(BaseModel):
    """Database statistics."""

    total_variants: int = 0
    total_pathogenic: int = 0
    total_benign: int = 0
    total_uncertain: int = 0
    total_jobs: int = 0
    genes_covered: int = 0
