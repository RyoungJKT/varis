"""Variant Record — the shared data structure for the entire Varis pipeline.

Every module reads from and writes to this structure. If a module fails,
its fields remain None AND the reason is recorded in null_reasons.

SCHEMA CONTRACT:
  - Every nullable field MUST have a corresponding entry in null_reasons when null
  - Reason codes: "not_attempted", "tool_missing", "tool_crashed", "license_unavailable",
    "low_confidence_structure", "rate_limited", "timed_out", "no_data_available",
    "intentionally_skipped", "coordinate_mapping_failed", "upstream_dependency_failed"
  - Schema version MUST be set on creation and never changed for that record
  - Feature availability flags (feature_available_*) are used by the ML model
    to distinguish "scientifically hard" from "pipeline broken"

Usage:
    record = create_variant_record("BRCA1", "p.Arg1699Trp")
    record = m1_ingestion.run(record)
    record = m2_structure.run(record)
    ...
"""

from __future__ import annotations

from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from typing import Optional
import json


# =============================================================================
# SCHEMA VERSION — Increment on any field addition, removal, or type change
# =============================================================================
RECORD_SCHEMA_VERSION = "1.0.0"


# =============================================================================
# NULL REASON CODES — Every None field must explain WHY it is None
# =============================================================================
class NullReason:
    """Standard reason codes for null fields.

    Every module MUST use these codes when setting a field to None.
    This makes missingness auditable and prevents the ML model from
    learning artifacts from pipeline failures vs. scientific difficulty.
    """
    NOT_ATTEMPTED = "not_attempted"
    TOOL_MISSING = "tool_missing"
    TOOL_CRASHED = "tool_crashed"
    LICENSE_UNAVAILABLE = "license_unavailable"
    LOW_CONFIDENCE_STRUCTURE = "low_confidence_structure"
    RATE_LIMITED = "rate_limited"
    TIMED_OUT = "timed_out"
    NO_DATA_AVAILABLE = "no_data_available"
    INTENTIONALLY_SKIPPED = "intentionally_skipped"
    COORDINATE_MAPPING_FAILED = "coordinate_mapping_failed"
    UPSTREAM_DEPENDENCY_FAILED = "upstream_dependency_failed"
    ISOFORM_MISMATCH = "isoform_mismatch"
    VALIDATION_FAILED = "validation_failed"


@dataclass
class VariantRecord:
    """The universal data container passed through every pipeline module.

    Every field defaults to None. Each module populates its own fields
    and leaves others untouched. The pipeline never crashes — if a module
    fails, its fields stay None and the reason is recorded in null_reasons.
    """

    # =========================================================================
    # SCHEMA METADATA — Record structure versioning
    # =========================================================================
    record_schema_version: str = RECORD_SCHEMA_VERSION
    feature_set_version: Optional[str] = None    # "lite_v1", "full_v1", "full_ddg_v1"

    # =========================================================================
    # INPUT — What the user provided
    # =========================================================================
    gene_symbol: Optional[str] = None
    hgvs_protein: Optional[str] = None
    hgvs_coding: Optional[str] = None
    hgvs_genomic: Optional[str] = None

    # =========================================================================
    # NORMALIZATION — Coordinate mapping and disambiguation (populated by M1)
    # =========================================================================
    canonical_transcript: Optional[str] = None        # MANE Select, e.g., "NM_007294.4"
    canonical_protein: Optional[str] = None           # e.g., "NP_009225.1"
    uniprot_residue_position: Optional[int] = None    # Position in UniProt canonical isoform
    structure_residue_position: Optional[int] = None  # Position in PDB/AlphaFold structure
    coordinate_mapping_method: Optional[str] = None   # "direct" / "sifts" / "blast_alignment"
    coordinate_mapping_confidence: Optional[str] = None  # "exact" / "high" / "low" / "failed"
    normalization_warnings: Optional[list[str]] = None
    input_notation_normalized: Optional[str] = None   # Canonical HGVS after normalization

    # =========================================================================
    # IDENTIFIERS — Cross-referenced IDs (populated by M1)
    # =========================================================================
    clinvar_id: Optional[str] = None
    clinvar_classification: Optional[str] = None
    clinvar_review_status: Optional[str] = None
    clinvar_conditions: Optional[list[str]] = None
    uniprot_id: Optional[str] = None
    ensembl_gene_id: Optional[str] = None
    ensembl_transcript_id: Optional[str] = None

    # =========================================================================
    # VARIANT DETAILS — Parsed from HGVS (populated by M1)
    # =========================================================================
    residue_position: Optional[int] = None
    ref_amino_acid: Optional[str] = None
    alt_amino_acid: Optional[str] = None
    ref_aa_single: Optional[str] = None
    alt_aa_single: Optional[str] = None
    charge_change: Optional[str] = None

    # =========================================================================
    # PROTEIN — Sequence and metadata (populated by M1)
    # =========================================================================
    protein_name: Optional[str] = None
    protein_sequence: Optional[str] = None
    protein_length: Optional[int] = None
    protein_function: Optional[str] = None

    # =========================================================================
    # POPULATION — Allele frequencies (populated by M1)
    # =========================================================================
    gnomad_frequency: Optional[float] = None
    gnomad_popmax: Optional[float] = None
    gnomad_homozygotes: Optional[int] = None

    # =========================================================================
    # EXTERNAL SCORES — Pre-computed predictions (populated by M1)
    # =========================================================================
    alphamissense_score: Optional[float] = None
    alphamissense_class: Optional[str] = None

    # =========================================================================
    # STRUCTURE — 3D structure data (populated by M2)
    # =========================================================================
    structure_source: Optional[str] = None
    pdb_path: Optional[str] = None
    pdb_fixed_path: Optional[str] = None
    plddt_score: Optional[float] = None
    plddt_mean: Optional[float] = None
    structure_resolution: Optional[str] = None

    # =========================================================================
    # STRUCTURAL FEATURES — From 3D analysis (populated by M3)
    # =========================================================================
    ddg_foldx: Optional[float] = None
    ddg_pyrosetta: Optional[float] = None
    ddg_mean: Optional[float] = None
    solvent_accessibility: Optional[float] = None
    burial_category: Optional[str] = None
    secondary_structure: Optional[str] = None
    secondary_structure_name: Optional[str] = None
    helix_disruption: Optional[bool] = None
    functional_site_distance: Optional[float] = None
    nearest_functional_site: Optional[str] = None
    hbonds_lost: Optional[int] = None
    contacts_changed: Optional[int] = None
    domain_name: Optional[str] = None
    domain_id: Optional[str] = None
    domain_criticality: Optional[str] = None

    # =========================================================================
    # CONSERVATION — Evolutionary analysis (populated by M4)
    # =========================================================================
    conservation_score: Optional[float] = None
    conservation_method: Optional[str] = None
    num_orthologs: Optional[int] = None
    position_entropy: Optional[float] = None
    conserved_across_mammals: Optional[bool] = None

    # =========================================================================
    # FEATURE AVAILABILITY — Explicit flags for ML model (populated by M3/M4)
    # These prevent the model from learning artifacts from missingness.
    # If ddg_available=False and ddg_missing_reason="tool_missing", the model
    # knows this is a pipeline limitation, not a scientific signal.
    # =========================================================================
    ddg_available: Optional[bool] = None
    ddg_missing_reason: Optional[str] = None
    sasa_available: Optional[bool] = None
    sasa_missing_reason: Optional[str] = None
    dssp_available: Optional[bool] = None
    dssp_missing_reason: Optional[str] = None
    conservation_available: Optional[bool] = None
    conservation_missing_reason: Optional[str] = None
    domain_available: Optional[bool] = None
    domain_missing_reason: Optional[str] = None
    contacts_available: Optional[bool] = None
    contacts_missing_reason: Optional[str] = None

    # =========================================================================
    # ML SCORING — Ensemble predictions (populated by M5)
    # =========================================================================
    score_ensemble: Optional[float] = None
    score_catboost: Optional[float] = None
    score_xgboost: Optional[float] = None
    score_lightgbm: Optional[float] = None
    confidence_lower: Optional[float] = None
    confidence_upper: Optional[float] = None
    classification: Optional[str] = None
    model_agreement: Optional[str] = None
    shap_top_features: Optional[list[dict]] = None
    ensemble_version: Optional[str] = None
    features_used: Optional[int] = None

    # =========================================================================
    # ACMG EVIDENCE — Suggested clinical classification codes (populated by M5)
    # NOTE: These are computational suggestions, not clinical adjudications.
    # =========================================================================
    acmg_codes: Optional[list[str]] = None
    acmg_pm1: Optional[bool] = None
    acmg_pm5: Optional[bool] = None
    acmg_pp3: Optional[bool] = None
    acmg_pp2: Optional[bool] = None
    acmg_ps3_proxy: Optional[bool] = None

    # =========================================================================
    # METADATA — Pipeline tracking
    # =========================================================================
    variant_id: Optional[str] = None
    pipeline_version: Optional[str] = None
    investigation_timestamp: Optional[str] = None
    modules_completed: Optional[list[str]] = None
    modules_failed: Optional[list[str]] = None
    tool_versions: Optional[dict] = None
    processing_time_seconds: Optional[float] = None
    published_to_varisdb: Optional[bool] = None
    submitted_to_clinvar: Optional[bool] = None

    # =========================================================================
    # NULL REASONS — Why each null field is null (populated by every module)
    # Keys = field names, values = NullReason codes.
    # Example: {"ddg_foldx": "tool_missing", "conservation_score": "rate_limited"}
    # =========================================================================
    null_reasons: Optional[dict[str, str]] = field(default_factory=dict)

    # =========================================================================
    # METHODS
    # =========================================================================

    def set_with_reason(self, field_name: str, value, reason: str = None) -> None:
        """Set a field value. If value is None, record the reason.

        This is the PREFERRED way to set fields. It ensures every None
        has an explanation in null_reasons.

        Args:
            field_name: Name of the field to set.
            value: Value to set (can be None).
            reason: Required if value is None. Must be a NullReason constant.
        """
        setattr(self, field_name, value)
        if self.null_reasons is None:
            self.null_reasons = {}
        if value is None and reason:
            self.null_reasons[field_name] = reason
        elif value is not None and field_name in self.null_reasons:
            del self.null_reasons[field_name]

    def set_feature_status(self, feature_group: str, available: bool,
                           reason: str = None) -> None:
        """Set feature availability flag and reason for a feature group.

        Args:
            feature_group: One of "ddg", "sasa", "dssp", "conservation",
                          "domain", "contacts".
            available: Whether the feature was successfully computed.
            reason: If not available, why. Must be a NullReason constant.
        """
        setattr(self, f"{feature_group}_available", available)
        if not available and reason:
            setattr(self, f"{feature_group}_missing_reason", reason)
        elif available:
            setattr(self, f"{feature_group}_missing_reason", None)

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return asdict(self)

    def to_json(self, indent: int = 2) -> str:
        """Serialize to JSON string."""
        return json.dumps(self.to_dict(), indent=indent, default=str)

    @classmethod
    def from_dict(cls, data: dict) -> "VariantRecord":
        """Create VariantRecord from dictionary, ignoring unknown keys."""
        valid_fields = {f.name for f in cls.__dataclass_fields__.values()}
        filtered = {k: v for k, v in data.items() if k in valid_fields}
        return cls(**filtered)

    @classmethod
    def from_json(cls, json_str: str) -> "VariantRecord":
        """Create VariantRecord from JSON string."""
        return cls.from_dict(json.loads(json_str))

    def save(self, path: str) -> None:
        """Save to JSON file."""
        with open(path, "w") as f:
            f.write(self.to_json())

    @classmethod
    def load(cls, path: str) -> "VariantRecord":
        """Load from JSON file."""
        with open(path, "r") as f:
            return cls.from_json(f.read())

    def get_structural_features(self) -> dict:
        """Extract structural feature dict for ML model input.

        Returns the 15 features the ML model uses, with None for missing.
        The gradient-boosted tree ensemble handles None natively.
        """
        return {
            "ddg_foldx": self.ddg_foldx,
            "ddg_pyrosetta": self.ddg_pyrosetta,
            "solvent_accessibility": self.solvent_accessibility,
            "secondary_structure": self.secondary_structure_name,
            "functional_site_distance": self.functional_site_distance,
            "hbonds_lost": self.hbonds_lost,
            "contacts_changed": self.contacts_changed,
            "domain_name": self.domain_name,
            "conservation_score": self.conservation_score,
            "plddt_score": self.plddt_score,
            "gnomad_frequency": self.gnomad_frequency,
            "alphamissense_score": self.alphamissense_score,
            "charge_change": self.charge_change,
            "burial_category": self.burial_category,
            "helix_disruption": self.helix_disruption,
        }

    def get_feature_availability_flags(self) -> dict:
        """Extract feature availability flags for ML model input.

        These let the model distinguish 'scientifically hard variant'
        from 'pipeline broken'. Include alongside structural features.
        """
        return {
            "ddg_available": self.ddg_available,
            "sasa_available": self.sasa_available,
            "dssp_available": self.dssp_available,
            "conservation_available": self.conservation_available,
            "domain_available": self.domain_available,
            "contacts_available": self.contacts_available,
        }

    def get_ml_features(self) -> dict:
        """Get complete feature set for ML: structural + availability flags."""
        features = self.get_structural_features()
        features.update(self.get_feature_availability_flags())
        return features

    def count_available_features(self) -> int:
        """Count how many ML features are non-null."""
        return sum(1 for v in self.get_structural_features().values()
                   if v is not None)

    def mark_module_completed(self, module: str) -> None:
        """Record that a module completed successfully."""
        if self.modules_completed is None:
            self.modules_completed = []
        if module not in self.modules_completed:
            self.modules_completed.append(module)

    def mark_module_failed(self, module: str, reason: str = None) -> None:
        """Record that a module (or sub-tool) failed, with optional reason."""
        if self.modules_failed is None:
            self.modules_failed = []
        if module not in self.modules_failed:
            self.modules_failed.append(module)

    def validate(self) -> list[str]:
        """Validate the record for common issues. Returns list of warnings."""
        warnings = []
        if self.record_schema_version != RECORD_SCHEMA_VERSION:
            warnings.append(
                f"Schema version mismatch: record={self.record_schema_version}, "
                f"current={RECORD_SCHEMA_VERSION}"
            )
        feature_fields = [
            "ddg_foldx", "ddg_pyrosetta", "solvent_accessibility",
            "secondary_structure", "conservation_score", "domain_name",
        ]
        if self.null_reasons is None:
            self.null_reasons = {}
        for f in feature_fields:
            val = getattr(self, f, None)
            if val is None and f not in self.null_reasons:
                if self.modules_completed and any(
                    m.startswith(("M3", "M4")) for m in self.modules_completed
                ):
                    warnings.append(
                        f"Field '{f}' is null but has no reason code")
        if self.structure_residue_position is None and self.pdb_path is not None:
            warnings.append(
                "Structure available but residue position not mapped")
        return warnings


def create_variant_record(gene_symbol: str, hgvs_protein: str) -> VariantRecord:
    """Create a new VariantRecord with input fields populated.

    Args:
        gene_symbol: Gene name, e.g., "BRCA1"
        hgvs_protein: Protein-level HGVS notation, e.g., "p.Arg1699Trp"

    Returns:
        A new VariantRecord with input fields set, schema version stamped,
        and all other fields None with reason "not_attempted".
    """
    return VariantRecord(
        gene_symbol=gene_symbol,
        hgvs_protein=hgvs_protein,
        variant_id=f"{gene_symbol}_{hgvs_protein}",
        record_schema_version=RECORD_SCHEMA_VERSION,
        investigation_timestamp=datetime.now(timezone.utc).isoformat(),
        modules_completed=[],
        modules_failed=[],
        null_reasons={},
    )
