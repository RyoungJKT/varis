"""Tests for ClinVar submission formatter and submitter."""
import os
import pytest
from unittest.mock import patch

from varis.models.variant_record import VariantRecord, create_variant_record
from varis.m6_platform.api.clinvar_submitter import (
    format_clinvar_submission,
    submit_to_clinvar,
)


# =========================================================================
# FIXTURES
# =========================================================================


@pytest.fixture
def pathogenic_record() -> VariantRecord:
    """BRCA1 p.Arg1699Trp — high-confidence likely pathogenic with full evidence."""
    record = create_variant_record("BRCA1", "p.Arg1699Trp")
    record.variant_id = "BRCA1_p.Arg1699Trp"
    record.canonical_transcript = "NM_007294.4"
    record.hgvs_coding = "c.5095C>T"
    record.hgvs_protein = "p.Arg1699Trp"
    record.reference_build = "GRCh38"
    record.clinvar_chrom = "17"
    record.clinvar_pos = 43057051
    record.clinvar_ref = "G"
    record.clinvar_alt = "A"
    record.clinvar_id = "VCV000055361"
    record.residue_position = 1699
    record.ref_amino_acid = "Arg"
    record.alt_amino_acid = "Trp"
    record.gene_symbol = "BRCA1"
    record.uniprot_id = "P38398"
    record.gnomad_frequency = 0.00003

    # Structural features
    record.ddg_mean = 3.9
    record.burial_category = "core"
    record.secondary_structure_name = "helix"
    record.domain_name = "BRCT"
    record.domain_criticality = "critical"
    record.conservation_score = 1.0

    # ML scoring
    record.score_ensemble = 0.92
    record.confidence_lower = 0.86
    record.confidence_upper = 0.97
    record.classification = "likely_pathogenic"
    record.model_agreement = "high"
    record.ensemble_version = "v1.2"
    record.features_used = 15

    # Evidence tags
    record.evidence_tags = [
        "computational_support",
        "rarity_evidence",
        "energetics_support",
        "domain_context",
    ]

    # Pipeline metadata
    record.modules_completed = ["M1", "M2", "M3", "M4", "M5"]
    record.pipeline_version = "v1.0"

    return record


@pytest.fixture
def benign_record() -> VariantRecord:
    """BRCA1 p.Lys1183Arg — high-confidence likely benign."""
    record = create_variant_record("BRCA1", "p.Lys1183Arg")
    record.variant_id = "BRCA1_p.Lys1183Arg"
    record.canonical_transcript = "NM_007294.4"
    record.hgvs_coding = "c.3548A>G"
    record.hgvs_protein = "p.Lys1183Arg"
    record.reference_build = "GRCh38"
    record.clinvar_chrom = "17"
    record.clinvar_pos = 43093220
    record.clinvar_ref = "T"
    record.clinvar_alt = "C"
    record.gene_symbol = "BRCA1"

    # ML scoring
    record.score_ensemble = 0.15
    record.confidence_lower = 0.08
    record.confidence_upper = 0.22
    record.classification = "likely_benign"
    record.model_agreement = "high"
    record.ensemble_version = "v1.2"
    record.features_used = 12
    record.evidence_tags = []

    # Pipeline metadata
    record.modules_completed = ["M1", "M2", "M3", "M4", "M5"]
    record.pipeline_version = "v1.0"

    return record


@pytest.fixture
def incomplete_record() -> VariantRecord:
    """Only M1 completed — no classification, no scoring."""
    record = create_variant_record("BRCA1", "p.Arg1699Trp")
    record.variant_id = "BRCA1_p.Arg1699Trp"
    record.gene_symbol = "BRCA1"
    record.modules_completed = ["M1"]
    return record


# =========================================================================
# TESTS — format_clinvar_submission
# =========================================================================


class TestFormatPathogenicSubmission:
    """Test formatting a likely pathogenic variant for ClinVar."""

    def test_format_pathogenic_submission(self, pathogenic_record: VariantRecord) -> None:
        """Valid pathogenic record produces submission dict with correct fields."""
        result = format_clinvar_submission(pathogenic_record)

        assert result is not None
        assert result["record_status"] == "novel"
        assert result["variant_id"] == "BRCA1_p.Arg1699Trp"
        assert result["gene_symbol"] == "BRCA1"
        assert result["hgvs_coding"] == "NM_007294.4:c.5095C>T"
        assert result["hgvs_protein"] == "p.Arg1699Trp"
        assert result["reference_build"] == "GRCh38"
        assert result["chromosome"] == "17"
        assert result["position"] == 43057051
        assert result["ref_allele"] == "G"
        assert result["alt_allele"] == "A"
        assert result["clinical_significance"] == "Likely pathogenic"
        assert result["collection_method"] == "research"
        assert result["allele_origin"] == "germline"
        assert result["affected_status"] == "not provided"
        assert result["submitter"] == "Russell Genetics"
        assert result["local_id"] == "BRCA1_p.Arg1699Trp"
        assert result["clinvar_accession"] == "VCV000055361"
        # date_last_evaluated should be a YYYY-MM-DD string
        assert len(result["date_last_evaluated"]) == 10
        assert "-" in result["date_last_evaluated"]


class TestFormatBenignSubmission:
    """Test formatting a likely benign variant for ClinVar."""

    def test_format_benign_submission(self, benign_record: VariantRecord) -> None:
        """Benign classification maps to ClinVar term 'Likely benign'."""
        result = format_clinvar_submission(benign_record)

        assert result is not None
        assert result["clinical_significance"] == "Likely benign"
        assert result["gene_symbol"] == "BRCA1"


class TestEvidenceSummary:
    """Test that the evidence comment includes required information."""

    def test_format_includes_evidence_summary(self, pathogenic_record: VariantRecord) -> None:
        """Comment includes DDG, domain, and 'computational' language."""
        result = format_clinvar_submission(pathogenic_record)

        assert result is not None
        comment = result["comment"]
        assert "3.9" in comment  # DDG value
        assert "BRCT" in comment  # Domain name
        assert "computational" in comment.lower()

    def test_format_includes_disclaimer(self, pathogenic_record: VariantRecord) -> None:
        """Comment has disclaimer and does NOT say 'assigned'."""
        result = format_clinvar_submission(pathogenic_record)

        assert result is not None
        comment = result["comment"]
        assert "DISCLAIMER" in comment
        assert "computational structural evidence" in comment
        assert "Russell Genetics" in comment
        assert "not replace it" in comment.lower() or "not replace" in comment.lower()
        # Must NOT claim evidence was "assigned" (it's "suggested")
        assert "assigned" not in comment.lower()


class TestGatingLogic:
    """Test that submissions are gated properly."""

    def test_incomplete_record_returns_none(self, incomplete_record: VariantRecord) -> None:
        """Record without M5 completed returns None."""
        result = format_clinvar_submission(incomplete_record)
        assert result is None

    def test_uncertain_classification_returns_none(self, pathogenic_record: VariantRecord) -> None:
        """Uncertain classification is not submitted."""
        pathogenic_record.classification = "uncertain"
        result = format_clinvar_submission(pathogenic_record)
        assert result is None

    def test_low_confidence_returns_none(self, pathogenic_record: VariantRecord) -> None:
        """Low model agreement is not submitted."""
        pathogenic_record.model_agreement = "low"
        result = format_clinvar_submission(pathogenic_record)
        assert result is None


class TestVariantCoordinates:
    """Test that variant coordinates are correctly included."""

    def test_format_contains_variant_coordinates(self, pathogenic_record: VariantRecord) -> None:
        """HGVS and gene symbol present in output."""
        result = format_clinvar_submission(pathogenic_record)

        assert result is not None
        assert "NM_007294.4:c.5095C>T" in result["hgvs_coding"]
        assert result["gene_symbol"] == "BRCA1"
        assert result["hgvs_protein"] == "p.Arg1699Trp"


class TestModelVersion:
    """Test that model version is included."""

    def test_format_contains_model_version(self, pathogenic_record: VariantRecord) -> None:
        """Ensemble version appears in comment or method_description."""
        result = format_clinvar_submission(pathogenic_record)

        assert result is not None
        version_in_comment = "v1.2" in result["comment"]
        version_in_method = "v1.2" in result["method_description"]
        assert version_in_comment or version_in_method


# =========================================================================
# TESTS — submit_to_clinvar
# =========================================================================


class TestSubmitDryRun:
    """Test dry run submission."""

    def test_submit_dry_run(self, pathogenic_record: VariantRecord) -> None:
        """Dry run returns status with payload."""
        submission = format_clinvar_submission(pathogenic_record)
        assert submission is not None

        result = submit_to_clinvar(submission, dry_run=True, use_test_endpoint=True)

        assert result["status"] == "dry_run"
        assert result["variant_id"] == "BRCA1_p.Arg1699Trp"
        assert result["payload"] == submission


class TestSubmitRequiresApiKey:
    """Test that live submission requires an API key."""

    def test_submit_requires_api_key(self, pathogenic_record: VariantRecord) -> None:
        """Without CLINVAR_API_KEY env var, returns error."""
        submission = format_clinvar_submission(pathogenic_record)
        assert submission is not None

        # Ensure env var is not set
        with patch.dict(os.environ, {}, clear=True):
            result = submit_to_clinvar(
                submission, dry_run=False, use_test_endpoint=True
            )

        assert result["status"] == "error"
        assert result["variant_id"] == "BRCA1_p.Arg1699Trp"
        assert "CLINVAR_API_KEY" in result["reason"]
