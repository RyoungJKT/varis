"""Tests for M6: VarisDB API."""
import pytest
import time
from unittest.mock import patch, MagicMock
from pydantic import ValidationError


class TestPydanticModels:
    """Tests for API response models."""

    def test_investigation_response_valid(self):
        from varis.m6_platform.api.models import InvestigationResponse
        data = {
            "variant_id": "BRCA1_p.Arg1699Trp",
            "gene": "BRCA1",
            "hgvs": "p.Arg1699Trp",
            "structure": {"source": "alphafold"},
            "features": [],
            "prediction": {"score": 0.91, "classification": "likely_pathogenic"},
            "explanation": [],
            "provenance": {
                "data_sources": ["ClinVar"],
                "modules_completed": ["M1"],
            },
        }
        resp = InvestigationResponse(**data)
        assert resp.variant_id == "BRCA1_p.Arg1699Trp"

    def test_feature_item_model(self):
        from varis.m6_platform.api.models import FeatureItem
        item = FeatureItem(name="conservation_score", value=0.95, units=None, evidence_tag=None, available=True)
        assert item.available is True

    def test_explanation_item_model(self):
        from varis.m6_platform.api.models import ExplanationItem
        item = ExplanationItem(feature="conservation_score", value=0.95, shap=0.24)
        assert item.shap == 0.24

    def test_job_status_model(self):
        from varis.m6_platform.api.models import JobStatusResponse
        job = JobStatusResponse(job_id="abc-123", status="running", variant_id="BRCA1_p.Arg1699Trp", current_step="M3: Structural Analysis")
        assert job.status == "running"


class TestDatabase:
    """Tests for database.py — SQLAlchemy ORM."""

    @pytest.fixture
    def db_session(self, tmp_path):
        """Create an in-memory SQLite session for testing."""
        from varis.m6_platform.api.database import init_db, get_session
        engine, SessionLocal = init_db(f"sqlite:///{tmp_path}/test.db")
        session = SessionLocal()
        yield session
        session.close()

    def test_save_and_get_variant(self, db_session, fully_populated_record):
        """Save a VariantRecord and retrieve it."""
        from varis.m6_platform.api.database import save_variant_record, get_variant_record
        variant_id = save_variant_record(db_session, fully_populated_record)
        assert variant_id is not None
        retrieved = get_variant_record(db_session, variant_id)
        assert retrieved is not None
        assert retrieved["gene_symbol"] == "BRCA1"

    def test_get_nonexistent_variant(self, db_session):
        """Unknown variant_id returns None."""
        from varis.m6_platform.api.database import get_variant_record
        result = get_variant_record(db_session, "NONEXISTENT")
        assert result is None

    def test_search_by_gene(self, db_session, fully_populated_record):
        """Search finds variant by gene name."""
        from varis.m6_platform.api.database import save_variant_record, search_variants
        save_variant_record(db_session, fully_populated_record)
        results = search_variants(db_session, "BRCA1")
        assert len(results) >= 1

    def test_create_and_get_job(self, db_session):
        """Create a job and retrieve its status."""
        from varis.m6_platform.api.database import create_job, get_job, update_job_status
        job_id = create_job(db_session, "BRCA1_p.Arg1699Trp")
        assert job_id is not None
        job = get_job(db_session, job_id)
        assert job["status"] == "queued"
        update_job_status(db_session, job_id, "running", current_step="M1: Ingestion")
        job = get_job(db_session, job_id)
        assert job["status"] == "running"
        assert job["current_step"] == "M1: Ingestion"

    def test_idempotent_variant_exists(self, db_session, fully_populated_record):
        """Check if variant already exists."""
        from varis.m6_platform.api.database import save_variant_record, variant_exists
        save_variant_record(db_session, fully_populated_record)
        assert variant_exists(db_session, fully_populated_record.variant_id) is True
        assert variant_exists(db_session, "NONEXISTENT") is False


class TestValidation:
    """Tests for input validation."""

    def test_valid_hgvs(self):
        from varis.m6_platform.api.validation import validate_variant_input
        result = validate_variant_input("BRCA1", "p.Arg1699Trp")
        assert result["valid"] is True

    def test_security_length_cap(self):
        from varis.m6_platform.api.validation import validate_variant_input
        result = validate_variant_input("A" * 300, "p.Arg1699Trp")
        assert result["valid"] is False
        assert "length" in result["error"].lower()

    def test_security_bad_chars(self):
        from varis.m6_platform.api.validation import validate_variant_input
        result = validate_variant_input("BRCA1; DROP TABLE", "p.Arg1699Trp")
        assert result["valid"] is False

    def test_friendly_error_message(self):
        from varis.m6_platform.api.validation import validate_variant_input
        result = validate_variant_input("BRCA1", "invalid_format")
        assert result["valid"] is False
        assert "Expected" in result["error"]


class TestInvestigationBuilder:
    """Tests for building InvestigationResponse from VariantRecord."""

    def test_build_from_full_record(self, fully_populated_record):
        """Full record produces complete response."""
        from varis.m6_platform.api.investigation import build_investigation_response
        resp = build_investigation_response(fully_populated_record)
        assert resp.variant_id is not None
        assert resp.structure.source == "alphafold"
        assert resp.prediction.score is not None
        assert len(resp.provenance.data_sources) > 0
        assert len(resp.provenance.modules_completed) > 0

    def test_build_features_list(self, fully_populated_record):
        """Features list includes availability flags."""
        from varis.m6_platform.api.investigation import build_investigation_response
        resp = build_investigation_response(fully_populated_record)
        assert len(resp.features) > 0
        for f in resp.features:
            assert hasattr(f, "available")

    def test_build_from_partial_record(self, m1_completed_record):
        """Partial record (M1 only) still produces valid response."""
        from varis.m6_platform.api.investigation import build_investigation_response
        resp = build_investigation_response(m1_completed_record)
        assert resp.variant_id is not None
        assert resp.prediction.score is None  # M5 not run yet
        assert len(resp.explanation) == 0  # No SHAP without M5

    def test_explanation_pre_sorted(self, fully_populated_record):
        """SHAP items come pre-sorted by |shap| descending."""
        from varis.m6_platform.api.investigation import build_investigation_response
        fully_populated_record.shap_top_features = [
            {"feature": "a", "value": 1.0, "shap": 0.1},
            {"feature": "b", "value": 2.0, "shap": -0.3},
            {"feature": "c", "value": 3.0, "shap": 0.2},
        ]
        resp = build_investigation_response(fully_populated_record)
        shap_abs = [abs(e.shap) for e in resp.explanation]
        assert shap_abs == sorted(shap_abs, reverse=True)
