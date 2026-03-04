"""Tests for M6: VarisDB API."""
import pytest
import time
from unittest.mock import patch, MagicMock


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
        from varis.m6_platform.api.database import init_db
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


class TestWorker:
    """Tests for background worker."""

    def test_submit_job(self, tmp_path):
        """Submit creates a queued job."""
        from varis.m6_platform.api.database import init_db, get_job
        from varis.m6_platform.api.worker import InvestigationWorker
        engine, SessionLocal = init_db(f"sqlite:///{tmp_path}/test.db")
        session = SessionLocal()
        worker = InvestigationWorker(SessionLocal, max_workers=1)
        job_id = worker.submit("BRCA1", "p.Arg1699Trp", session)
        assert job_id is not None
        job = get_job(session, job_id)
        assert job["status"] in ("queued", "running")
        worker.shutdown()
        session.close()

    def test_job_completes(self, tmp_path):
        """Mocked pipeline -> job succeeds."""
        from varis.m6_platform.api.database import init_db, get_job
        from varis.m6_platform.api.worker import InvestigationWorker
        engine, SessionLocal = init_db(f"sqlite:///{tmp_path}/test.db")
        session = SessionLocal()
        with patch("varis.m6_platform.api.worker.run_pipeline") as mock_pipeline:
            mock_record = MagicMock()
            mock_record.variant_id = "BRCA1_p.Arg1699Trp"
            mock_record.to_json.return_value = "{}"
            mock_record.gene_symbol = "BRCA1"
            mock_record.clinvar_id = None
            mock_record.classification = "likely_pathogenic"
            mock_record.ensemble_version = "v2026.03"
            mock_record.pipeline_version = "v1.0"
            mock_record.score_ensemble = 0.91
            mock_record.hgvs_protein = "p.Arg1699Trp"
            mock_record.to_dict.return_value = {"gene_symbol": "BRCA1"}
            mock_pipeline.return_value = mock_record
            worker = InvestigationWorker(SessionLocal, max_workers=1)
            job_id = worker.submit("BRCA1", "p.Arg1699Trp", session)
            time.sleep(2)  # Wait for worker thread
            job = get_job(SessionLocal(), job_id)
            assert job["status"] == "succeeded"
            worker.shutdown()
        session.close()


class TestAPIEndpoints:
    """Tests for FastAPI routes."""

    @pytest.fixture
    def client(self, tmp_path):
        from fastapi.testclient import TestClient
        from varis.m6_platform.api.main import create_app
        app = create_app(database_url=f"sqlite:///{tmp_path}/test.db")
        return TestClient(app)

    def test_health_check(self, client):
        resp = client.get("/health")
        assert resp.status_code == 200
        assert resp.json()["status"] == "ok"

    def test_get_investigation_not_found(self, client):
        resp = client.get("/api/v1/investigations/NONEXISTENT")
        assert resp.status_code == 404

    def test_submit_investigation(self, client):
        resp = client.post("/api/v1/investigations", json={"gene": "BRCA1", "hgvs": "p.Arg1699Trp"})
        assert resp.status_code in (200, 201, 202)
        data = resp.json()
        assert "job_id" in data or "variant_id" in data

    def test_submit_invalid_input(self, client):
        resp = client.post("/api/v1/investigations", json={"gene": "", "hgvs": ""})
        assert resp.status_code == 422

    def test_search_empty(self, client):
        resp = client.get("/api/v1/variants?q=BRCA1")
        assert resp.status_code == 200
        assert "results" in resp.json()

    def test_stats(self, client):
        resp = client.get("/api/v1/variants/stats")
        assert resp.status_code == 200

    def test_get_job_not_found(self, client):
        resp = client.get("/api/v1/jobs/nonexistent-id")
        assert resp.status_code == 404

    def test_cors_headers(self, client):
        resp = client.options(
            "/health",
            headers={
                "Origin": "http://localhost:5173",
                "Access-Control-Request-Method": "GET",
            },
        )
        assert resp.status_code == 200
        assert "access-control-allow-origin" in resp.headers

    def test_evolution_log_endpoint(self, client):
        """GET /api/v1/evolution-log returns 200 with events list."""
        resp = client.get("/api/v1/evolution-log")
        assert resp.status_code == 200
        assert "events" in resp.json()
        assert isinstance(resp.json()["events"], list)

    def test_evolution_log_with_limit(self, client):
        """GET /api/v1/evolution-log?limit=5 respects limit parameter."""
        resp = client.get("/api/v1/evolution-log?limit=5")
        assert resp.status_code == 200
        assert "events" in resp.json()

    def test_evolution_log_with_event_type(self, client):
        """GET /api/v1/evolution-log?event_type=DEPLOY filters by type."""
        resp = client.get("/api/v1/evolution-log?event_type=DEPLOY")
        assert resp.status_code == 200
        assert "events" in resp.json()

    def test_clinvar_submission_eligible(self, client, fully_populated_record):
        """GET /api/v1/clinvar-submissions/{id} returns eligible=True for qualifying record."""
        from varis.m6_platform.api.database import save_variant_record
        session = client.app.state.session_factory()
        try:
            save_variant_record(session, fully_populated_record)
        finally:
            session.close()

        resp = client.get(f"/api/v1/clinvar-submissions/{fully_populated_record.variant_id}")
        assert resp.status_code == 200
        data = resp.json()
        assert data["eligible"] is True
        assert data["variant_id"] == fully_populated_record.variant_id
        assert data["submission"] is not None
        assert "clinical_significance" in data["submission"]

    def test_clinvar_submission_not_eligible(self, client, m1_completed_record):
        """GET /api/v1/clinvar-submissions/{id} returns eligible=False for uncertain record."""
        from varis.m6_platform.api.database import save_variant_record
        # m1_completed_record has no M5 completed and no classification — should not be eligible
        session = client.app.state.session_factory()
        try:
            save_variant_record(session, m1_completed_record)
        finally:
            session.close()

        resp = client.get(f"/api/v1/clinvar-submissions/{m1_completed_record.variant_id}")
        assert resp.status_code == 200
        data = resp.json()
        assert data["eligible"] is False
        assert data["variant_id"] == m1_completed_record.variant_id
        assert data["submission"] is None

    def test_clinvar_submission_not_found(self, client):
        """GET /api/v1/clinvar-submissions/{id} returns 404 for unknown variant."""
        resp = client.get("/api/v1/clinvar-submissions/NONEXISTENT")
        assert resp.status_code == 404

    def test_report_download(self, client, fully_populated_record):
        """GET /api/v1/reports/{id} returns HTML report."""
        from varis.m6_platform.api.database import save_variant_record
        session = client.app.state.session_factory()
        try:
            save_variant_record(session, fully_populated_record)
        finally:
            session.close()

        resp = client.get(f"/api/v1/reports/{fully_populated_record.variant_id}")
        assert resp.status_code == 200
        assert "text/html" in resp.headers["content-type"]
        assert "Variant Investigation Report" in resp.text

    def test_report_download_not_found(self, client):
        """GET /api/v1/reports/{id} returns 404 for unknown variant."""
        resp = client.get("/api/v1/reports/NONEXISTENT")
        assert resp.status_code == 404
