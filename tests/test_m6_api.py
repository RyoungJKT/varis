"""Tests for M6: VarisDB API."""
import pytest
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
