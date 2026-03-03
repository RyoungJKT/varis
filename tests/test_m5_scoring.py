"""Tests for M5: ML Scoring Engine."""
import pytest
import numpy as np
from unittest.mock import MagicMock, patch
from varis.models.variant_record import (
    VariantRecord, create_variant_record, RECORD_SCHEMA_VERSION,
)


class TestSchemaV140:
    """Verify schema v1.4.0 — evidence tags replace ACMG codes."""

    def test_schema_version(self):
        assert RECORD_SCHEMA_VERSION == "1.4.0"

    def test_evidence_tag_fields_exist(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert hasattr(record, "evidence_tags")
        assert hasattr(record, "evidence_computational_support")
        assert hasattr(record, "evidence_rarity")
        assert hasattr(record, "evidence_energetics")
        assert hasattr(record, "evidence_domain_context")

    def test_old_acmg_fields_removed(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert not hasattr(record, "acmg_codes")
        assert not hasattr(record, "acmg_pp3")
        assert not hasattr(record, "acmg_pp2")
        assert not hasattr(record, "acmg_ps3_proxy")
        assert not hasattr(record, "acmg_pm1")


class TestFeatureExtractor:
    """Tests for feature_extractor.py — one-hot, strict alignment."""

    def test_extract_from_full_record(self, fully_populated_record):
        """Full record extracts all features as ordered dict."""
        from varis.m5_scoring.feature_extractor import extract_features
        features = extract_features(fully_populated_record)
        assert isinstance(features, dict)
        assert len(features) > 15  # Expanded by one-hot

    def test_extract_from_partial_record(self, m1_completed_record):
        """Missing features are None, not absent."""
        from varis.m5_scoring.feature_extractor import extract_features
        features = extract_features(m1_completed_record)
        assert features["ddg_foldx"] is None
        assert "ddg_available" in features  # Key present (value is None for M1-only)

    def test_one_hot_encoding(self, fully_populated_record):
        """Categoricals expanded to binary columns."""
        from varis.m5_scoring.feature_extractor import extract_features
        features = extract_features(fully_populated_record)
        # secondary_structure_name="helix" → one-hot columns
        assert "secondary_structure_helix" in features
        assert "secondary_structure_sheet" in features
        assert "secondary_structure_coil" in features
        assert features["secondary_structure_helix"] == 1
        assert features["secondary_structure_sheet"] == 0

    def test_in_domain_flag(self, fully_populated_record):
        """domain_name converted to binary in_domain flag."""
        from varis.m5_scoring.feature_extractor import extract_features
        features = extract_features(fully_populated_record)
        assert "in_domain" in features
        assert features["in_domain"] == 1  # domain_name is "BRCT"
        assert "domain_name" not in features  # Raw name removed

    def test_feature_extraction_allows_null(self, m1_completed_record):
        """None values pass through (real missingness)."""
        from varis.m5_scoring.feature_extractor import extract_features
        features = extract_features(m1_completed_record)
        # Structural features should be None
        assert features["solvent_accessibility_relative"] is None

    def test_extract_with_column_order(self, fully_populated_record):
        """Features extracted in deterministic order."""
        from varis.m5_scoring.feature_extractor import extract_features
        f1 = extract_features(fully_populated_record)
        f2 = extract_features(fully_populated_record)
        assert list(f1.keys()) == list(f2.keys())

    def test_build_feature_vector(self, fully_populated_record):
        """Converts feature dict to ordered numpy array for model input."""
        from varis.m5_scoring.feature_extractor import extract_features, build_feature_vector
        features = extract_features(fully_populated_record)
        columns = list(features.keys())
        vector = build_feature_vector(features, columns)
        assert len(vector) == len(columns)

    def test_build_feature_vector_fails_on_missing_column(self, fully_populated_record):
        """Raises ValueError if a required column is missing."""
        from varis.m5_scoring.feature_extractor import extract_features, build_feature_vector
        features = extract_features(fully_populated_record)
        bad_columns = list(features.keys()) + ["nonexistent_column"]
        with pytest.raises(ValueError, match="Missing feature"):
            build_feature_vector(features, bad_columns)

    def test_missingness_simulation(self, fully_populated_record):
        """Drops feature blocks, sets availability flags."""
        from varis.m5_scoring.feature_extractor import simulate_missingness
        import copy
        record = copy.deepcopy(fully_populated_record)
        # Force a specific group to be dropped
        modified = simulate_missingness(record, groups=["sasa"])
        features = modified.get_ml_features()
        assert features["sasa_available"] is False
        assert modified.solvent_accessibility_relative is None
