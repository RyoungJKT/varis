"""Tests for M5: ML Scoring Engine."""

class TestFeatureExtractor:
    def test_extract_from_full_record(self, fully_populated_record):
        from varis.m5_scoring.feature_extractor import extract_features
        features = extract_features(fully_populated_record)
        assert len(features) == 15

    def test_extract_from_partial_record(self, m1_completed_record):
        from varis.m5_scoring.feature_extractor import extract_features
        features = extract_features(m1_completed_record)
        assert features["ddg_foldx"] is None
        assert features["alphamissense_score"] is not None

class TestACMGMapper:
    def test_ps3_proxy_with_high_ddg(self):
        pass
    def test_pm1_in_critical_domain(self):
        pass
