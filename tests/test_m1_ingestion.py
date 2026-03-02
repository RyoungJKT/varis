"""Tests for M1: Data Ingestion module."""

class TestHGVSParser:
    def test_parse_three_letter(self):
        """Parse standard three-letter HGVS: p.Arg1699Trp"""
        pass
    def test_parse_extracts_position(self):
        pass
    def test_parse_charge_change(self):
        """Arg (+) to Trp (0) should give '+ve → neutral'."""
        pass
    def test_invalid_hgvs_returns_none(self):
        """Invalid notation should set fields to None, not crash."""
        pass

class TestClinVarClient:
    def test_fetch_known_variant(self):
        pass
    def test_unknown_variant_returns_none(self):
        pass

class TestGnomADClient:
    def test_fetch_rare_variant(self):
        pass
    def test_absent_variant_returns_none(self):
        pass
