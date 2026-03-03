"""Tests for M4: Conservation Analysis."""
import pytest
import math
from unittest.mock import MagicMock, patch
from varis.models.variant_record import (
    VariantRecord, create_variant_record, NullReason, RECORD_SCHEMA_VERSION,
)


class TestSchemaV130:
    """Verify schema v1.3.0 fields exist on VariantRecord."""

    def test_schema_version_is_1_3_0(self):
        assert RECORD_SCHEMA_VERSION == "1.3.0"

    def test_new_msa_fields_exist(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert hasattr(record, "msa_num_sequences")
        assert hasattr(record, "msa_gap_fraction_at_site")
        assert hasattr(record, "msa_column_index")

    def test_insufficient_data_reason_exists(self):
        assert hasattr(NullReason, "INSUFFICIENT_DATA")
        assert NullReason.INSUFFICIENT_DATA == "insufficient_data"


class TestConservationScorer:
    """Tests for conservation_scorer.py — entropy and position mapping."""

    def test_entropy_fully_conserved(self):
        """All same AA -> entropy=0, score=1.0."""
        from varis.m4_conservation.conservation_scorer import _shannon_entropy
        column = ["R"] * 20
        assert _shannon_entropy(column) == pytest.approx(0.0, abs=0.001)

    def test_entropy_maximally_variable(self):
        """All different AAs -> entropy ~ log2(20), score ~ 0.0."""
        from varis.m4_conservation.conservation_scorer import _shannon_entropy
        aas = list("ACDEFGHIKLMNPQRSTVWY")
        assert _shannon_entropy(aas) == pytest.approx(math.log2(20), abs=0.01)

    def test_entropy_ignores_gaps(self):
        """Gaps excluded from frequency, but present in input."""
        from varis.m4_conservation.conservation_scorer import _shannon_entropy
        column = ["R"] * 18 + ["-", "-"]
        assert _shannon_entropy(column) == pytest.approx(0.0, abs=0.001)

    def test_position_mapping_no_gaps(self):
        """Direct mapping when query has no gaps."""
        from varis.m4_conservation.conservation_scorer import _map_position_to_column
        query_row = "MKRST"
        # Position is 1-indexed: M=pos1(col0), K=pos2(col1), R=pos3(col2)
        col = _map_position_to_column(query_row, position=3)
        assert col == 2

    def test_position_mapping_with_gaps(self):
        """Gaps in query shift the column index."""
        from varis.m4_conservation.conservation_scorer import _map_position_to_column
        query_row = "M-KR-ST"
        # M=pos1(col0), K=pos2(col2), R=pos3(col3)
        col = _map_position_to_column(query_row, position=3)
        assert col == 3

    def test_position_mapping_validates_ref_aa(self):
        """Returns None if ref AA doesn't match."""
        from varis.m4_conservation.conservation_scorer import _map_position_to_column
        query_row = "MKRST"
        col = _map_position_to_column(query_row, position=3, expected_aa="W")
        assert col is None  # R != W

    def test_score_conservation_full(self, m1_completed_record):
        """Full scorer with a small test alignment — all conserved."""
        from varis.m4_conservation.conservation_scorer import score_conservation
        alignment = {
            "sequences": {
                "query": "MKRST",
                **{f"orth{i}": "MKRST" for i in range(1, 11)},
            },
            "query_id": "query",
            "taxonomy": {
                "orth1": 9606, "orth2": 9615, "orth3": 10090,
                "orth4": 9913, "orth5": 9823,
                "orth6": 7955, "orth7": 8364, "orth8": 9031,
                "orth9": 28377, "orth10": 13616,
            },
        }
        m1_completed_record.residue_position = 3
        m1_completed_record.ref_aa_single = "R"
        result = score_conservation(m1_completed_record, alignment)
        assert result.conservation_score == pytest.approx(1.0, abs=0.01)
        assert result.position_entropy == pytest.approx(0.0, abs=0.01)
        assert result.msa_column_index == 2
        assert result.conservation_available is True

    def test_mammal_conservation_threshold(self, m1_completed_record):
        """>=90% mammals have ref AA -> True."""
        from varis.m4_conservation.conservation_scorer import score_conservation
        alignment = {
            "sequences": {
                "query": "MKRST",
                "m1": "MKRST", "m2": "MKRST", "m3": "MKRST",
                "m4": "MKRST", "m5": "MKRST",
                "o1": "MKWST",
            },
            "query_id": "query",
            "taxonomy": {
                "m1": 9606, "m2": 9615, "m3": 10090, "m4": 9913, "m5": 9823,
                "o1": 7955,
            },
        }
        m1_completed_record.residue_position = 3
        m1_completed_record.ref_aa_single = "R"
        result = score_conservation(m1_completed_record, alignment)
        assert result.conserved_across_mammals is True

    def test_mammal_conservation_too_few(self, m1_completed_record):
        """<5 mammals -> conserved_across_mammals=None."""
        from varis.m4_conservation.conservation_scorer import score_conservation
        alignment = {
            "sequences": {
                "query": "MKRST",
                "m1": "MKRST", "m2": "MKRST",
                **{f"o{i}": "MKRST" for i in range(1, 9)},
            },
            "query_id": "query",
            "taxonomy": {"m1": 9606, "m2": 9615,
                          **{f"o{i}": 7955 + i for i in range(1, 9)}},
        }
        m1_completed_record.residue_position = 3
        m1_completed_record.ref_aa_single = "R"
        result = score_conservation(m1_completed_record, alignment)
        assert result.conserved_across_mammals is None
