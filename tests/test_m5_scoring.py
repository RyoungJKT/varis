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
