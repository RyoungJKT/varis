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
