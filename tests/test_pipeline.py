"""Tests for the full pipeline orchestration."""

import pytest


class TestPipelineGracefulDegradation:
    def test_pipeline_never_crashes(self):
        """Pipeline should always return a VariantRecord, never raise."""
        from varis.pipeline import investigate
        record = investigate("FAKEGENE", "p.Ala1Val")
        assert record is not None
        assert record.gene_symbol == "FAKEGENE"

    @pytest.mark.timeout(120)
    def test_pipeline_tracks_failures(self):
        """Pipeline with bad gene should track failures but still complete."""
        from varis.pipeline import investigate
        record = investigate("FAKEGENE999", "p.Ala1Val")
        assert record is not None
        assert record.modules_completed is not None
        assert record.uniprot_id is None

    @pytest.mark.timeout(10)
    def test_pipeline_records_timing(self):
        """Pipeline should record processing time."""
        from varis.pipeline import investigate
        record = investigate("FAKEGENE", "p.Ala1Val")
        assert record.processing_time_seconds is not None
        assert record.processing_time_seconds >= 0
