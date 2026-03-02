"""Tests for the full pipeline orchestration."""

class TestPipelineGracefulDegradation:
    def test_pipeline_never_crashes(self):
        """Pipeline should always return a VariantRecord, never raise."""
        from varis.pipeline import investigate
        record = investigate("FAKEGENE", "p.Ala1Val")
        assert record is not None
        assert record.gene_symbol == "FAKEGENE"

    def test_pipeline_tracks_failures(self):
        pass

    def test_pipeline_records_timing(self):
        pass
