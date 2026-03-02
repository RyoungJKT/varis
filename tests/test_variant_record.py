"""Tests for the VariantRecord data model."""
from varis.models.variant_record import VariantRecord, create_variant_record

class TestVariantRecordCreation:
    def test_create_sets_input_fields(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert record.gene_symbol == "BRCA1"
        assert record.hgvs_protein == "p.Arg1699Trp"
        assert record.variant_id == "BRCA1_p.Arg1699Trp"

    def test_create_initializes_tracking(self):
        record = create_variant_record("TP53", "p.Arg175His")
        assert record.modules_completed == []
        assert record.modules_failed == []
        assert record.investigation_timestamp is not None

    def test_all_fields_default_to_none(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert record.clinvar_id is None
        assert record.ddg_foldx is None
        assert record.score_ensemble is None

class TestSerialization:
    def test_to_json_and_back(self, empty_record):
        json_str = empty_record.to_json()
        restored = VariantRecord.from_json(json_str)
        assert restored.gene_symbol == empty_record.gene_symbol

    def test_to_dict(self, empty_record):
        d = empty_record.to_dict()
        assert isinstance(d, dict)
        assert d["gene_symbol"] == "BRCA1"

    def test_from_dict_ignores_unknown_keys(self):
        data = {"gene_symbol": "BRCA1", "unknown_field": "ignored"}
        record = VariantRecord.from_dict(data)
        assert record.gene_symbol == "BRCA1"

    def test_save_and_load(self, empty_record, tmp_path):
        path = str(tmp_path / "test.json")
        empty_record.save(path)
        loaded = VariantRecord.load(path)
        assert loaded.gene_symbol == empty_record.gene_symbol

class TestFeatureExtraction:
    def test_get_structural_features(self, fully_populated_record):
        features = fully_populated_record.get_structural_features()
        assert features["ddg_foldx"] == 3.7
        assert features["conservation_score"] == 1.0

    def test_count_available_features_full(self, fully_populated_record):
        assert fully_populated_record.count_available_features() == 15

    def test_count_available_features_empty(self, empty_record):
        assert empty_record.count_available_features() == 0

class TestModuleTracking:
    def test_mark_module_completed(self, empty_record):
        empty_record.mark_module_completed("M1")
        assert "M1" in empty_record.modules_completed

    def test_mark_module_failed(self, empty_record):
        empty_record.mark_module_failed("M3.foldx")
        assert "M3.foldx" in empty_record.modules_failed

    def test_no_duplicate_tracking(self, empty_record):
        empty_record.mark_module_completed("M1")
        empty_record.mark_module_completed("M1")
        assert empty_record.modules_completed.count("M1") == 1
