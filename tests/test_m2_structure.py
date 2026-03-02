"""Tests for M2: Structure Retrieval and Preparation."""

from varis.models.variant_record import (
    VariantRecord,
    create_variant_record,
    RECORD_SCHEMA_VERSION,
)


class TestSchemaV120:
    """Verify VariantRecord schema v1.2.0 field additions and removals."""

    def test_schema_version(self):
        """Schema version should be 1.2.0."""
        assert RECORD_SCHEMA_VERSION == "1.2.0"
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert record.record_schema_version == "1.2.0"

    # --- New M2 fields ---

    def test_m2_field_mutation_site_present(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert hasattr(record, "mutation_site_present")
        assert record.mutation_site_present is None

    def test_m2_field_mutation_site_plddt(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert hasattr(record, "mutation_site_plddt")
        assert record.mutation_site_plddt is None

    def test_m2_field_plddt_available(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert hasattr(record, "plddt_available")
        assert record.plddt_available is None

    def test_m2_field_mutation_site_confidence_bucket(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert hasattr(record, "mutation_site_confidence_bucket")
        assert record.mutation_site_confidence_bucket is None

    def test_m2_field_numbering_scheme(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert hasattr(record, "numbering_scheme")
        assert record.numbering_scheme is None

    def test_m2_field_structure_quality_summary(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert hasattr(record, "structure_quality_summary")
        assert record.structure_quality_summary is None

    def test_m2_field_preparation_steps(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert hasattr(record, "preparation_steps")
        assert record.preparation_steps is None

    def test_m2_field_pdb_hash(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert hasattr(record, "pdb_hash")
        assert record.pdb_hash is None

    def test_m2_field_structure_source_url(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert hasattr(record, "structure_source_url")
        assert record.structure_source_url is None

    # --- New M3 fields ---

    def test_m3_field_solvent_accessibility_relative(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert hasattr(record, "solvent_accessibility_relative")
        assert record.solvent_accessibility_relative is None

    def test_m3_field_contacts_wt(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert hasattr(record, "contacts_wt")
        assert record.contacts_wt is None

    def test_m3_field_hbonds_wt(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert hasattr(record, "hbonds_wt")
        assert record.hbonds_wt is None

    def test_m3_field_packing_density(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert hasattr(record, "packing_density")
        assert record.packing_density is None

    def test_m3_field_domain_start(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert hasattr(record, "domain_start")
        assert record.domain_start is None

    def test_m3_field_domain_end(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert hasattr(record, "domain_end")
        assert record.domain_end is None

    # --- Removed fields must NOT exist ---

    def test_removed_helix_disruption(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert not hasattr(record, "helix_disruption")

    def test_removed_hbonds_lost(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert not hasattr(record, "hbonds_lost")

    def test_removed_contacts_changed(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert not hasattr(record, "contacts_changed")

    def test_removed_structure_resolution(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert not hasattr(record, "structure_resolution")

    def test_removed_solvent_accessibility(self):
        """Old field name (absolute SASA) should be gone; replaced by _relative."""
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert not hasattr(record, "solvent_accessibility")

    def test_removed_functional_site_distance(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert not hasattr(record, "functional_site_distance")

    def test_removed_nearest_functional_site(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert not hasattr(record, "nearest_functional_site")

    def test_removed_plddt_score(self):
        """plddt_score is replaced by mutation_site_plddt."""
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert not hasattr(record, "plddt_score")

    # --- Serialization round-trip for new fields ---

    def test_new_fields_serialize(self):
        """New M2/M3 fields should survive to_dict/from_dict round-trip."""
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record.mutation_site_present = True
        record.mutation_site_plddt = 92.4
        record.plddt_available = True
        record.mutation_site_confidence_bucket = "high"
        record.numbering_scheme = "uniprot"
        record.structure_quality_summary = "good"
        record.preparation_steps = ["download", "fix_residues"]
        record.pdb_hash = "abc123"
        record.structure_source_url = "https://alphafold.ebi.ac.uk/files/AF-P38398-F1-model_v4.pdb"
        record.solvent_accessibility_relative = 0.03
        record.contacts_wt = 12
        record.hbonds_wt = 3
        record.packing_density = 0.72
        record.domain_start = 1646
        record.domain_end = 1736

        d = record.to_dict()
        restored = VariantRecord.from_dict(d)
        assert restored.mutation_site_present is True
        assert restored.mutation_site_plddt == 92.4
        assert restored.plddt_available is True
        assert restored.mutation_site_confidence_bucket == "high"
        assert restored.numbering_scheme == "uniprot"
        assert restored.structure_quality_summary == "good"
        assert restored.preparation_steps == ["download", "fix_residues"]
        assert restored.pdb_hash == "abc123"
        assert restored.structure_source_url == "https://alphafold.ebi.ac.uk/files/AF-P38398-F1-model_v4.pdb"
        assert restored.solvent_accessibility_relative == 0.03
        assert restored.contacts_wt == 12
        assert restored.hbonds_wt == 3
        assert restored.packing_density == 0.72
        assert restored.domain_start == 1646
        assert restored.domain_end == 1736
