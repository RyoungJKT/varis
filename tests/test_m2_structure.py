"""Tests for M2: Structure Retrieval and Preparation."""

import json
from unittest.mock import MagicMock

import pytest
from Bio.PDB import PDBParser

from varis.config import STRUCTURES_DIR
from varis.m2_structure.esmfold_predictor import predict_esmfold
from varis.m2_structure.structure_validator import validate_structure
from varis.models.variant_record import (
    VariantRecord,
    create_variant_record,
    RECORD_SCHEMA_VERSION,
)

BRCA1_PDB = STRUCTURES_DIR / "AF-P38398-F1-model_v6.pdb"


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


class TestStructureValidator:
    """Tests for structure_validator.validate_structure()."""

    def test_structure_sanity_fixture(self):
        """Parse the real BRCA1 PDB with BioPython and verify basic properties."""
        assert BRCA1_PDB.exists(), f"PDB fixture missing: {BRCA1_PDB}"
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("brca1", str(BRCA1_PDB))
        model = structure[0]
        chain = next(iter(model.get_chains()))
        standard_residues = [r for r in chain.get_residues() if r.id[0] == " "]
        assert len(standard_residues) >= 1699, (
            f"Expected >=1699 standard residues, got {len(standard_residues)}"
        )
        # Numbering starts at 1
        first_resid = standard_residues[0].id[1]
        assert first_resid == 1, f"Expected numbering to start at 1, got {first_resid}"

    def test_validate_existing_pdb(self, m1_completed_record):
        """Validate the real BRCA1 PDB — mutation site should be present."""
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.structure_source = "alphafold"
        result = validate_structure(m1_completed_record)
        assert result.mutation_site_present is True
        assert result.numbering_scheme == "uniprot_canonical"

    def test_validate_missing_residue(self, m1_completed_record):
        """Position 99999 cannot exist — mutation_site_present should be False."""
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.structure_source = "alphafold"
        m1_completed_record.residue_position = 99999
        result = validate_structure(m1_completed_record)
        assert result.mutation_site_present is False
        assert "mutation_site_present" in result.null_reasons

    def test_site_out_of_range(self, m1_completed_record):
        """Position 5000 is beyond BRCA1 length — should not be found."""
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.structure_source = "alphafold"
        m1_completed_record.residue_position = 5000
        result = validate_structure(m1_completed_record)
        assert result.mutation_site_present is False

    def test_plddt_only_for_predicted(self, m1_completed_record):
        """Experimental PDB structures should not have pLDDT extracted."""
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.structure_source = "pdb"
        result = validate_structure(m1_completed_record)
        assert result.plddt_available is False
        assert result.mutation_site_plddt is None

    def test_plddt_range_valid(self, m1_completed_record):
        """pLDDT values from AlphaFold should be in [0, 100]."""
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.structure_source = "alphafold"
        result = validate_structure(m1_completed_record)
        assert result.mutation_site_plddt is not None
        assert 0 <= result.mutation_site_plddt <= 100
        assert result.plddt_mean is not None
        assert 0 <= result.plddt_mean <= 100

    def test_plddt_site_bfactor_value(self, m1_completed_record):
        """Validator pLDDT should match the raw B-factor of residue 1699 CA."""
        # Get expected value directly from BioPython
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("brca1", str(BRCA1_PDB))
        model = structure[0]
        chain = next(iter(model.get_chains()))
        residue = chain[(" ", 1699, " ")]
        expected_bfactor = residue["CA"].get_bfactor()

        # Get value from validator
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.structure_source = "alphafold"
        result = validate_structure(m1_completed_record)

        assert result.mutation_site_plddt == pytest.approx(expected_bfactor, abs=0.01)

    def test_structure_quality_summary(self, m1_completed_record):
        """Quality summary should contain the required keys."""
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.structure_source = "alphafold"
        result = validate_structure(m1_completed_record)
        assert result.structure_quality_summary is not None
        summary = json.loads(result.structure_quality_summary)
        assert "plddt_mean" in summary
        assert "plddt_site" in summary
        assert "percent_low_confidence" in summary

    def test_pdb_hash_computed(self, m1_completed_record):
        """PDB hash should be a 64-character hex string (SHA-256)."""
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.structure_source = "alphafold"
        result = validate_structure(m1_completed_record)
        assert result.pdb_hash is not None
        assert len(result.pdb_hash) == 64
        # Verify it's valid hex
        int(result.pdb_hash, 16)


class TestESMFoldPredictor:
    """Tests for esmfold_predictor.predict_esmfold()."""

    def test_sequence_too_long(self):
        """Sequences over 400aa should be skipped with INTENTIONALLY_SKIPPED."""
        record = create_variant_record("TEST", "p.Met1Val")
        record.protein_sequence = "M" * 401
        record.pdb_path = None

        result = predict_esmfold(record)

        assert result.pdb_path is None
        assert "pdb_path" in result.null_reasons

    def test_no_sequence_available(self):
        """Missing protein_sequence should be skipped with UPSTREAM_DEPENDENCY_FAILED."""
        record = create_variant_record("TEST", "p.Met1Val")
        record.protein_sequence = None
        record.pdb_path = None

        result = predict_esmfold(record)

        assert result.pdb_path is None

    def test_skips_if_structure_exists(self):
        """If pdb_path is already set, ESMFold should not overwrite it."""
        record = create_variant_record("TEST", "p.Met1Val")
        record.pdb_path = "/some/existing.pdb"

        result = predict_esmfold(record)

        assert result.pdb_path == "/some/existing.pdb"

    def test_successful_prediction(self, tmp_path, monkeypatch):
        """Mock a successful ESMFold API call and verify PDB is saved."""
        # Set up the record
        record = create_variant_record("TEST", "p.Met1Val")
        record.protein_sequence = "MKFLILLFNILCLFPVLAADNHGVS"  # 25 aa
        record.pdb_path = None
        record.uniprot_id = "TEST123"

        # Redirect STRUCTURES_DIR to tmp_path
        monkeypatch.setattr(
            "varis.m2_structure.esmfold_predictor.STRUCTURES_DIR", tmp_path
        )

        # Mock the httpx client
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = (
            "ATOM      1  N   MET A   1"
            "       0.000   0.000   0.000\nEND\n"
        )
        mock_client = MagicMock()
        mock_client.post.return_value = mock_response

        result = predict_esmfold(record, client=mock_client)

        assert result.pdb_path is not None
        assert result.structure_source == "esmfold"


class TestPDBFixer:
    """Tests for pdb_fixer.fix_structure()."""

    def test_pdb_fixer_skipped_when_clean(self, m1_completed_record):
        """Clean AlphaFold PDB -> no fixer run, preparation_steps=['validated']."""
        from varis.m2_structure.pdb_fixer import fix_structure

        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.preparation_steps = ["validated"]
        result = fix_structure(m1_completed_record)
        assert result.pdb_fixed_path is None or result.pdb_fixed_path == result.pdb_path
        assert "validated" in result.preparation_steps

    def test_pdb_fixer_no_structure(self, m1_completed_record):
        """No pdb_path -> skip gracefully."""
        from varis.m2_structure.pdb_fixer import fix_structure

        m1_completed_record.pdb_path = None
        result = fix_structure(m1_completed_record)
        assert result.pdb_fixed_path is None

    def test_pdb_fixer_output_parses(self, m1_completed_record):
        """If fixer runs, output PDB still parses successfully."""
        from varis.m2_structure.pdb_fixer import fix_structure

        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.preparation_steps = ["validated"]
        result = fix_structure(m1_completed_record)
        pdb_to_check = result.pdb_fixed_path or result.pdb_path
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("test", pdb_to_check)
        assert structure is not None

    def test_site_preserved_after_fix(self, m1_completed_record):
        """Target residue still present after PDBFixer."""
        from varis.m2_structure.pdb_fixer import fix_structure

        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.preparation_steps = ["validated"]
        result = fix_structure(m1_completed_record)
        pdb_to_check = result.pdb_fixed_path or result.pdb_path
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("test", pdb_to_check)
        chain = next(iter(structure[0].get_chains()))
        residue_ids = {r.id[1] for r in chain.get_residues() if r.id[0] == " "}
        assert 1699 in residue_ids


class TestM2Integration:
    """Integration tests for the full M2 orchestrator."""

    def test_m2_full_pipeline_brca1(self, m1_completed_record):
        """Full M2 on BRCA1 AlphaFold structure."""
        from varis.m2_structure import run

        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.structure_source = "alphafold"
        result = run(m1_completed_record)
        assert "M2" in result.modules_completed
        assert result.mutation_site_present is True
        assert result.pdb_hash is not None
        assert result.numbering_scheme == "uniprot_canonical"
        assert result.plddt_available is True
        assert result.mutation_site_plddt is not None
        assert 0.0 <= result.mutation_site_plddt <= 100.0

    def test_m2_no_structure(self, m1_completed_record):
        """M2 with no structure -> marks failed."""
        from varis.m2_structure import run

        m1_completed_record.pdb_path = None
        m1_completed_record.structure_source = None
        m1_completed_record.protein_sequence = "M" * 500  # Too long for ESMFold
        result = run(m1_completed_record)
        assert "M2" in result.modules_failed

    def test_m2_golden_record_keys(self, m1_completed_record):
        """Verify all expected keys are set after M2."""
        from varis.m2_structure import run

        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.structure_source = "alphafold"
        result = run(m1_completed_record)
        expected_keys = [
            "mutation_site_present", "pdb_hash", "numbering_scheme",
            "preparation_steps", "plddt_available",
        ]
        for key in expected_keys:
            assert getattr(result, key) is not None, f"{key} should not be None"
