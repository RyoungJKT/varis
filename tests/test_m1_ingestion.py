"""Tests for M1: Data Ingestion module."""

import pytest
from varis.models.variant_record import create_variant_record


class TestHGVSParser:
    def test_parse_three_letter(self):
        """Parse standard three-letter HGVS: p.Arg1699Trp"""
        from varis.m1_ingestion.hgvs_parser import parse_hgvs
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record = parse_hgvs(record)
        assert record.ref_amino_acid == "Arg"
        assert record.alt_amino_acid == "Trp"
        assert record.ref_aa_single == "R"
        assert record.alt_aa_single == "W"

    def test_parse_extracts_position(self):
        """Position should be extracted as an integer."""
        from varis.m1_ingestion.hgvs_parser import parse_hgvs
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record = parse_hgvs(record)
        assert record.residue_position == 1699

    def test_parse_single_letter(self):
        """Also handle single-letter notation: p.R1699W"""
        from varis.m1_ingestion.hgvs_parser import parse_hgvs
        record = create_variant_record("BRCA1", "p.R1699W")
        record = parse_hgvs(record)
        assert record.residue_position == 1699
        assert record.ref_amino_acid == "Arg"
        assert record.alt_amino_acid == "Trp"
        assert record.ref_aa_single == "R"
        assert record.alt_aa_single == "W"

    def test_parse_charge_change(self):
        """Arg (+) to Trp (0) should produce a charge change string."""
        from varis.m1_ingestion.hgvs_parser import parse_hgvs
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record = parse_hgvs(record)
        assert record.charge_change is not None
        assert "positive" in record.charge_change.lower()

    def test_parse_no_charge_change(self):
        """Ala (0) to Val (0) should indicate no charge change."""
        from varis.m1_ingestion.hgvs_parser import parse_hgvs
        record = create_variant_record("FAKEGENE", "p.Ala100Val")
        record = parse_hgvs(record)
        assert record.charge_change is not None
        assert "no change" in record.charge_change.lower()

    def test_invalid_hgvs_returns_none(self):
        """Invalid notation should set fields to None with reason, not crash."""
        from varis.m1_ingestion.hgvs_parser import parse_hgvs
        record = create_variant_record("BRCA1", "not_valid_hgvs")
        record = parse_hgvs(record)
        assert record.residue_position is None
        assert "residue_position" in record.null_reasons

    def test_normalized_notation(self):
        """Should populate input_notation_normalized."""
        from varis.m1_ingestion.hgvs_parser import parse_hgvs
        record = create_variant_record("BRCA1", "p.R1699W")
        record = parse_hgvs(record)
        assert record.input_notation_normalized == "p.Arg1699Trp"


class TestUniProtClient:
    @pytest.mark.timeout(30)
    def test_fetch_brca1(self):
        """Fetch BRCA1 protein data from UniProt — real API call."""
        from varis.m1_ingestion.uniprot_client import fetch_uniprot
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record = fetch_uniprot(record)
        assert record.uniprot_id == "P38398"
        assert record.protein_sequence is not None
        assert len(record.protein_sequence) > 1000
        assert record.protein_length == len(record.protein_sequence)
        assert record.protein_name is not None

    @pytest.mark.timeout(30)
    def test_unknown_gene_returns_none(self):
        """Unknown gene should set fields to None with reason, not crash."""
        from varis.m1_ingestion.uniprot_client import fetch_uniprot
        record = create_variant_record("FAKEGENE999", "p.Ala1Val")
        record = fetch_uniprot(record)
        assert record.uniprot_id is None
        assert "uniprot_id" in record.null_reasons


class TestVariantNormalizer:
    def test_validates_correct_position(self):
        """When ref AA matches UniProt sequence, confidence should be 'high'."""
        from varis.m1_ingestion.variant_normalizer import normalize_variant
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record.residue_position = 1699
        record.ref_aa_single = "R"
        record.protein_sequence = "M" + "A" * 1697 + "R" + "A" * 164
        record = normalize_variant(record)
        assert record.coordinate_mapping_confidence == "high"
        assert record.uniprot_residue_position == 1699
        assert record.coordinate_mapping_method == "direct"

    def test_detects_mismatch(self):
        """When ref AA does NOT match, confidence should be 'failed'."""
        from varis.m1_ingestion.variant_normalizer import normalize_variant
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record.residue_position = 1699
        record.ref_aa_single = "R"
        record.protein_sequence = "M" + "A" * 1697 + "G" + "A" * 164
        record = normalize_variant(record)
        assert record.coordinate_mapping_confidence == "failed"
        assert record.normalization_warnings is not None
        assert len(record.normalization_warnings) > 0

    def test_no_sequence_upstream_failure(self):
        """If protein_sequence is None (UniProt failed), set upstream_dependency_failed."""
        from varis.m1_ingestion.variant_normalizer import normalize_variant
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record.residue_position = 1699
        record.ref_aa_single = "R"
        record.protein_sequence = None
        record = normalize_variant(record)
        assert record.coordinate_mapping_confidence is None
        assert "coordinate_mapping_confidence" in record.null_reasons

    def test_position_out_of_range(self):
        """Position beyond sequence length should set confidence to 'failed'."""
        from varis.m1_ingestion.variant_normalizer import normalize_variant
        record = create_variant_record("BRCA1", "p.Arg9999Trp")
        record.residue_position = 9999
        record.ref_aa_single = "R"
        record.protein_sequence = "MAAAA"
        record = normalize_variant(record)
        assert record.coordinate_mapping_confidence == "failed"


class TestClinVarClient:
    @pytest.mark.timeout(30)
    def test_fetch_known_variant(self):
        """Fetch BRCA1 p.Arg1699Trp from ClinVar — real API call."""
        from varis.m1_ingestion.clinvar_client import fetch_clinvar
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record = fetch_clinvar(record)
        assert record.clinvar_id is not None
        assert record.clinvar_classification is not None

    @pytest.mark.timeout(30)
    def test_extracts_genomic_coordinates(self):
        """ClinVar should provide genomic coordinates for gnomAD."""
        from varis.m1_ingestion.clinvar_client import fetch_clinvar
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record = fetch_clinvar(record)
        if record.clinvar_id is not None:
            assert record.reference_build is not None
            assert record.clinvar_chrom is not None
            assert record.clinvar_pos is not None
            assert record.coordinate_source == "clinvar"

    @pytest.mark.timeout(30)
    def test_unknown_variant_returns_none(self):
        """Unknown variant should set fields to None, not crash."""
        from varis.m1_ingestion.clinvar_client import fetch_clinvar
        record = create_variant_record("FAKEGENE999", "p.Ala1Val")
        record = fetch_clinvar(record)
        assert record.clinvar_id is None
        assert "clinvar_id" in record.null_reasons


class TestGnomADClient:
    @pytest.mark.timeout(30)
    def test_fetch_with_genomic_coords(self):
        """Fetch gnomAD frequency using ClinVar-provided genomic coordinates."""
        from varis.m1_ingestion.gnomad_client import fetch_gnomad
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record.reference_build = "GRCh38"
        record.clinvar_chrom = "17"
        record.clinvar_pos = 43057051
        record.clinvar_ref = "G"
        record.clinvar_alt = "A"
        record.coordinate_source = "clinvar"
        record = fetch_gnomad(record)
        if record.gnomad_frequency is not None:
            assert 0 <= record.gnomad_frequency <= 1

    @pytest.mark.timeout(10)
    def test_skips_without_genomic_coords(self):
        """Should skip and set reason when no genomic coordinates available."""
        from varis.m1_ingestion.gnomad_client import fetch_gnomad
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record = fetch_gnomad(record)
        assert record.gnomad_frequency is None
        assert record.null_reasons.get("gnomad_frequency") == "no_genomic_coordinates"

    @pytest.mark.timeout(10)
    def test_skips_wrong_build(self):
        """Should skip if reference build doesn't match gnomAD endpoint."""
        from varis.m1_ingestion.gnomad_client import fetch_gnomad
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record.reference_build = "GRCh37"
        record.clinvar_chrom = "17"
        record.clinvar_pos = 43057051
        record.clinvar_ref = "G"
        record.clinvar_alt = "A"
        record = fetch_gnomad(record)
        assert record.gnomad_frequency is None


class TestAlphaFoldClient:
    @pytest.mark.timeout(60)
    def test_download_brca1_structure(self):
        """Download BRCA1 AlphaFold structure — real API call."""
        from pathlib import Path
        from varis.m1_ingestion.alphafold_client import fetch_alphafold_structure
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record.uniprot_id = "P38398"
        record = fetch_alphafold_structure(record)
        assert record.structure_source == "alphafold"
        assert record.pdb_path is not None
        assert Path(record.pdb_path).exists()
        assert Path(record.pdb_path).stat().st_size > 0

    @pytest.mark.timeout(60)
    def test_caches_existing_file(self):
        """Second download should use cached file."""
        from varis.m1_ingestion.alphafold_client import fetch_alphafold_structure
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record.uniprot_id = "P38398"
        record = fetch_alphafold_structure(record)
        first_path = record.pdb_path

        record2 = create_variant_record("BRCA1", "p.Arg1699Trp")
        record2.uniprot_id = "P38398"
        record2 = fetch_alphafold_structure(record2)
        assert record2.pdb_path == first_path

    @pytest.mark.timeout(10)
    def test_no_uniprot_id(self):
        """Should skip when uniprot_id is not available."""
        from varis.m1_ingestion.alphafold_client import fetch_alphafold_structure
        record = create_variant_record("FAKEGENE", "p.Ala1Val")
        record = fetch_alphafold_structure(record)
        assert record.pdb_path is None
        assert "pdb_path" in record.null_reasons


class TestAlphaMissenseClient:
    @pytest.mark.timeout(30)
    def test_fetch_brca1_score(self):
        """Fetch AlphaMissense score for BRCA1 p.Arg1699Trp — real API call."""
        from varis.m1_ingestion.alphamissense_client import fetch_alphamissense
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record.uniprot_id = "P38398"
        record.residue_position = 1699
        record.ref_aa_single = "R"
        record.alt_aa_single = "W"
        record = fetch_alphamissense(record)
        if record.alphamissense_score is not None:
            assert 0.8 < record.alphamissense_score < 1.0
            assert record.alphamissense_class == "likely_pathogenic"

    @pytest.mark.timeout(10)
    def test_unknown_variant(self):
        """Unknown variant should not crash."""
        from varis.m1_ingestion.alphamissense_client import fetch_alphamissense
        record = create_variant_record("FAKEGENE", "p.Ala1Val")
        record.uniprot_id = "FAKE123"
        record.residue_position = 1
        record.ref_aa_single = "A"
        record.alt_aa_single = "V"
        record = fetch_alphamissense(record)
        assert record is not None


class TestM1Integration:
    @pytest.mark.timeout(120)
    def test_full_m1_brca1(self):
        """Run full M1 pipeline on BRCA1 p.Arg1699Trp with real API calls."""
        import json
        from varis.m1_ingestion import run

        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record = run(record)

        # HGVS parsing
        assert record.residue_position == 1699
        assert record.ref_amino_acid == "Arg"
        assert record.alt_amino_acid == "Trp"

        # UniProt
        assert record.uniprot_id == "P38398"
        assert record.protein_sequence is not None
        assert record.protein_length > 1000

        # Normalizer validation
        assert record.coordinate_mapping_confidence == "high"

        # ClinVar (should find this well-known variant)
        assert record.clinvar_id is not None

        # AlphaFold (should download structure)
        assert record.structure_source == "alphafold"
        assert record.pdb_path is not None

        # Module tracking
        assert "M1" in record.modules_completed

        # Should serialize to valid JSON
        json_str = record.to_json()
        data = json.loads(json_str)
        assert data["gene_symbol"] == "BRCA1"
        assert data["variant_id"] == "BRCA1_p.Arg1699Trp"
