"""Phase 2 Integration Test — Full M2 -> M3 pipeline on BRCA1 p.Arg1699Trp.

Uses the real AlphaFold PDB file already downloaded by M1.
Mocks external API calls (InterPro) for offline testing.
"""
import pytest
from unittest.mock import MagicMock, patch

from varis.config import STRUCTURES_DIR
from varis.models.variant_record import create_variant_record

BRCA1_PDB = STRUCTURES_DIR / "AF-P38398-F1-model_v6.pdb"


def _make_interpro_mock_client(with_domain: bool = True) -> MagicMock:
    """Create a mock httpx.Client that returns BRCT domain data for InterPro.

    Args:
        with_domain: If True, returns a BRCT domain hit at 1646-1736.
            If False, returns empty results (no domain found).

    Returns:
        A MagicMock configured as an httpx.Client context manager.
    """
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.raise_for_status = MagicMock()

    if with_domain:
        mock_response.json.return_value = {
            "results": [{
                "metadata": {
                    "accession": "PF00533",
                    "name": "BRCT",
                    "type": "domain",
                },
                "proteins": [{
                    "entry_protein_locations": [{
                        "fragments": [{"start": 1646, "end": 1736}]
                    }]
                }]
            }]
        }
    else:
        mock_response.json.return_value = {"results": []}

    mock_client = MagicMock()
    mock_client.get.return_value = mock_response
    return mock_client


class TestPhase2Integration:
    """End-to-end M2 -> M3 pipeline test."""

    @pytest.fixture
    def brca1_m1_record(self):
        """Create a record that has completed M1 with real AlphaFold PDB."""
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        record.residue_position = 1699
        record.ref_amino_acid = "Arg"
        record.alt_amino_acid = "Trp"
        record.ref_aa_single = "R"
        record.alt_aa_single = "W"
        record.charge_change = "+ve \u2192 neutral"
        record.uniprot_id = "P38398"
        record.protein_name = "Breast cancer type 1 susceptibility protein"
        record.protein_length = 1863
        record.protein_sequence = "M" * 1863  # placeholder sequence
        record.structure_source = "alphafold"
        record.pdb_path = str(BRCA1_PDB)
        record.mark_module_completed("M1")
        return record

    def test_m2_then_m3_pipeline(self, brca1_m1_record):
        """Full M2 -> M3 pipeline produces expected outputs."""
        from varis.m2_structure import run as run_m2
        from varis.m3_structural_analysis import run as run_m3

        # ---- M2: Structure validation and quality extraction ----
        record = run_m2(brca1_m1_record)
        assert "M2" in record.modules_completed
        assert record.mutation_site_present is True
        assert record.pdb_hash is not None
        assert len(record.pdb_hash) == 64  # SHA-256 hex digest
        assert record.numbering_scheme == "uniprot_canonical"
        assert record.plddt_available is True
        assert record.mutation_site_plddt is not None
        assert 0.0 <= record.mutation_site_plddt <= 100.0

        # ---- M3: Structural feature extraction (InterPro mocked) ----
        mock_client = _make_interpro_mock_client(with_domain=True)

        with patch(
            "varis.m3_structural_analysis.interpro_client.httpx.Client",
            return_value=mock_client,
        ):
            record = run_m3(record)

        assert "M3" in record.modules_completed

        # Verify all feature availability flags are set (not None)
        assert record.sasa_available is not None
        assert record.dssp_available is not None
        assert record.contacts_available is not None
        assert record.ddg_available is not None
        assert record.domain_available is not None

        # If FreeSASA is installed, check results make physical sense
        if record.sasa_available:
            assert 0.0 <= record.solvent_accessibility_relative <= 1.0
            assert record.burial_category in ("core", "surface")

        # If DSSP ran, check results
        if record.dssp_available:
            assert record.secondary_structure is not None
            assert record.secondary_structure_name in ("helix", "sheet", "coil")

        # Contacts should work (pure BioPython, always available)
        assert record.contacts_available is True
        assert record.contacts_wt > 0
        assert record.hbonds_wt >= 0
        assert record.packing_density > 0.0

        # DDG should be unavailable (no FoldX/PyRosetta in test environment)
        assert record.ddg_available is False

        # Domain should be available (mocked InterPro returned BRCT)
        assert record.domain_available is True
        assert record.domain_name == "BRCT"
        assert record.domain_id == "PF00533"
        assert record.domain_criticality == "critical"
        assert record.domain_start == 1646
        assert record.domain_end == 1736

    def test_golden_record_schema(self, brca1_m1_record):
        """After M2+M3, record serializes to JSON with all expected keys."""
        from varis.m2_structure import run as run_m2
        from varis.m3_structural_analysis import run as run_m3

        record = run_m2(brca1_m1_record)

        # Mock InterPro with empty results (no domain)
        mock_client = _make_interpro_mock_client(with_domain=False)

        with patch(
            "varis.m3_structural_analysis.interpro_client.httpx.Client",
            return_value=mock_client,
        ):
            record = run_m3(record)

        data = record.to_dict()

        # M2 keys must be present
        m2_keys = [
            "mutation_site_present", "pdb_hash", "numbering_scheme",
            "preparation_steps", "plddt_available",
        ]
        for key in m2_keys:
            assert key in data, f"Missing M2 key: {key}"
            assert data[key] is not None, f"M2 key {key} should not be None"

        # Feature availability keys must be present and non-None
        avail_keys = [
            "sasa_available", "dssp_available", "contacts_available",
            "ddg_available", "domain_available",
        ]
        for key in avail_keys:
            assert key in data, f"Missing availability key: {key}"
            assert data[key] is not None, f"{key} should not be None after M3"

        # Schema version and modules must be present
        assert data["record_schema_version"] is not None
        assert "M1" in data["modules_completed"]
        assert "M2" in data["modules_completed"]
        assert "M3" in data["modules_completed"]

    def test_golden_record_json_roundtrip(self, brca1_m1_record):
        """Record survives JSON serialization and deserialization."""
        from varis.m2_structure import run as run_m2
        from varis.m3_structural_analysis import run as run_m3
        from varis.models.variant_record import VariantRecord

        record = run_m2(brca1_m1_record)

        mock_client = _make_interpro_mock_client(with_domain=True)

        with patch(
            "varis.m3_structural_analysis.interpro_client.httpx.Client",
            return_value=mock_client,
        ):
            record = run_m3(record)

        # Round-trip through JSON
        json_str = record.to_json()
        restored = VariantRecord.from_json(json_str)

        # Core identity preserved
        assert restored.gene_symbol == "BRCA1"
        assert restored.hgvs_protein == "p.Arg1699Trp"
        assert restored.residue_position == 1699

        # Module completion preserved
        assert "M1" in restored.modules_completed
        assert "M2" in restored.modules_completed
        assert "M3" in restored.modules_completed

        # Feature availability preserved
        assert restored.contacts_available is True
        assert restored.ddg_available is False
        assert restored.domain_available is True

        # Structural features preserved (contacts always work)
        assert restored.contacts_wt == record.contacts_wt
        assert restored.hbonds_wt == record.hbonds_wt

    def test_m3_skips_site_tools_when_site_absent(self, brca1_m1_record):
        """When M2 reports mutation_site_present=False, site tools are skipped."""
        from varis.m3_structural_analysis import run as run_m3

        # Simulate M2 finding site is absent
        brca1_m1_record.mutation_site_present = False
        brca1_m1_record.mark_module_completed("M2")

        mock_client = _make_interpro_mock_client(with_domain=False)

        with patch(
            "varis.m3_structural_analysis.interpro_client.httpx.Client",
            return_value=mock_client,
        ):
            record = run_m3(brca1_m1_record)

        # M3 still completes
        assert "M3" in record.modules_completed

        # Site-dependent tools are explicitly skipped
        assert record.sasa_available is False
        assert record.dssp_available is False
        assert record.contacts_available is False

        # DDG stubs still run (site-independent in orchestrator)
        assert record.ddg_available is False

    def test_ml_features_populated(self, brca1_m1_record):
        """get_ml_features() returns non-empty dict after M2+M3."""
        from varis.m2_structure import run as run_m2
        from varis.m3_structural_analysis import run as run_m3

        record = run_m2(brca1_m1_record)

        mock_client = _make_interpro_mock_client(with_domain=True)

        with patch(
            "varis.m3_structural_analysis.interpro_client.httpx.Client",
            return_value=mock_client,
        ):
            record = run_m3(record)

        ml_features = record.get_ml_features()

        # Should contain structural features + availability flags
        assert "contacts_wt" in ml_features
        assert "ddg_available" in ml_features
        assert "sasa_available" in ml_features
        assert "domain_available" in ml_features
        assert "contacts_available" in ml_features
        assert "mutation_site_plddt" in ml_features
        assert "charge_change" in ml_features

        # Contacts should have a real value
        assert ml_features["contacts_wt"] is not None
        assert ml_features["contacts_wt"] > 0

        # Availability flags should be booleans
        assert isinstance(ml_features["ddg_available"], bool)
        assert isinstance(ml_features["contacts_available"], bool)

        # Count available features -- at least contacts + domain + pLDDT + charge
        available_count = record.count_available_features()
        assert available_count >= 4, (
            f"Expected at least 4 available features, got {available_count}"
        )
