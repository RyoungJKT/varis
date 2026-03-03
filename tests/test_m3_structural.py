"""Tests for M3: Structural Analysis — Feature extraction from 3D structure."""
import pytest
import shutil
from pathlib import Path
from unittest.mock import MagicMock
from varis.config import STRUCTURES_DIR

BRCA1_PDB = STRUCTURES_DIR / "AF-P38398-F1-model_v6.pdb"


class TestFreeSASA:
    """Tests for freesasa_wrapper.py — solvent accessibility."""

    def test_freesasa_computes_sasa(self, m1_completed_record):
        """Returns float, sasa_available=True."""
        from varis.m3_structural_analysis.freesasa_wrapper import run_freesasa
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        result = run_freesasa(m1_completed_record)
        assert result.sasa_available is True
        assert isinstance(result.solvent_accessibility_relative, float)

    def test_freesasa_relative_bounds(self, m1_completed_record):
        """Relative SASA in [0.0, 1.0]."""
        from varis.m3_structural_analysis.freesasa_wrapper import run_freesasa
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        result = run_freesasa(m1_completed_record)
        assert 0.0 <= result.solvent_accessibility_relative <= 1.0

    def test_freesasa_burial_category(self, m1_completed_record):
        """Returns 'core' or 'surface' only (no 'interface')."""
        from varis.m3_structural_analysis.freesasa_wrapper import run_freesasa
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        result = run_freesasa(m1_completed_record)
        assert result.burial_category in ("core", "surface")

    def test_freesasa_no_structure(self, m1_completed_record):
        """No PDB -> sasa_available=False with reason."""
        from varis.m3_structural_analysis.freesasa_wrapper import run_freesasa
        m1_completed_record.pdb_path = None
        result = run_freesasa(m1_completed_record)
        assert result.sasa_available is False
        assert result.sasa_missing_reason is not None


class TestDSSP:
    """Tests for dssp_wrapper.py — secondary structure assignment."""

    def test_dssp_returns_valid_code(self, m1_completed_record):
        """Result in allowed DSSP code set."""
        if not shutil.which("mkdssp"):
            pytest.skip("mkdssp binary not installed")
        from varis.m3_structural_analysis.dssp_wrapper import run_dssp
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        result = run_dssp(m1_completed_record)
        if result.dssp_available:
            valid_codes = {"H", "B", "E", "G", "I", "T", "S", "-", "C"}
            assert result.secondary_structure in valid_codes
            assert result.secondary_structure_name in ("helix", "sheet", "coil")

    def test_dssp_mkdssp_missing(self, m1_completed_record, monkeypatch):
        """dssp_available=False when mkdssp not found."""
        from varis.m3_structural_analysis.dssp_wrapper import run_dssp
        monkeypatch.setattr(shutil, "which", lambda x: None)
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        result = run_dssp(m1_completed_record)
        assert result.dssp_available is False
        assert result.dssp_missing_reason is not None

    def test_dssp_no_structure(self, m1_completed_record):
        """No PDB → dssp_available=False."""
        from varis.m3_structural_analysis.dssp_wrapper import run_dssp
        m1_completed_record.pdb_path = None
        result = run_dssp(m1_completed_record)
        assert result.dssp_available is False


class TestContacts:
    """Tests for biopython_contacts.py — WT local environment."""

    def test_contacts_wt_valid(self, m1_completed_record):
        """Integers >= 0, nonzero for folded region."""
        from varis.m3_structural_analysis.biopython_contacts import run_contacts
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        result = run_contacts(m1_completed_record)
        if result.contacts_available:
            assert isinstance(result.contacts_wt, int)
            assert result.contacts_wt >= 0
            assert isinstance(result.hbonds_wt, int)
            assert result.hbonds_wt >= 0
            assert result.contacts_wt > 0  # Folded region should have contacts

    def test_contacts_heavy_atoms_only(self, m1_completed_record):
        """Contacts count uses heavy atoms only (no hydrogens)."""
        from varis.m3_structural_analysis.biopython_contacts import run_contacts
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        result = run_contacts(m1_completed_record)
        if result.contacts_available:
            assert result.contacts_wt < 200  # Sanity upper bound for heavy atoms

    def test_contacts_packing_density(self, m1_completed_record):
        """Packing density is a float >= 0."""
        from varis.m3_structural_analysis.biopython_contacts import run_contacts
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        result = run_contacts(m1_completed_record)
        if result.contacts_available:
            assert isinstance(result.packing_density, float)
            assert result.packing_density >= 0.0

    def test_contacts_no_structure(self, m1_completed_record):
        """No PDB → contacts_available=False."""
        from varis.m3_structural_analysis.biopython_contacts import run_contacts
        m1_completed_record.pdb_path = None
        result = run_contacts(m1_completed_record)
        assert result.contacts_available is False


class TestFoldXStub:
    """Tests for foldx_wrapper.py — ΔΔG stub."""

    def test_foldx_no_binary(self, m1_completed_record):
        """FoldX not installed → ddg_available=False, reason tool_missing."""
        from varis.m3_structural_analysis.foldx_wrapper import run_foldx
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        result = run_foldx(m1_completed_record)
        assert result.ddg_available is False
        assert result.ddg_missing_reason == "tool_missing"


class TestPyRosettaStub:
    """Tests for pyrosetta_wrapper.py — ΔΔG stub."""

    def test_pyrosetta_not_installed(self, m1_completed_record):
        """PyRosetta not installed → ddg_pyrosetta is None."""
        from varis.m3_structural_analysis.pyrosetta_wrapper import run_pyrosetta
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        result = run_pyrosetta(m1_completed_record)
        assert result.ddg_pyrosetta is None


class TestInterPro:
    """Tests for interpro_client.py — domain identification via InterPro API."""

    def test_interpro_brca1_domain(self, m1_completed_record):
        """Returns BRCT domain with boundaries for position 1699 (mocked)."""
        from varis.m3_structural_analysis.interpro_client import run_interpro
        mock_response = MagicMock()
        mock_response.status_code = 200
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
        mock_response.raise_for_status = MagicMock()
        mock_client = MagicMock()
        mock_client.get.return_value = mock_response
        result = run_interpro(m1_completed_record, client=mock_client)
        assert result.domain_available is True
        assert result.domain_name == "BRCT"
        assert result.domain_id == "PF00533"
        assert result.domain_start == 1646
        assert result.domain_end == 1736

    def test_interpro_position_outside_domain(self, m1_completed_record):
        """Position outside any domain -> domain_available=False (mocked)."""
        from varis.m3_structural_analysis.interpro_client import run_interpro
        m1_completed_record.residue_position = 50
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "results": [{
                "metadata": {"accession": "PF00533", "name": "BRCT", "type": "domain"},
                "proteins": [{"entry_protein_locations": [{"fragments": [{"start": 1646, "end": 1736}]}]}]
            }]
        }
        mock_response.raise_for_status = MagicMock()
        mock_client = MagicMock()
        mock_client.get.return_value = mock_response
        result = run_interpro(m1_completed_record, client=mock_client)
        assert result.domain_available is False

    def test_interpro_no_uniprot_id(self, m1_completed_record):
        """No UniProt ID -> domain_available=False."""
        from varis.m3_structural_analysis.interpro_client import run_interpro
        m1_completed_record.uniprot_id = None
        result = run_interpro(m1_completed_record)
        assert result.domain_available is False

    def test_interpro_stores_boundaries(self, m1_completed_record):
        """Domain start/end boundaries are stored."""
        from varis.m3_structural_analysis.interpro_client import run_interpro
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "results": [{
                "metadata": {"accession": "PF00533", "name": "BRCT", "type": "domain"},
                "proteins": [{"entry_protein_locations": [{"fragments": [{"start": 1646, "end": 1736}]}]}]
            }]
        }
        mock_response.raise_for_status = MagicMock()
        mock_client = MagicMock()
        mock_client.get.return_value = mock_response
        result = run_interpro(m1_completed_record, client=mock_client)
        assert result.domain_start is not None
        assert result.domain_end is not None
        assert result.domain_start < result.domain_end
