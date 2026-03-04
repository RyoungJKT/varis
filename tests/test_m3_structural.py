"""Tests for M3: Structural Analysis — Feature extraction from 3D structure."""
import subprocess
import pytest
import shutil
from pathlib import Path
from unittest.mock import MagicMock, patch
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


class TestEvoEF2:
    """Tests for evoef2_wrapper.py — EvoEF2 ΔΔG prediction."""

    def test_evoef2_no_binary(self, m1_completed_record):
        """EvoEF2 not installed → ddg_evoef2=None, reason tool_missing."""
        from varis.m3_structural_analysis.evoef2_wrapper import run_evoef2
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        with patch("varis.m3_structural_analysis.evoef2_wrapper._find_binary", return_value=None):
            result = run_evoef2(m1_completed_record)
        assert result.ddg_evoef2 is None
        assert result.null_reasons.get("ddg_evoef2") == "tool_missing"

    def test_evoef2_no_structure(self, m1_completed_record):
        """No PDB → ddg_evoef2=None, reason upstream_dependency_failed."""
        from varis.m3_structural_analysis.evoef2_wrapper import run_evoef2
        m1_completed_record.pdb_path = None
        with patch("varis.m3_structural_analysis.evoef2_wrapper._find_binary", return_value="/usr/local/bin/EvoEF2"):
            result = run_evoef2(m1_completed_record)
        assert result.ddg_evoef2 is None
        assert result.null_reasons.get("ddg_evoef2") == "upstream_dependency_failed"

    def test_evoef2_mutation_string_format(self):
        """Verify mutation string format: RA1699W;"""
        from varis.m3_structural_analysis.evoef2_wrapper import _build_mutation_string
        result = _build_mutation_string("R", "A", 1699, "W")
        assert result == "RA1699W;"

    def test_evoef2_parse_total_energy(self):
        """Parse 'Total = -456.78' from ComputeStability output."""
        from varis.m3_structural_analysis.evoef2_wrapper import _parse_total_energy
        stdout = "Reference ALA:   12.34\nTotal           =     -456.78\nDone."
        result = _parse_total_energy(stdout)
        assert result == pytest.approx(-456.78)

    def test_evoef2_parse_total_energy_not_found(self):
        """Returns None when Total line is missing."""
        from varis.m3_structural_analysis.evoef2_wrapper import _parse_total_energy
        assert _parse_total_energy("no energy line here") is None

    def test_evoef2_timeout(self, m1_completed_record):
        """Subprocess timeout → ddg_evoef2=None, reason timed_out."""
        from varis.m3_structural_analysis.evoef2_wrapper import run_evoef2
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        with patch("varis.m3_structural_analysis.evoef2_wrapper._find_binary", return_value="/usr/local/bin/EvoEF2"):
            with patch("varis.m3_structural_analysis.evoef2_wrapper.subprocess.run",
                        side_effect=subprocess.TimeoutExpired("EvoEF2", 120)):
                result = run_evoef2(m1_completed_record)
        assert result.ddg_evoef2 is None
        assert result.null_reasons.get("ddg_evoef2") == "timed_out"

    def test_evoef2_successful_run(self, m1_completed_record):
        """Mock subprocess calls → verify ddg_evoef2 is computed."""
        from varis.m3_structural_analysis.evoef2_wrapper import run_evoef2

        m1_completed_record.pdb_path = str(BRCA1_PDB)

        # Create a real PDB file if it doesn't exist (for copy step)
        pdb_dir = Path(m1_completed_record.pdb_path).parent
        pdb_dir.mkdir(parents=True, exist_ok=True)
        if not Path(m1_completed_record.pdb_path).exists():
            Path(m1_completed_record.pdb_path).write_text("ATOM mock PDB\nEND\n")

        def mock_subprocess_run(cmd, **kwargs):
            """Mock EvoEF2 subprocess calls."""
            cwd = Path(kwargs.get("cwd", "."))
            cmd_str = " ".join(cmd)
            mock_result = MagicMock()
            mock_result.returncode = 0
            mock_result.stderr = ""

            if "RepairStructure" in cmd_str:
                # Create the repaired PDB file
                for f in cwd.glob("*.pdb"):
                    repaired = cwd / f"{f.stem}_Repair.pdb"
                    repaired.write_text("ATOM repaired\nEND\n")
                    break
                mock_result.stdout = "Repair done"
            elif "BuildMutant" in cmd_str:
                # Create the mutant PDB file
                for f in cwd.glob("*_Repair.pdb"):
                    mutant = cwd / f"{f.stem}_Model_0001.pdb"
                    mutant.write_text("ATOM mutant\nEND\n")
                    break
                mock_result.stdout = "Mutant built"
            elif "ComputeStability" in cmd_str:
                if "Model_0001" in cmd_str:
                    mock_result.stdout = "Total           =     -450.00\n"
                else:
                    mock_result.stdout = "Total           =     -453.50\n"
            return mock_result

        with patch("varis.m3_structural_analysis.evoef2_wrapper._find_binary", return_value="/usr/local/bin/EvoEF2"):
            with patch("varis.m3_structural_analysis.evoef2_wrapper.subprocess.run", side_effect=mock_subprocess_run):
                result = run_evoef2(m1_completed_record)

        assert result.ddg_evoef2 is not None
        # DDG = mutant (-450.00) - wt (-453.50) = 3.50
        assert result.ddg_evoef2 == pytest.approx(3.5, abs=0.01)


class TestFoldXStub:
    """Tests for foldx_wrapper.py — ΔΔG stub."""

    def test_foldx_no_binary(self, m1_completed_record):
        """FoldX not installed → ddg_foldx=None, reason tool_missing."""
        from varis.m3_structural_analysis.foldx_wrapper import run_foldx
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        result = run_foldx(m1_completed_record)
        assert result.ddg_foldx is None
        assert result.null_reasons.get("ddg_foldx") == "tool_missing"


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


class TestM3Orchestrator:
    """Tests for M3 orchestration — site-gating and full pipeline."""

    def test_m3_skips_when_site_absent(self, m1_completed_record):
        """All site-dependent features skipped when mutation_site_present=False."""
        from varis.m3_structural_analysis import run
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.mutation_site_present = False
        result = run(m1_completed_record)
        # Site-dependent features should be explicitly skipped
        assert result.sasa_available is False
        assert result.dssp_available is False
        assert result.contacts_available is False
        # ddg_available is set by orchestrator based on ddg_mean
        assert result.ddg_available is not None
        # M3 should still complete
        assert "M3" in result.modules_completed

    def test_m3_no_structure(self, m1_completed_record):
        """M3 with no structure marks failed."""
        from varis.m3_structural_analysis import run
        m1_completed_record.pdb_path = None
        result = run(m1_completed_record)
        assert "M3" in result.modules_failed

    def test_m3_integration_golden_record(self, m1_completed_record):
        """Full M3 pipeline — verify all feature availability flags are set."""
        from varis.m3_structural_analysis import run
        m1_completed_record.pdb_path = str(BRCA1_PDB)
        m1_completed_record.mutation_site_present = True
        result = run(m1_completed_record)
        assert "M3" in result.modules_completed
        for flag in ["sasa_available", "dssp_available", "contacts_available",
                      "ddg_available", "domain_available"]:
            val = getattr(result, flag)
            assert val is not None, f"{flag} should be True or False, not None"
