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
