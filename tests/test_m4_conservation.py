"""Tests for M4: Conservation Analysis."""
import json
import pytest
import math
from unittest.mock import MagicMock, patch
from varis.models.variant_record import (
    create_variant_record, NullReason, RECORD_SCHEMA_VERSION,
)


class TestSchemaV130:
    """Verify schema v1.4.0 fields exist on VariantRecord."""

    def test_schema_version_is_1_4_0(self):
        assert RECORD_SCHEMA_VERSION == "1.5.0"

    def test_new_msa_fields_exist(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert hasattr(record, "msa_num_sequences")
        assert hasattr(record, "msa_gap_fraction_at_site")
        assert hasattr(record, "msa_column_index")

    def test_insufficient_data_reason_exists(self):
        assert hasattr(NullReason, "INSUFFICIENT_DATA")
        assert NullReason.INSUFFICIENT_DATA == "insufficient_data"


class TestConservationScorer:
    """Tests for conservation_scorer.py — entropy and position mapping."""

    def test_entropy_fully_conserved(self):
        """All same AA -> entropy=0, score=1.0."""
        from varis.m4_conservation.conservation_scorer import _shannon_entropy
        column = ["R"] * 20
        assert _shannon_entropy(column) == pytest.approx(0.0, abs=0.001)

    def test_entropy_maximally_variable(self):
        """All different AAs -> entropy ~ log2(20), score ~ 0.0."""
        from varis.m4_conservation.conservation_scorer import _shannon_entropy
        aas = list("ACDEFGHIKLMNPQRSTVWY")
        assert _shannon_entropy(aas) == pytest.approx(math.log2(20), abs=0.01)

    def test_entropy_ignores_gaps(self):
        """Gaps excluded from frequency, but present in input."""
        from varis.m4_conservation.conservation_scorer import _shannon_entropy
        column = ["R"] * 18 + ["-", "-"]
        assert _shannon_entropy(column) == pytest.approx(0.0, abs=0.001)

    def test_position_mapping_no_gaps(self):
        """Direct mapping when query has no gaps."""
        from varis.m4_conservation.conservation_scorer import _map_position_to_column
        query_row = "MKRST"
        # Position is 1-indexed: M=pos1(col0), K=pos2(col1), R=pos3(col2)
        col = _map_position_to_column(query_row, position=3)
        assert col == 2

    def test_position_mapping_with_gaps(self):
        """Gaps in query shift the column index."""
        from varis.m4_conservation.conservation_scorer import _map_position_to_column
        query_row = "M-KR-ST"
        # M=pos1(col0), K=pos2(col2), R=pos3(col3)
        col = _map_position_to_column(query_row, position=3)
        assert col == 3

    def test_position_mapping_validates_ref_aa(self):
        """Returns None if ref AA doesn't match."""
        from varis.m4_conservation.conservation_scorer import _map_position_to_column
        query_row = "MKRST"
        col = _map_position_to_column(query_row, position=3, expected_aa="W")
        assert col is None  # R != W

    def test_score_conservation_full(self, m1_completed_record):
        """Full scorer with a small test alignment — all conserved."""
        from varis.m4_conservation.conservation_scorer import score_conservation
        alignment = {
            "sequences": {
                "query": "MKRST",
                **{f"orth{i}": "MKRST" for i in range(1, 11)},
            },
            "query_id": "query",
            "taxonomy": {
                "orth1": 9606, "orth2": 9615, "orth3": 10090,
                "orth4": 9913, "orth5": 9823,
                "orth6": 7955, "orth7": 8364, "orth8": 9031,
                "orth9": 28377, "orth10": 13616,
            },
        }
        m1_completed_record.residue_position = 3
        m1_completed_record.ref_aa_single = "R"
        result = score_conservation(m1_completed_record, alignment)
        assert result.conservation_score == pytest.approx(1.0, abs=0.01)
        assert result.position_entropy == pytest.approx(0.0, abs=0.01)
        assert result.msa_column_index == 2
        assert result.conservation_available is True

    def test_mammal_conservation_threshold(self, m1_completed_record):
        """>=90% mammals have ref AA -> True."""
        from varis.m4_conservation.conservation_scorer import score_conservation
        alignment = {
            "sequences": {
                "query": "MKRST",
                "m1": "MKRST", "m2": "MKRST", "m3": "MKRST",
                "m4": "MKRST", "m5": "MKRST",
                "o1": "MKWST",
            },
            "query_id": "query",
            "taxonomy": {
                "m1": 9606, "m2": 9615, "m3": 10090, "m4": 9913, "m5": 9823,
                "o1": 7955,
            },
        }
        m1_completed_record.residue_position = 3
        m1_completed_record.ref_aa_single = "R"
        result = score_conservation(m1_completed_record, alignment)
        assert result.conserved_across_mammals is True

    def test_mammal_conservation_too_few(self, m1_completed_record):
        """<5 mammals -> conserved_across_mammals=None."""
        from varis.m4_conservation.conservation_scorer import score_conservation
        alignment = {
            "sequences": {
                "query": "MKRST",
                "m1": "MKRST", "m2": "MKRST",
                **{f"o{i}": "MKRST" for i in range(1, 9)},
            },
            "query_id": "query",
            "taxonomy": {"m1": 9606, "m2": 9615,
                          **{f"o{i}": 7955 + i for i in range(1, 9)}},
        }
        m1_completed_record.residue_position = 3
        m1_completed_record.ref_aa_single = "R"
        result = score_conservation(m1_completed_record, alignment)
        assert result.conserved_across_mammals is None


class TestUniProtOrthologs:
    """Tests for uniprot_orthologs.py — fetch ortholog sequences."""

    def test_fetch_orthologs_success(self, m1_completed_record):
        """Mocked: returns sequences with taxonomy."""
        from varis.m4_conservation.uniprot_orthologs import fetch_orthologs
        # Build mock FASTA with 15 sequences
        fasta_entries = []
        taxon_ids = [9606, 9598, 9544, 9615, 9913, 9823, 10090, 10116,
                     9986, 9685, 7955, 8364, 9031, 28377, 13616]
        for i, taxon in enumerate(taxon_ids):
            fasta_entries.append(f">sp|P{i:05d}|ORTH{i} OS=Species OX={taxon}\nMKRST\n")
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = "".join(fasta_entries)
        mock_client = MagicMock()
        mock_client.get.return_value = mock_response
        record, orthologs = fetch_orthologs(m1_completed_record, client=mock_client)
        assert orthologs is not None
        assert len(orthologs["sequences"]) >= 10
        assert "query_id" in orthologs
        assert "taxonomy" in orthologs

    def test_fetch_orthologs_too_few(self, m1_completed_record):
        """<10 sequences -> returns None."""
        from varis.m4_conservation.uniprot_orthologs import fetch_orthologs
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = ">sp|P00001|ORTH1 OX=9606\nMKRST\n"
        mock_client = MagicMock()
        mock_client.get.return_value = mock_response
        record, orthologs = fetch_orthologs(m1_completed_record, client=mock_client)
        assert orthologs is None

    def test_fetch_orthologs_no_uniprot(self, m1_completed_record):
        """No uniprot_id -> skip."""
        from varis.m4_conservation.uniprot_orthologs import fetch_orthologs
        m1_completed_record.uniprot_id = None
        record, orthologs = fetch_orthologs(m1_completed_record)
        assert orthologs is None

    def test_fetch_orthologs_caps_at_100(self, m1_completed_record):
        """Max 100 orthologs returned (plus query)."""
        from varis.m4_conservation.uniprot_orthologs import fetch_orthologs
        fasta = "".join(f">sp|P{i:05d}|ORTH{i} OX={9606+i}\nMKRST\n" for i in range(150))
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.text = fasta
        mock_client = MagicMock()
        mock_client.get.return_value = mock_response
        record, orthologs = fetch_orthologs(m1_completed_record, client=mock_client)
        assert orthologs is not None
        # 100 orthologs + 1 query
        assert len(orthologs["sequences"]) <= 101


class TestClustalClient:
    """Tests for clustal_client.py — EBI Clustal Omega API."""

    def test_clustal_alignment(self, m1_completed_record):
        """Mocked: submit + poll returns valid alignment."""
        from varis.m4_conservation.clustal_client import run_alignment
        orthologs = {
            "sequences": {
                "query": "MKRST",
                "orth1": "MKRST",
                "orth2": "MKRAT",
            },
            "query_id": "query",
            "taxonomy": {"orth1": 9606, "orth2": 10090},
        }
        # Mock submit (POST -> job ID)
        mock_submit = MagicMock()
        mock_submit.status_code = 200
        mock_submit.text = "clustalo-R20260303-123456"
        # Mock status (GET -> FINISHED)
        mock_status = MagicMock()
        mock_status.status_code = 200
        mock_status.text = "FINISHED"
        # Mock result (GET -> aligned FASTA)
        mock_result = MagicMock()
        mock_result.status_code = 200
        mock_result.text = ">query\nMKRST\n>orth1\nMKRST\n>orth2\nMKRAT\n"

        mock_client = MagicMock()
        mock_client.post.return_value = mock_submit
        mock_client.get.side_effect = [mock_status, mock_result]

        record, alignment = run_alignment(m1_completed_record, orthologs, client=mock_client)
        assert alignment is not None
        assert "sequences" in alignment
        assert len(alignment["sequences"]) == 3
        assert alignment["taxonomy"] == {"orth1": 9606, "orth2": 10090}

    def test_clustal_timeout(self, m1_completed_record):
        """Poll exceeds max retries -> returns None."""
        from varis.m4_conservation.clustal_client import run_alignment
        orthologs = {
            "sequences": {"query": "MKRST", "orth1": "MKRST"},
            "query_id": "query",
            "taxonomy": {},
        }
        mock_submit = MagicMock()
        mock_submit.status_code = 200
        mock_submit.text = "clustalo-R20260303-123456"
        mock_status = MagicMock()
        mock_status.status_code = 200
        mock_status.text = "RUNNING"  # Never finishes

        mock_client = MagicMock()
        mock_client.post.return_value = mock_submit
        mock_client.get.return_value = mock_status

        record, alignment = run_alignment(
            m1_completed_record, orthologs, client=mock_client,
            max_polls=2, poll_interval=0,
        )
        assert alignment is None

    def test_clustal_no_sequences(self, m1_completed_record):
        """No orthologs -> skip."""
        from varis.m4_conservation.clustal_client import run_alignment
        record, alignment = run_alignment(m1_completed_record, None)
        assert alignment is None


class TestConSurfFallback:
    """Tests for consurf_fallback.py — pre-computed conservation."""

    def test_consurf_known_protein(self, m1_completed_record):
        """Mocked: returns grade, mapped to score."""
        from varis.m4_conservation.consurf_fallback import fetch_consurf
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "grades": {str(i): {"grade": 5} for i in range(1, 1864)},
        }
        mock_response.json.return_value["grades"]["1699"] = {"grade": 9}
        mock_client = MagicMock()
        mock_client.get.return_value = mock_response
        result = fetch_consurf(m1_completed_record, client=mock_client)
        assert result.conservation_score is not None
        assert result.conservation_score == pytest.approx(1.0, abs=0.01)  # Grade 9 -> (9-1)/8 = 1.0
        assert result.conservation_method == "consurf"
        assert result.conservation_available is True

    def test_consurf_unknown_protein(self, m1_completed_record):
        """Mocked: 404 -> conservation_available=False."""
        from varis.m4_conservation.consurf_fallback import fetch_consurf
        mock_response = MagicMock()
        mock_response.status_code = 404
        mock_client = MagicMock()
        mock_client.get.return_value = mock_response
        result = fetch_consurf(m1_completed_record, client=mock_client)
        assert result.conservation_available is False

    def test_consurf_no_uniprot_id(self, m1_completed_record):
        """No uniprot_id -> skip."""
        from varis.m4_conservation.consurf_fallback import fetch_consurf
        m1_completed_record.uniprot_id = None
        result = fetch_consurf(m1_completed_record)
        assert result.conservation_available is False


class TestM4Orchestrator:
    """Tests for M4 orchestration — caching, fallback, integration."""

    def test_m4_no_sequence(self, m1_completed_record, tmp_path):
        """No protein_sequence and no uniprot_id -> M4 fails gracefully."""
        from varis.m4_conservation import run
        m1_completed_record.protein_sequence = None
        m1_completed_record.uniprot_id = None
        with patch("varis.m4_conservation._CACHE_DIR", tmp_path / "conservation"):
            result = run(m1_completed_record)
        assert "M4" in result.modules_failed

    def test_m4_fallback_to_consurf(self, m1_completed_record, tmp_path):
        """Primary fails -> ConSurf runs."""
        from varis.m4_conservation import run
        cache_dir = tmp_path / "conservation"
        with patch("varis.m4_conservation._CACHE_DIR", cache_dir):
            with patch("varis.m4_conservation.uniprot_orthologs.fetch_orthologs",
                        return_value=(m1_completed_record, None)):
                mock_response = MagicMock()
                mock_response.status_code = 200
                mock_response.json.return_value = {
                    "grades": {"1699": {"grade": 8}},
                }
                mock_client = MagicMock()
                mock_client.get.return_value = mock_response
                with patch("varis.m4_conservation.consurf_fallback.httpx.Client",
                            return_value=mock_client):
                    result = run(m1_completed_record)
        assert result.conservation_score is not None
        assert result.conservation_method == "consurf"
        assert "M4" in result.modules_completed

    def test_m4_integration(self, m1_completed_record, tmp_path):
        """Full pipeline with mocks — golden record check."""
        from varis.m4_conservation import run
        # Build sequence with R at position 1699 (0-indexed: 1698)
        seq = "M" * 1698 + "R" + "M" * (1863 - 1699)
        orthologs = {
            "sequences": {
                "query": seq,
                **{f"orth{i}": seq for i in range(15)},
            },
            "query_id": "query",
            "taxonomy": {f"orth{i}": 9606 + i for i in range(15)},
        }
        cache_dir = tmp_path / "conservation"
        with patch("varis.m4_conservation._CACHE_DIR", cache_dir):
            with patch("varis.m4_conservation.uniprot_orthologs.fetch_orthologs",
                        return_value=(m1_completed_record, orthologs)):
                with patch("varis.m4_conservation.clustal_client.run_alignment",
                            return_value=(m1_completed_record, orthologs)):
                    result = run(m1_completed_record)
        assert "M4" in result.modules_completed
        assert result.conservation_available is True
        assert result.conservation_score is not None
        assert 0.0 <= result.conservation_score <= 1.0
        assert result.conservation_method == "clustal_omega"

    def test_m4_cache_hit(self, m1_completed_record, tmp_path):
        """Cached scores are loaded without running the pipeline."""
        from varis.m4_conservation import run
        # Prepare a cache file
        cache_dir = tmp_path / "conservation"
        cache_dir.mkdir()
        cache_data = {
            "method": "clustal_omega",
            "num_orthologs": 45,
            "msa_num_sequences": 46,
            "positions": {
                "1699": {
                    "conservation_score": 0.95,
                    "position_entropy": 0.22,
                    "msa_column_index": 1699,
                    "msa_gap_fraction_at_site": 0.02,
                    "conserved_across_mammals": True,
                }
            },
        }
        cache_file = cache_dir / "P38398_scores.json"
        with open(cache_file, "w") as f:
            json.dump(cache_data, f)

        with patch("varis.m4_conservation._CACHE_DIR", cache_dir):
            result = run(m1_completed_record)

        assert "M4" in result.modules_completed
        assert result.conservation_score == pytest.approx(0.95)
        assert result.conservation_method == "clustal_omega"
        assert result.num_orthologs == 45
        assert result.msa_num_sequences == 46
        assert result.position_entropy == pytest.approx(0.22)

    def test_m4_cache_save(self, m1_completed_record, tmp_path):
        """After successful scoring, cache is written."""
        from varis.m4_conservation import run
        # Build sequence with R at position 1699 (0-indexed: 1698)
        seq = "M" * 1698 + "R" + "M" * (1863 - 1699)
        orthologs = {
            "sequences": {
                "query": seq,
                **{f"orth{i}": seq for i in range(15)},
            },
            "query_id": "query",
            "taxonomy": {f"orth{i}": 9606 + i for i in range(15)},
        }
        cache_dir = tmp_path / "conservation"

        with patch("varis.m4_conservation._CACHE_DIR", cache_dir):
            with patch("varis.m4_conservation.uniprot_orthologs.fetch_orthologs",
                        return_value=(m1_completed_record, orthologs)):
                with patch("varis.m4_conservation.clustal_client.run_alignment",
                            return_value=(m1_completed_record, orthologs)):
                    result = run(m1_completed_record)

        assert "M4" in result.modules_completed
        cache_file = cache_dir / "P38398_scores.json"
        assert cache_file.exists()
        with open(cache_file) as f:
            saved = json.load(f)
        assert "1699" in saved["positions"]
        assert saved["positions"]["1699"]["conservation_score"] is not None
