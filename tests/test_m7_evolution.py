"""Tests for M7: Self-Evolution — Model Archive, Evolution Log, Deploy Gate, Lock, Retrain Loop, and CLI."""

import json
import os
import pytest
from pathlib import Path
from unittest.mock import patch, MagicMock

from varis.m7_evolution.auto_retrain import (
    evaluate_deploy_gate,
    acquire_lock,
    release_lock,
    run_retrain_loop,
)
from varis.m7_evolution.__main__ import build_parser
from varis.m7_evolution.model_archive import (
    archive_version,
    deploy_version,
    rollback,
    get_current_version,
    list_versions,
    mark_rejected,
    cleanup_old_versions,
    _compute_file_hashes,
)


def _create_fake_model_dir(path: Path) -> Path:
    """Create a directory with fake model files for testing."""
    path.mkdir(parents=True, exist_ok=True)
    (path / "model.cbm").write_text("fake catboost model")
    (path / "model.json").write_text("fake xgboost model")
    (path / "model.txt").write_text("fake lightgbm model")
    (path / "model.pkl").write_text("fake sklearn model")
    return path


class TestModelArchive:
    """Tests for model version archive, deploy, rollback, and cleanup."""

    def test_archive_version(self, tmp_path: Path) -> None:
        """Archive a version: verify dir exists and metadata has correct fields."""
        source = _create_fake_model_dir(tmp_path / "source")
        archive_root = tmp_path / "archive_root"

        result = archive_version(source, "v2026.03", {"roc_auc": 0.92}, archive_root)

        assert result is True
        version_dir = archive_root / "archive" / "v2026.03"
        assert version_dir.is_dir()
        assert (version_dir / "model.cbm").exists()
        assert (version_dir / "model.json").exists()
        assert (version_dir / "model.txt").exists()
        assert (version_dir / "model.pkl").exists()

        # Verify metadata
        meta_path = version_dir / "version_metadata.json"
        assert meta_path.exists()
        with open(meta_path, "r") as f:
            metadata = json.load(f)

        assert metadata["version"] == "v2026.03"
        assert metadata["status"] == "archived"
        assert metadata["metrics"]["roc_auc"] == 0.92
        assert "file_hashes" in metadata
        assert "reproducibility" in metadata
        assert "archived_at" in metadata
        assert "python_version" in metadata["reproducibility"]

    def test_list_versions(self, tmp_path: Path) -> None:
        """Archive 2 versions, list returns both."""
        archive_root = tmp_path / "archive_root"

        source1 = _create_fake_model_dir(tmp_path / "source1")
        source2 = _create_fake_model_dir(tmp_path / "source2")

        archive_version(source1, "v2026.01", {"roc_auc": 0.88}, archive_root)
        archive_version(source2, "v2026.02", {"roc_auc": 0.91}, archive_root)

        versions = list_versions(archive_root)

        assert len(versions) == 2
        version_names = [v["version"] for v in versions]
        assert "v2026.01" in version_names
        assert "v2026.02" in version_names

    def test_file_hashes(self, tmp_path: Path) -> None:
        """Verify SHA256 hashes are stored with 'sha256:' prefix."""
        source = _create_fake_model_dir(tmp_path / "source")
        archive_root = tmp_path / "archive_root"

        archive_version(source, "v2026.03", {"roc_auc": 0.92}, archive_root)

        meta_path = archive_root / "archive" / "v2026.03" / "version_metadata.json"
        with open(meta_path, "r") as f:
            metadata = json.load(f)

        file_hashes = metadata["file_hashes"]
        assert len(file_hashes) == 4  # .cbm, .json, .txt, .pkl

        for filename, hash_value in file_hashes.items():
            assert hash_value.startswith("sha256:"), f"Hash for {filename} missing prefix"
            hex_part = hash_value.split(":", 1)[1]
            assert len(hex_part) == 64, f"Hash for {filename} has wrong length"

        # Verify hashes match direct computation
        version_dir = archive_root / "archive" / "v2026.03"
        direct_hashes = _compute_file_hashes(version_dir)
        for filename, hash_value in file_hashes.items():
            assert direct_hashes[filename] == hash_value

    def test_atomic_deploy(self, tmp_path: Path) -> None:
        """Deploy creates symlink (or dir on Windows) and metadata shows 'production'."""
        source = _create_fake_model_dir(tmp_path / "source")
        archive_root = tmp_path / "archive_root"

        archive_version(source, "v2026.03", {"roc_auc": 0.92}, archive_root)
        result = deploy_version("v2026.03", archive_root)

        assert result is True

        # Check current link exists
        current_link = archive_root / "current"
        assert current_link.exists()

        # Check metadata updated to production
        meta_path = archive_root / "archive" / "v2026.03" / "version_metadata.json"
        with open(meta_path, "r") as f:
            metadata = json.load(f)
        assert metadata["status"] == "production"
        assert "deployed_at" in metadata

        # Check get_current_version
        assert get_current_version(archive_root) == "v2026.03"

    def test_rollback(self, tmp_path: Path) -> None:
        """Deploy v1, deploy v2, rollback -> current is v1 again."""
        archive_root = tmp_path / "archive_root"

        source1 = _create_fake_model_dir(tmp_path / "source1")
        source2 = _create_fake_model_dir(tmp_path / "source2")

        # Archive and deploy v1
        archive_version(source1, "v2026.01", {"roc_auc": 0.88}, archive_root)
        deploy_version("v2026.01", archive_root)
        assert get_current_version(archive_root) == "v2026.01"

        # Archive and deploy v2 (v1 becomes archived)
        archive_version(source2, "v2026.02", {"roc_auc": 0.91}, archive_root)
        deploy_version("v2026.02", archive_root)
        assert get_current_version(archive_root) == "v2026.02"

        # Rollback to v1
        result = rollback("Regression detected in production", archive_root)
        assert result is True
        assert get_current_version(archive_root) == "v2026.01"

        # Verify v2 is now rolled_back
        meta_path = archive_root / "archive" / "v2026.02" / "version_metadata.json"
        with open(meta_path, "r") as f:
            metadata = json.load(f)
        assert metadata["status"] == "rolled_back"
        assert metadata["rollback_reason"] == "Regression detected in production"

    def test_retention_cleanup(self, tmp_path: Path) -> None:
        """Create 8 candidates, cleanup to 5."""
        archive_root = tmp_path / "archive_root"

        # Create 8 archived versions
        for i in range(1, 9):
            source = _create_fake_model_dir(tmp_path / f"source_{i}")
            version = f"v2026.{i:02d}"
            archive_version(source, version, {"roc_auc": 0.80 + i * 0.01}, archive_root)

        # Verify all 8 exist
        versions_before = list_versions(archive_root)
        assert len(versions_before) == 8

        # Cleanup to max_candidates=5
        cleanup_old_versions(archive_root, max_candidates=5)

        # Verify only 5 remain (the newest 5)
        versions_after = list_versions(archive_root)
        assert len(versions_after) == 5

        remaining_names = sorted([v["version"] for v in versions_after])
        assert remaining_names == ["v2026.04", "v2026.05", "v2026.06", "v2026.07", "v2026.08"]

    def test_mark_rejected(self, tmp_path: Path) -> None:
        """Mark a version as rejected, verify status and reason."""
        source = _create_fake_model_dir(tmp_path / "source")
        archive_root = tmp_path / "archive_root"

        archive_version(source, "v2026.03", {"roc_auc": 0.70}, archive_root)
        mark_rejected("v2026.03", "Metrics below threshold", archive_root)

        versions = list_versions(archive_root)
        assert len(versions) == 1
        assert versions[0]["status"] == "rejected"
        assert versions[0]["rejection_reason"] == "Metrics below threshold"

    def test_deploy_nonexistent_version(self, tmp_path: Path) -> None:
        """Deploying a version that doesn't exist returns False."""
        archive_root = tmp_path / "archive_root"
        (archive_root / "archive").mkdir(parents=True)

        result = deploy_version("v9999.01", archive_root)
        assert result is False

    def test_rollback_no_previous(self, tmp_path: Path) -> None:
        """Rollback with only one version returns False."""
        archive_root = tmp_path / "archive_root"

        source = _create_fake_model_dir(tmp_path / "source")
        archive_version(source, "v2026.01", {"roc_auc": 0.88}, archive_root)
        deploy_version("v2026.01", archive_root)

        result = rollback("Test reason", archive_root)
        assert result is False

    def test_archive_nonexistent_source(self, tmp_path: Path) -> None:
        """Archiving from a nonexistent directory returns False."""
        archive_root = tmp_path / "archive_root"
        result = archive_version(
            tmp_path / "nonexistent", "v2026.01", {"roc_auc": 0.5}, archive_root,
        )
        assert result is False

    def test_cleanup_preserves_production(self, tmp_path: Path) -> None:
        """Cleanup does not remove the production version."""
        archive_root = tmp_path / "archive_root"

        # Create 3 versions, deploy the middle one
        for i in range(1, 4):
            source = _create_fake_model_dir(tmp_path / f"source_{i}")
            archive_version(source, f"v2026.{i:02d}", {"roc_auc": 0.80 + i * 0.01}, archive_root)

        deploy_version("v2026.02", archive_root)

        # Cleanup with max_candidates=0 — should remove non-production but keep production
        cleanup_old_versions(archive_root, max_candidates=0)

        versions = list_versions(archive_root)
        version_names = [v["version"] for v in versions]
        assert "v2026.02" in version_names  # Production preserved


class TestEvolutionLog:
    """Tests for SQLAlchemy-based evolution log persistence."""

    @pytest.fixture
    def log_db(self, tmp_path):
        """Create a temporary SQLite-backed evolution log session factory."""
        from varis.m7_evolution.evolution_log import init_evolution_log

        db_url = f"sqlite:///{tmp_path}/test_evo.db"
        return init_evolution_log(db_url)

    def test_log_event(self, log_db) -> None:
        """Log a DEPLOY event and verify it appears in get_log."""
        from varis.m7_evolution.evolution_log import get_log, log_event

        log_event(log_db, "DEPLOY", model_version="v2026.03", details={"roc_auc": 0.87})
        events = get_log(log_db)
        assert len(events) == 1
        assert events[0]["event_type"] == "DEPLOY"
        assert events[0]["model_version"] == "v2026.03"

    def test_get_log_filtered(self, log_db) -> None:
        """Filter events by type: 2 DEPLOYs and 1 REJECT."""
        from varis.m7_evolution.evolution_log import get_log, log_event

        log_event(log_db, "DEPLOY", model_version="v2026.03")
        log_event(log_db, "REJECT", model_version="v2026.04")
        log_event(log_db, "DEPLOY", model_version="v2026.05")
        assert len(get_log(log_db, event_type="DEPLOY")) == 2
        assert len(get_log(log_db, event_type="REJECT")) == 1

    def test_log_event_with_details(self, log_db) -> None:
        """Verify details dict is round-tripped through JSON serialization."""
        from varis.m7_evolution.evolution_log import get_log, log_event

        details = {"old_roc_auc": 0.85, "new_roc_auc": 0.87, "delta": 0.02}
        log_event(log_db, "DEPLOY", model_version="v2026.03", details=details)
        events = get_log(log_db)
        assert events[0]["details"]["delta"] == 0.02

    def test_get_current_model_version(self, log_db) -> None:
        """Most recent DEPLOY event determines current model version."""
        from varis.m7_evolution.evolution_log import get_current_model_version, log_event

        log_event(log_db, "DEPLOY", model_version="v2026.03")
        log_event(log_db, "DEPLOY", model_version="v2026.04")
        assert get_current_model_version(log_db) == "v2026.04"

    def test_get_current_model_version_no_deploys(self, log_db) -> None:
        """No DEPLOY events returns None."""
        from varis.m7_evolution.evolution_log import get_current_model_version, log_event

        log_event(log_db, "REJECT", model_version="v2026.03")
        assert get_current_model_version(log_db) is None

    def test_log_event_returns_dict(self, log_db) -> None:
        """log_event returns a dict with all expected keys."""
        from varis.m7_evolution.evolution_log import log_event

        result = log_event(log_db, "RETRAIN_START", model_version="v2026.03")
        assert isinstance(result, dict)
        assert "id" in result
        assert result["event_type"] == "RETRAIN_START"
        assert result["model_version"] == "v2026.03"
        assert result["timestamp"] is not None

    def test_get_log_limit(self, log_db) -> None:
        """Verify get_log respects the limit parameter."""
        from varis.m7_evolution.evolution_log import get_log, log_event

        for i in range(10):
            log_event(log_db, "ERROR", details={"index": i})
        events = get_log(log_db, limit=3)
        assert len(events) == 3

    def test_get_log_ordering(self, log_db) -> None:
        """Events are returned newest-first (descending timestamp)."""
        from varis.m7_evolution.evolution_log import get_log, log_event

        log_event(log_db, "DEPLOY", model_version="v2026.01")
        log_event(log_db, "DEPLOY", model_version="v2026.02")
        log_event(log_db, "DEPLOY", model_version="v2026.03")
        events = get_log(log_db)
        assert events[0]["model_version"] == "v2026.03"
        assert events[2]["model_version"] == "v2026.01"


class TestDeployGate:
    """Tests for evaluate_deploy_gate threshold logic."""

    def _base_metrics(self) -> dict:
        """Return a baseline metrics dict."""
        return {"roc_auc": 0.90, "pr_auc": 0.85, "calibration_ece": 0.05}

    def _base_benchmarks(self) -> dict:
        """Return a baseline benchmarks dict."""
        return {"BRCA1_R1699W": 0.92, "TP53_R175H": 0.10, "MLH1_V384D": 0.50}

    def test_deploy_gate_passes(self) -> None:
        """All metrics improve -> DEPLOY."""
        current = self._base_metrics()
        candidate = {"roc_auc": 0.92, "pr_auc": 0.87, "calibration_ece": 0.03}
        benchmarks = self._base_benchmarks()

        result = evaluate_deploy_gate(current, candidate, benchmarks, benchmarks)

        assert result["decision"] == "DEPLOY"
        assert result["metric_deltas"]["roc_auc"] > 0
        assert result["metric_deltas"]["pr_auc"] > 0

    def test_deploy_gate_regression(self) -> None:
        """ROC-AUC drops >0.005 -> REJECT with 'roc_auc' in reason."""
        current = self._base_metrics()
        candidate = {"roc_auc": 0.89, "pr_auc": 0.87, "calibration_ece": 0.03}
        benchmarks = self._base_benchmarks()

        result = evaluate_deploy_gate(current, candidate, benchmarks, benchmarks)

        assert result["decision"] == "REJECT"
        assert "roc_auc" in result["reason"]

    def test_deploy_gate_no_improvement(self) -> None:
        """No metric improves >=0.01 -> REJECT with 'no meaningful improvement'."""
        current = self._base_metrics()
        # Tiny improvements, not enough to pass soft gate
        candidate = {"roc_auc": 0.905, "pr_auc": 0.855, "calibration_ece": 0.045}
        benchmarks = self._base_benchmarks()

        result = evaluate_deploy_gate(current, candidate, benchmarks, benchmarks)

        assert result["decision"] == "REJECT"
        assert "no meaningful improvement" in result["reason"]

    def test_deploy_gate_benchmark_flip(self) -> None:
        """Score goes 0.10->0.90 (benign->pathogenic) -> REJECT with 'classification flip'."""
        current = self._base_metrics()
        candidate = {"roc_auc": 0.92, "pr_auc": 0.87, "calibration_ece": 0.03}
        current_benchmarks = {"BRCA1_R1699W": 0.92, "TP53_R175H": 0.10}
        candidate_benchmarks = {"BRCA1_R1699W": 0.92, "TP53_R175H": 0.90}

        result = evaluate_deploy_gate(
            current, candidate, current_benchmarks, candidate_benchmarks,
        )

        assert result["decision"] == "REJECT"
        assert "classification flip" in result["reason"]

    def test_deploy_gate_score_drift(self) -> None:
        """Benchmark score drifts 0.25 -> REJECT with 'score drift'."""
        current = self._base_metrics()
        candidate = {"roc_auc": 0.92, "pr_auc": 0.87, "calibration_ece": 0.03}
        current_benchmarks = {"MLH1_V384D": 0.50}
        # Drift of 0.25 exceeds the 0.15 limit
        candidate_benchmarks = {"MLH1_V384D": 0.75}

        result = evaluate_deploy_gate(
            current, candidate, current_benchmarks, candidate_benchmarks,
        )

        assert result["decision"] == "REJECT"
        assert "score drift" in result["reason"]

    def test_deploy_gate_ece_improvement_only(self) -> None:
        """PR-AUC flat but ECE improves >=0.01 -> DEPLOY (soft gate passes)."""
        current = self._base_metrics()
        # PR-AUC unchanged, but ECE improves by 0.02 (lower is better)
        candidate = {"roc_auc": 0.90, "pr_auc": 0.855, "calibration_ece": 0.03}
        benchmarks = self._base_benchmarks()

        result = evaluate_deploy_gate(current, candidate, benchmarks, benchmarks)

        assert result["decision"] == "DEPLOY"

    def test_deploy_gate_within_uncertain(self) -> None:
        """Benchmark score changes within uncertain range (0.3->0.5) -> allowed."""
        current = self._base_metrics()
        candidate = {"roc_auc": 0.92, "pr_auc": 0.87, "calibration_ece": 0.03}
        current_benchmarks = {"MLH1_V384D": 0.30}
        # Drift of 0.10 is within the 0.15 limit, and both are "uncertain"
        candidate_benchmarks = {"MLH1_V384D": 0.40}

        result = evaluate_deploy_gate(
            current, candidate, current_benchmarks, candidate_benchmarks,
        )

        # No classification flip (both uncertain), drift within limit
        assert result["decision"] == "DEPLOY"


class TestConcurrencyLock:
    """Tests for filesystem-based PID concurrency lock."""

    def test_acquire_and_release(self, tmp_path: Path) -> None:
        """Acquire succeeds, file exists, release removes it."""
        lock_file = tmp_path / "retrain.lock"

        assert acquire_lock(lock_file) is True
        assert lock_file.exists()

        # Verify lock content
        lock_data = json.loads(lock_file.read_text())
        assert lock_data["pid"] == os.getpid()
        assert "timestamp" in lock_data

        release_lock(lock_file)
        assert not lock_file.exists()

    def test_lock_blocks_second_run(self, tmp_path: Path) -> None:
        """Acquire twice with same PID -> second returns False."""
        lock_file = tmp_path / "retrain.lock"

        assert acquire_lock(lock_file) is True
        # Second acquire should fail because our own PID is still alive
        assert acquire_lock(lock_file) is False

        release_lock(lock_file)

    def test_stale_lock_cleaned(self, tmp_path: Path) -> None:
        """Write lock with PID 99999999, acquire detects stale and succeeds."""
        lock_file = tmp_path / "retrain.lock"

        # Write a lock file with a PID that almost certainly doesn't exist
        stale_lock = {
            "pid": 99999999,
            "timestamp": "2025-01-01T00:00:00+00:00",
        }
        lock_file.write_text(json.dumps(stale_lock))

        # Acquire should detect the stale lock and succeed
        assert acquire_lock(lock_file) is True

        # Verify our PID is now in the lock
        lock_data = json.loads(lock_file.read_text())
        assert lock_data["pid"] == os.getpid()

        release_lock(lock_file)


class TestRetrainLoop:
    """Tests for the full auto-retrain orchestration loop."""

    def test_retrain_loop_mocked(self, tmp_path: Path) -> None:
        """Full loop with mocked training -- completes."""
        from varis.m7_evolution.evolution_log import init_evolution_log
        from varis.m7_evolution.model_archive import archive_version, deploy_version

        # Set up: archive a "current" model first
        src = tmp_path / "src"
        src.mkdir()
        for f in [
            "catboost_model.cbm",
            "xgboost_model.json",
            "lightgbm_model.txt",
            "calibrator.pkl",
            "feature_columns.json",
        ]:
            (src / f).write_text("fake")
        archive_version(
            src,
            "v2026.03",
            {"roc_auc": 0.85, "pr_auc": 0.82, "calibration_ece": 0.05},
            archive_root=tmp_path,
        )
        deploy_version("v2026.03", archive_root=tmp_path)

        log_db = init_evolution_log(f"sqlite:///{tmp_path}/evo.db")

        # Mock the training pipeline
        with patch("varis.m7_evolution.auto_retrain.select_training_variants") as mock_sel, \
             patch("varis.m7_evolution.auto_retrain.compute_features") as mock_comp, \
             patch("varis.m7_evolution.auto_retrain.load_cached_records") as mock_load, \
             patch("varis.m7_evolution.auto_retrain.train_and_evaluate") as mock_train, \
             patch("varis.m7_evolution.auto_retrain.run_benchmarks_for_model") as mock_bench:

            mock_sel.return_value = MagicMock()
            mock_comp.return_value = (100, 50, 50, 0)
            mock_load.return_value = (
                [],
                [1] * 50 + [0] * 50,
                ["BRCA1"] * 50 + ["TP53"] * 50,
            )
            # Return value that won't trip "not records" check:
            # The load returns empty records list, but we need at least one
            # record for the "if not records" check to pass.
            mock_load.return_value = (
                [MagicMock()],
                [1] * 50 + [0] * 50,
                ["BRCA1"] * 50 + ["TP53"] * 50,
            )
            mock_train.return_value = {
                "cv_results": {
                    "roc_auc_mean": 0.90,
                    "pr_auc_mean": 0.88,
                    "pr_auc_std": 0.02,
                    "roc_auc_std": 0.02,
                },
                "model_version": "v2026.04",
            }

            # Create fake candidate model files
            cand_dir = tmp_path / "candidates" / "v2026.04_candidate"
            cand_dir.mkdir(parents=True)
            for f in [
                "catboost_model.cbm",
                "xgboost_model.json",
                "lightgbm_model.txt",
                "calibrator.pkl",
                "feature_columns.json",
            ]:
                (cand_dir / f).write_text("fake candidate")

            mock_bench.return_value = {"BRCA1_p.Arg1699Trp": 0.91}

            result = run_retrain_loop(
                archive_root=tmp_path,
                log_db=log_db,
                lock_path=tmp_path / ".lock",
            )

        assert result["completed"] is True
        assert result["decision"] in ("DEPLOY", "REJECT")
        assert result["duration"] > 0.0

    def test_retrain_loop_locked(self, tmp_path: Path) -> None:
        """Lock held -- exits without training."""
        from varis.m7_evolution.evolution_log import init_evolution_log

        log_db = init_evolution_log(f"sqlite:///{tmp_path}/evo.db")
        lock_path = tmp_path / ".lock"
        acquire_lock(lock_path)

        result = run_retrain_loop(
            archive_root=tmp_path,
            log_db=log_db,
            lock_path=lock_path,
        )

        assert result["completed"] is False
        assert "lock" in result.get("reason", "").lower()

        # Clean up the lock for other tests
        release_lock(lock_path)


class TestCLI:
    """Tests for the M7 CLI entry point (__main__.py)."""

    def test_build_parser(self) -> None:
        """Parser builds successfully with all subcommands."""
        from varis.m7_evolution.__main__ import build_parser

        parser = build_parser()
        assert parser is not None

        # Verify retrain subcommand parses
        args = parser.parse_args(["retrain"])
        assert args.command == "retrain"

        # Verify rollback subcommand parses
        args = parser.parse_args(["rollback", "--reason", "test"])
        assert args.command == "rollback"
        assert args.reason == "test"

        # Verify status subcommand parses
        args = parser.parse_args(["status"])
        assert args.command == "status"

        # Verify log subcommand parses with defaults
        args = parser.parse_args(["log"])
        assert args.command == "log"
        assert args.limit == 20
        assert args.type is None

        # Verify log subcommand parses with options
        args = parser.parse_args(["log", "--limit", "5", "--type", "DEPLOY"])
        assert args.limit == 5
        assert args.type == "DEPLOY"

    def test_verbose_flag(self) -> None:
        """Verbose flag is parsed correctly."""
        from varis.m7_evolution.__main__ import build_parser

        parser = build_parser()
        args = parser.parse_args(["-v", "status"])
        assert args.verbose is True

        args = parser.parse_args(["status"])
        assert args.verbose is False

    def test_cmd_status(self, tmp_path: Path, capsys) -> None:
        """Status command prints version info."""
        from varis.m7_evolution.__main__ import cmd_status, build_parser

        # Create an archive with one version
        source = _create_fake_model_dir(tmp_path / "source")
        archive_root = tmp_path / "archive_root"
        archive_version(source, "v2026.03", {"roc_auc": 0.92}, archive_root)
        deploy_version("v2026.03", archive_root)

        parser = build_parser()
        args = parser.parse_args(["--archive-root", str(archive_root), "status"])
        cmd_status(args)

        captured = capsys.readouterr()
        assert "v2026.03" in captured.out
        assert "production" in captured.out

    def test_cmd_log(self, tmp_path: Path, capsys) -> None:
        """Log command prints events."""
        from varis.m7_evolution.__main__ import cmd_log, build_parser
        from varis.m7_evolution.evolution_log import init_evolution_log, log_event

        db_url = f"sqlite:///{tmp_path}/test_cli.db"
        log_db = init_evolution_log(db_url)
        log_event(log_db, "DEPLOY", model_version="v2026.03", details={"roc_auc": 0.92})

        parser = build_parser()
        args = parser.parse_args(["log", "--limit", "10"])

        with patch("varis.m7_evolution.__main__._get_log_db", return_value=log_db):
            cmd_log(args)

        captured = capsys.readouterr()
        assert "DEPLOY" in captured.out
        assert "v2026.03" in captured.out

    def test_cmd_rollback(self, tmp_path: Path, capsys) -> None:
        """Rollback command succeeds with two deployed versions."""
        from varis.m7_evolution.__main__ import cmd_rollback, build_parser

        archive_root = tmp_path / "archive_root"

        source1 = _create_fake_model_dir(tmp_path / "source1")
        source2 = _create_fake_model_dir(tmp_path / "source2")

        archive_version(source1, "v2026.01", {"roc_auc": 0.88}, archive_root)
        deploy_version("v2026.01", archive_root)
        archive_version(source2, "v2026.02", {"roc_auc": 0.91}, archive_root)
        deploy_version("v2026.02", archive_root)

        parser = build_parser()
        args = parser.parse_args([
            "--archive-root", str(archive_root),
            "rollback", "--reason", "test rollback",
        ])
        cmd_rollback(args)

        captured = capsys.readouterr()
        assert "Rollback successful" in captured.out

        # Verify we're back to v2026.01
        assert get_current_version(archive_root) == "v2026.01"

    def test_scout_command_parses(self):
        """CLI parser accepts 'scout' subcommand."""
        parser = build_parser()
        args = parser.parse_args(["scout"])
        assert args.command == "scout"

    def test_integrate_command_parses(self):
        """CLI parser accepts 'integrate --tool X' subcommand."""
        parser = build_parser()
        args = parser.parse_args(["integrate", "--tool", "some-package"])
        assert args.command == "integrate"
        assert args.tool == "some-package"


class TestModuleImports:
    """Verify M7 modules can be imported and have expected public API."""

    def test_tool_scout_api(self):
        """tool_scout has all expected public functions."""
        from varis.m7_evolution import tool_scout
        assert callable(tool_scout.scan_sources)
        assert callable(tool_scout.score_candidate)
        assert callable(tool_scout.deduplicate)
        assert callable(tool_scout.run_scout_loop)

    def test_auto_integrator_api(self):
        """auto_integrator has all expected public functions."""
        from varis.m7_evolution import auto_integrator
        assert callable(auto_integrator.attempt_integration)
        assert callable(auto_integrator.attempt_install)
        assert callable(auto_integrator.probe_tool)
        assert callable(auto_integrator.benchmark_new_feature)


class TestM7Run:
    """Tests for the M7 module-level run() function."""

    def test_run_function_exists(self) -> None:
        """M7 has a run() function for pipeline integration."""
        from varis.m7_evolution import run
        assert callable(run)

    def test_run_returns_dict_on_error(self) -> None:
        """run() returns dict with completed=False on import error."""
        from varis.m7_evolution import run

        with patch(
            "varis.m7_evolution.auto_retrain.run_retrain_loop",
            side_effect=RuntimeError("test error"),
        ):
            result = run()

        assert isinstance(result, dict)
        assert result["completed"] is False


class TestToolScout:
    """Tests for Tool Scout (Loop 2) — discovery, scoring, and deduplication."""

    # --- Scoring tests ---

    def test_score_candidate_high_relevance(self) -> None:
        """Candidate with DDG/SASA/conservation keywords scores >= 3."""
        from varis.m7_evolution.tool_scout import score_candidate

        candidate = {
            "name": "ddg-predictor",
            "description": "Predicts delta-g and SASA for protein stability analysis with conservation scores",
            "keywords": ["ddg", "sasa", "conservation"],
            "stars": 50,
            "updated": "2020-01-01",
        }
        score = score_candidate(candidate)
        assert score >= 3

    def test_score_candidate_low_relevance(self) -> None:
        """Candidate with unrelated keywords (web/flask) scores < 3."""
        from varis.m7_evolution.tool_scout import score_candidate

        candidate = {
            "name": "flask-admin",
            "description": "Admin panel for Flask web applications",
            "keywords": ["web", "flask", "admin"],
            "stars": 5000,
            "updated": "2020-01-01",
        }
        score = score_candidate(candidate)
        assert score < 3

    def test_score_candidate_recency_bonus(self) -> None:
        """Recent candidate scores higher than an old one with same keywords."""
        from datetime import datetime, timedelta
        from varis.m7_evolution.tool_scout import score_candidate

        recent_date = (datetime.now() - timedelta(days=30)).strftime("%Y-%m-%d")
        old_date = "2020-01-01"

        base = {
            "name": "stability-tool",
            "description": "Protein stability predictor",
            "keywords": ["stability"],
            "stars": 50,
        }

        recent_candidate = {**base, "updated": recent_date}
        old_candidate = {**base, "updated": old_date}

        assert score_candidate(recent_candidate) > score_candidate(old_candidate)

    # --- Scanner tests ---

    def test_scan_pypi_returns_candidates(self) -> None:
        """Mock httpx, verify scan_pypi returns list of candidate dicts."""
        from varis.m7_evolution.tool_scout import scan_pypi

        # First response: HTML search page with package links
        mock_search_response = MagicMock()
        mock_search_response.status_code = 200
        mock_search_response.text = '<a href="/project/bio-stability/">bio-stability</a>'

        # Second response: JSON API for package details
        mock_detail_response = MagicMock()
        mock_detail_response.status_code = 200
        mock_detail_response.json.return_value = {
            "info": {
                "name": "bio-stability",
                "summary": "Protein stability tool",
                "keywords": "protein,stability",
                "project_url": "https://pypi.org/project/bio-stability/",
                "version": "1.0.0",
            },
            "releases": {
                "1.0.0": [{"upload_time_iso_8601": "2026-01-15T00:00:00Z"}],
            },
            "urls": [],
        }

        def side_effect(url, **kwargs):
            if "/search" in url:
                return mock_search_response
            return mock_detail_response

        with patch("varis.m7_evolution.tool_scout.httpx.get", side_effect=side_effect):
            results = scan_pypi("protein stability")

        assert isinstance(results, list)
        assert len(results) >= 1
        assert results[0]["name"] == "bio-stability"
        assert results[0]["source"] == "pypi"

    def test_scan_github_returns_candidates(self) -> None:
        """Mock httpx, verify scan_github returns candidate with stars/source."""
        from varis.m7_evolution.tool_scout import scan_github

        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "items": [
                {
                    "name": "variant-predictor",
                    "full_name": "user/variant-predictor",
                    "html_url": "https://github.com/user/variant-predictor",
                    "description": "Variant effect predictor tool",
                    "stargazers_count": 250,
                    "topics": ["bioinformatics", "variant"],
                    "updated_at": "2026-02-01T00:00:00Z",
                },
            ],
        }

        with patch("varis.m7_evolution.tool_scout.httpx.get", return_value=mock_response):
            results = scan_github("variant effect")

        assert isinstance(results, list)
        assert len(results) >= 1
        assert results[0]["source"] == "github"
        assert results[0]["stars"] == 250
        assert results[0]["name"] == "variant-predictor"

    def test_scan_sources_aggregates(self) -> None:
        """Mock all three scanners, verify aggregation from all sources."""
        from varis.m7_evolution.tool_scout import scan_sources, SCOUT_SOURCES

        pypi_result = [{"name": "tool-a", "source": "pypi", "description": "A"}]
        github_result = [{"name": "tool-b", "source": "github", "description": "B"}]
        biorxiv_result = [{"name": "tool-c", "source": "biorxiv", "description": "C"}]

        # Each scanner is called once per keyword in SCOUT_SOURCES
        n_pypi = len(SCOUT_SOURCES["pypi"]["keywords"])
        n_github = len(SCOUT_SOURCES["github"]["keywords"])
        n_biorxiv = len(SCOUT_SOURCES["biorxiv"]["keywords"])
        expected_total = n_pypi + n_github + n_biorxiv

        with patch("varis.m7_evolution.tool_scout.scan_pypi", return_value=pypi_result), \
             patch("varis.m7_evolution.tool_scout.scan_github", return_value=github_result), \
             patch("varis.m7_evolution.tool_scout.scan_biorxiv", return_value=biorxiv_result):
            results = scan_sources()

        assert isinstance(results, list)
        assert len(results) == expected_total
        # All three source types should be represented
        sources = {r["source"] for r in results}
        assert sources == {"pypi", "github", "biorxiv"}

    # --- Orchestrator tests ---

    def test_run_scout_loop(self, tmp_path: Path) -> None:
        """Mock scan_sources, verify scoring + logging to evolution log."""
        from varis.m7_evolution.evolution_log import init_evolution_log, get_log
        from varis.m7_evolution.tool_scout import run_scout_loop

        db_url = f"sqlite:///{tmp_path}/scout_test.db"
        log_db = init_evolution_log(db_url)

        candidates = [
            {
                "name": "ddg-analyzer",
                "source": "pypi",
                "url": "https://pypi.org/project/ddg-analyzer/",
                "description": "Delta-G stability prediction with conservation and SASA analysis",
                "keywords": ["ddg", "stability", "conservation", "sasa"],
                "stars": 200,
                "updated": "2026-02-15",
            },
            {
                "name": "flask-login",
                "source": "pypi",
                "url": "https://pypi.org/project/flask-login/",
                "description": "User authentication for Flask",
                "keywords": ["web", "flask"],
                "stars": 3000,
                "updated": "2026-01-01",
            },
        ]

        with patch("varis.m7_evolution.tool_scout.scan_sources", return_value=candidates):
            result = run_scout_loop(log_db=log_db)

        assert result["scanned"] == 2
        # Only the high-scoring ddg-analyzer should be logged
        assert result["logged"] >= 1
        # Verify it was actually written to the evolution log
        events = get_log(log_db, event_type="TOOL_DISCOVERY")
        assert len(events) >= 1
        assert any("ddg-analyzer" in str(e.get("details", "")) for e in events)

    def test_run_scout_loop_deduplicates(self, tmp_path: Path) -> None:
        """Run twice with same candidates, verify only logged once."""
        from varis.m7_evolution.evolution_log import init_evolution_log, get_log
        from varis.m7_evolution.tool_scout import run_scout_loop

        db_url = f"sqlite:///{tmp_path}/scout_dedup.db"
        log_db = init_evolution_log(db_url)

        candidates = [
            {
                "name": "stability-predictor",
                "source": "github",
                "url": "https://github.com/user/stability-predictor",
                "description": "Protein stability and delta-g prediction with conservation",
                "keywords": ["stability", "ddg", "conservation"],
                "stars": 150,
                "updated": "2026-02-01",
            },
        ]

        with patch("varis.m7_evolution.tool_scout.scan_sources", return_value=candidates):
            result1 = run_scout_loop(log_db=log_db)
            result2 = run_scout_loop(log_db=log_db)

        # First run should log it
        assert result1["logged"] >= 1
        # Second run should deduplicate and not log again
        assert result2["logged"] == 0

        # Verify only one TOOL_DISCOVERY event in the DB
        events = get_log(log_db, event_type="TOOL_DISCOVERY")
        stability_events = [
            e for e in events if "stability-predictor" in str(e.get("details", ""))
        ]
        assert len(stability_events) == 1


class TestAutoIntegrator:
    """Tests for Auto-Integrator (Loop 3) — install, probe, benchmark, orchestrator."""

    # --- Install tests ---

    def test_attempt_install_success(self) -> None:
        """Mock subprocess.run returning 0, verify True and 2 calls (dry-run + install)."""
        from varis.m7_evolution.auto_integrator import attempt_install

        mock_result = MagicMock()
        mock_result.returncode = 0

        with patch("varis.m7_evolution.auto_integrator.subprocess.run", return_value=mock_result) as mock_run:
            result = attempt_install("some-package")

        assert result is True
        assert mock_run.call_count == 2  # dry-run + actual install

    def test_attempt_install_failure(self) -> None:
        """Mock subprocess.run returning 1 on dry-run, verify False."""
        from varis.m7_evolution.auto_integrator import attempt_install

        mock_result = MagicMock()
        mock_result.returncode = 1

        with patch("varis.m7_evolution.auto_integrator.subprocess.run", return_value=mock_result):
            result = attempt_install("nonexistent-package")

        assert result is False

    # --- Probe tests ---

    def test_probe_tool_finds_callable(self) -> None:
        """Mock importlib.import_module with fake module that has `predict`, verify info dict."""
        from varis.m7_evolution.auto_integrator import probe_tool

        fake_module = MagicMock()
        fake_module.__dir__ = lambda self: ["predict", "SomeClass", "__version__"]
        fake_module.predict = lambda x: x  # callable
        fake_module.SomeClass = type("SomeClass", (), {})
        fake_module.__version__ = "1.0.0"

        with patch("varis.m7_evolution.auto_integrator.importlib.import_module", return_value=fake_module):
            info = probe_tool("some_package")

        assert info is not None
        assert info["package"] == "some_package"
        assert "predict" in info["callables"]
        assert "module_attrs" in info

    def test_probe_tool_import_fails(self) -> None:
        """Mock ImportError, verify None is returned."""
        from varis.m7_evolution.auto_integrator import probe_tool

        with patch(
            "varis.m7_evolution.auto_integrator.importlib.import_module",
            side_effect=ImportError("No module named 'nonexistent'"),
        ):
            info = probe_tool("nonexistent")

        assert info is None

    # --- Benchmark tests ---

    def test_benchmark_new_feature_improves(self) -> None:
        """Current roc_auc=0.850, candidate=0.860, verify INTEGRATE."""
        from varis.m7_evolution.auto_integrator import benchmark_new_feature

        current = {"roc_auc": 0.850, "pr_auc": 0.820}
        candidate = {"roc_auc": 0.860, "pr_auc": 0.825}

        result = benchmark_new_feature(current, candidate)

        assert result["decision"] == "INTEGRATE"
        assert "metric_deltas" in result
        assert result["metric_deltas"]["roc_auc"] == pytest.approx(0.010, abs=1e-6)

    def test_benchmark_new_feature_regresses(self) -> None:
        """Candidate worse, verify REJECT."""
        from varis.m7_evolution.auto_integrator import benchmark_new_feature

        current = {"roc_auc": 0.850, "pr_auc": 0.820}
        candidate = {"roc_auc": 0.840, "pr_auc": 0.810}

        result = benchmark_new_feature(current, candidate)

        assert result["decision"] == "REJECT"
        assert "metric_deltas" in result

    # --- Orchestrator tests ---

    def test_attempt_integration_success(self, tmp_path: Path) -> None:
        """Mock install+probe+benchmark, verify INTEGRATE + evolution log event."""
        from varis.m7_evolution.auto_integrator import attempt_integration
        from varis.m7_evolution.evolution_log import init_evolution_log, get_log

        db_url = f"sqlite:///{tmp_path}/integrator_test.db"
        log_db = init_evolution_log(db_url)

        proposal = {"name": "test-tool", "package": "test_tool"}
        current_metrics = {"roc_auc": 0.850, "pr_auc": 0.820}
        candidate_metrics = {"roc_auc": 0.860, "pr_auc": 0.825}

        with patch("varis.m7_evolution.auto_integrator.attempt_install", return_value=True), \
             patch("varis.m7_evolution.auto_integrator.probe_tool", return_value={
                 "package": "test_tool",
                 "callables": ["predict"],
                 "module_attrs": 5,
             }):
            result = attempt_integration(
                proposal,
                log_db=log_db,
                current_metrics=current_metrics,
                candidate_metrics=candidate_metrics,
            )

        assert result["decision"] == "INTEGRATE"
        assert result["name"] == "test-tool"

        # Verify evolution log event was written
        events = get_log(log_db, event_type="TOOL_INTEGRATION")
        assert len(events) >= 1
        assert any("test-tool" in str(e.get("details", "")) for e in events)

    def test_attempt_integration_install_fails(self, tmp_path: Path) -> None:
        """Mock install failure, verify REJECT with 'install' in reason."""
        from varis.m7_evolution.auto_integrator import attempt_integration
        from varis.m7_evolution.evolution_log import init_evolution_log

        db_url = f"sqlite:///{tmp_path}/integrator_fail.db"
        log_db = init_evolution_log(db_url)

        proposal = {"name": "bad-tool", "package": "bad_tool"}

        with patch("varis.m7_evolution.auto_integrator.attempt_install", return_value=False):
            result = attempt_integration(proposal, log_db=log_db)

        assert result["decision"] == "REJECT"
        assert "install" in result["reason"].lower()
