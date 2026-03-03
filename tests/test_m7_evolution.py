"""Tests for M7: Self-Evolution — Model Archive and Evolution Log."""

import json
import pytest
from pathlib import Path

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
