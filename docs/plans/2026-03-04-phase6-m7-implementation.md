# Phase 6: M7 Model Lifecycle Management Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement model version archiving, evolution logging, deploy gate with threshold-based governance, auto-retrain loop, and rollback — so Varis can safely retrain monthly and maintain an auditable history of every model change.

**Architecture:** Three new files in `m7_evolution/`: `model_archive.py` (version storage, atomic deploy, rollback), `evolution_log.py` (SQLAlchemy event log), `auto_retrain.py` (retrain loop with deploy gate). CLI via `__main__.py`. Reuses existing `m5_scoring/train.py` for actual training. Filesystem lock for concurrency.

**Tech Stack:** SQLAlchemy (evolution_log table), pathlib (archive management), hashlib (file hashes), os/signal (PID-based lock), argparse (CLI)

---

### Task 1: Model Archive — Save, List, Deploy, Rollback

**Files:**
- Create: `varis/m7_evolution/model_archive.py`
- Create: `tests/test_m7_evolution.py`

**Step 1: Write failing tests**

Create `tests/test_m7_evolution.py`:

```python
"""Tests for M7: Model Lifecycle Management."""
import json
import os
import pytest
from pathlib import Path


class TestModelArchive:
    """Tests for model_archive.py — version storage and deployment."""

    @pytest.fixture
    def archive_dir(self, tmp_path):
        """Create a temporary archive directory structure."""
        (tmp_path / "archive").mkdir()
        (tmp_path / "candidates").mkdir()
        return tmp_path

    def _create_fake_model(self, model_dir, version="v2026.03"):
        """Create fake model files for testing."""
        model_dir.mkdir(parents=True, exist_ok=True)
        for fname in ["catboost_model.cbm", "xgboost_model.json",
                       "lightgbm_model.txt", "calibrator.pkl", "feature_columns.json"]:
            (model_dir / fname).write_text(f"fake {fname} content for {version}")

    def test_archive_version(self, archive_dir):
        """Save a model version to the archive."""
        from varis.m7_evolution.model_archive import archive_version
        src = archive_dir / "source_model"
        self._create_fake_model(src)
        metrics = {"roc_auc": 0.85, "pr_auc": 0.82}
        result = archive_version(src, "v2026.03", metrics, archive_root=archive_dir)
        assert result is True
        assert (archive_dir / "archive" / "v2026.03").exists()
        meta_path = archive_dir / "archive" / "v2026.03" / "version_metadata.json"
        assert meta_path.exists()
        meta = json.loads(meta_path.read_text())
        assert meta["version"] == "v2026.03"
        assert meta["status"] == "archived"
        assert "file_hashes" in meta

    def test_list_versions(self, archive_dir):
        """List all archived versions."""
        from varis.m7_evolution.model_archive import archive_version, list_versions
        src = archive_dir / "src"
        self._create_fake_model(src, "v2026.03")
        archive_version(src, "v2026.03", {"roc_auc": 0.85}, archive_root=archive_dir)
        self._create_fake_model(src, "v2026.04")
        archive_version(src, "v2026.04", {"roc_auc": 0.87}, archive_root=archive_dir)
        versions = list_versions(archive_root=archive_dir)
        assert len(versions) >= 2

    def test_file_hashes(self, archive_dir):
        """SHA256 hashes stored correctly."""
        from varis.m7_evolution.model_archive import archive_version
        src = archive_dir / "src"
        self._create_fake_model(src)
        archive_version(src, "v2026.03", {}, archive_root=archive_dir)
        meta = json.loads(
            (archive_dir / "archive" / "v2026.03" / "version_metadata.json").read_text()
        )
        assert all(h.startswith("sha256:") for h in meta["file_hashes"].values())

    def test_atomic_deploy(self, archive_dir):
        """Deploy creates symlink current → archive version."""
        from varis.m7_evolution.model_archive import archive_version, deploy_version
        src = archive_dir / "src"
        self._create_fake_model(src)
        archive_version(src, "v2026.03", {}, archive_root=archive_dir)
        deploy_version("v2026.03", archive_root=archive_dir)
        current = archive_dir / "current"
        assert current.is_symlink() or current.is_dir()
        meta = json.loads(
            (archive_dir / "archive" / "v2026.03" / "version_metadata.json").read_text()
        )
        assert meta["status"] == "production"

    def test_rollback(self, archive_dir):
        """Rollback restores previous production version."""
        from varis.m7_evolution.model_archive import (
            archive_version, deploy_version, rollback, get_current_version,
        )
        src = archive_dir / "src"
        self._create_fake_model(src, "v2026.03")
        archive_version(src, "v2026.03", {"roc_auc": 0.85}, archive_root=archive_dir)
        deploy_version("v2026.03", archive_root=archive_dir)
        self._create_fake_model(src, "v2026.04")
        archive_version(src, "v2026.04", {"roc_auc": 0.87}, archive_root=archive_dir)
        deploy_version("v2026.04", archive_root=archive_dir)
        assert get_current_version(archive_root=archive_dir) == "v2026.04"
        rollback("testing rollback", archive_root=archive_dir)
        assert get_current_version(archive_root=archive_dir) == "v2026.03"

    def test_retention_cleanup(self, archive_dir):
        """Old candidates pruned to last 5."""
        from varis.m7_evolution.model_archive import archive_version, cleanup_old_versions
        src = archive_dir / "src"
        for i in range(8):
            version = f"v2026.{i:02d}_candidate"
            self._create_fake_model(src, version)
            cand_dir = archive_dir / "candidates" / version
            cand_dir.mkdir(parents=True, exist_ok=True)
            for f in src.iterdir():
                (cand_dir / f.name).write_bytes(f.read_bytes())
        cleanup_old_versions(archive_root=archive_dir, max_candidates=5)
        remaining = list((archive_dir / "candidates").iterdir())
        assert len(remaining) <= 5
```

**Step 2: Implement model_archive.py**

Key functions:
- `archive_version(source_dir, version, metrics, archive_root)` — copy model files to `archive/{version}/`, compute SHA256 hashes, write `version_metadata.json`
- `deploy_version(version, archive_root)` — validate hashes, create/update `current` symlink, update metadata status to "production"
- `rollback(reason, archive_root)` — find previous production version, validate hashes, symlink flip, update metadata
- `get_current_version(archive_root)` — read symlink target or scan metadata
- `list_versions(archive_root)` — list all archived versions with metadata
- `mark_rejected(version, reason, archive_root)` — update status to "rejected"
- `cleanup_old_versions(archive_root, max_candidates=5, max_months=12)` — prune old candidates and archived versions
- `_compute_file_hashes(directory)` — SHA256 per model file
- `_get_reproducibility_info()` — git commit, python version, library versions

**Step 3: Run tests, commit**

```
feat(m7): implement model archive with atomic deploy and rollback

Version storage with SHA256 hashes, symlink-based atomic deploy,
rollback to previous production, retention cleanup.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 2: Evolution Log — SQLAlchemy Event Persistence

**Files:**
- Replace: `varis/m7_evolution/evolution_log.py`
- Modify: `tests/test_m7_evolution.py`

**Step 1: Write failing tests**

Add to `tests/test_m7_evolution.py`:

```python
class TestEvolutionLog:
    """Tests for evolution_log.py — event persistence."""

    @pytest.fixture
    def log_db(self, tmp_path):
        """Create a temporary database for logging."""
        from varis.m7_evolution.evolution_log import init_evolution_log
        db_url = f"sqlite:///{tmp_path}/test_evo.db"
        session_factory = init_evolution_log(db_url)
        return session_factory

    def test_log_event(self, log_db):
        """Event persists with correct fields."""
        from varis.m7_evolution.evolution_log import log_event, get_log
        log_event(log_db, "DEPLOY", model_version="v2026.03",
                  details={"roc_auc": 0.87})
        events = get_log(log_db)
        assert len(events) == 1
        assert events[0]["event_type"] == "DEPLOY"
        assert events[0]["model_version"] == "v2026.03"

    def test_get_log_filtered(self, log_db):
        """Filter by event_type works."""
        from varis.m7_evolution.evolution_log import log_event, get_log
        log_event(log_db, "DEPLOY", model_version="v2026.03")
        log_event(log_db, "REJECT", model_version="v2026.04")
        log_event(log_db, "DEPLOY", model_version="v2026.05")
        deploys = get_log(log_db, event_type="DEPLOY")
        assert len(deploys) == 2
        rejects = get_log(log_db, event_type="REJECT")
        assert len(rejects) == 1

    def test_log_event_with_details(self, log_db):
        """Details stored as JSON and retrievable."""
        from varis.m7_evolution.evolution_log import log_event, get_log
        details = {"old_roc_auc": 0.85, "new_roc_auc": 0.87, "delta": 0.02}
        log_event(log_db, "DEPLOY", model_version="v2026.03", details=details)
        events = get_log(log_db)
        assert events[0]["details"]["delta"] == 0.02
```

**Step 2: Implement evolution_log.py**

- SQLAlchemy model `EvolutionLogRow` with columns from design
- `init_evolution_log(database_url)` → session factory (creates table if needed)
- `log_event(session_factory, event_type, model_version=None, details=None)` — insert row
- `get_log(session_factory, limit=50, event_type=None)` → list of dicts
- Uses the same `DATABASE_URL` config as M6 when called from the API

**Step 3: Run tests, commit**

```
feat(m7): implement evolution log with SQLAlchemy persistence

Event types: RETRAIN_START, RETRAIN_COMPLETE, DEPLOY, REJECT,
ROLLBACK, ERROR. Filterable by type. Details stored as JSON.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 3: Deploy Gate — Threshold-Based Governance

**Files:**
- Modify: `tests/test_m7_evolution.py`
- Create logic in `varis/m7_evolution/auto_retrain.py` (deploy gate portion)

**Step 1: Write failing tests**

```python
class TestDeployGate:
    """Tests for deploy gate — threshold-based governance."""

    def test_deploy_gate_passes(self):
        """All metrics improve → DEPLOY."""
        from varis.m7_evolution.auto_retrain import evaluate_deploy_gate
        current = {"roc_auc": 0.85, "pr_auc": 0.82, "calibration_ece": 0.05}
        candidate = {"roc_auc": 0.87, "pr_auc": 0.84, "calibration_ece": 0.03}
        current_benchmarks = {"BRCA1_p.Arg1699Trp": 0.90}
        candidate_benchmarks = {"BRCA1_p.Arg1699Trp": 0.91}
        result = evaluate_deploy_gate(current, candidate, current_benchmarks, candidate_benchmarks)
        assert result["decision"] == "DEPLOY"

    def test_deploy_gate_regression(self):
        """ROC-AUC drops > 0.005 → REJECT."""
        from varis.m7_evolution.auto_retrain import evaluate_deploy_gate
        current = {"roc_auc": 0.87, "pr_auc": 0.84, "calibration_ece": 0.03}
        candidate = {"roc_auc": 0.86, "pr_auc": 0.85, "calibration_ece": 0.03}
        result = evaluate_deploy_gate(current, candidate, {}, {})
        assert result["decision"] == "REJECT"
        assert "roc_auc" in result["reason"]

    def test_deploy_gate_no_improvement(self):
        """No meaningful improvement → REJECT."""
        from varis.m7_evolution.auto_retrain import evaluate_deploy_gate
        current = {"roc_auc": 0.85, "pr_auc": 0.84, "calibration_ece": 0.03}
        candidate = {"roc_auc": 0.85, "pr_auc": 0.844, "calibration_ece": 0.029}
        result = evaluate_deploy_gate(current, candidate, {}, {})
        assert result["decision"] == "REJECT"
        assert "no meaningful improvement" in result["reason"].lower()

    def test_deploy_gate_benchmark_flip(self):
        """Classification flip likely_benign → likely_pathogenic → REJECT."""
        from varis.m7_evolution.auto_retrain import evaluate_deploy_gate
        current = {"roc_auc": 0.85, "pr_auc": 0.82, "calibration_ece": 0.05}
        candidate = {"roc_auc": 0.87, "pr_auc": 0.84, "calibration_ece": 0.03}
        # Score flip: was 0.1 (benign), now 0.9 (pathogenic)
        current_bm = {"BRCA1_p.Lys1183Arg": 0.10}
        candidate_bm = {"BRCA1_p.Lys1183Arg": 0.90}
        result = evaluate_deploy_gate(current, candidate, current_bm, candidate_bm)
        assert result["decision"] == "REJECT"
        assert "classification flip" in result["reason"].lower()

    def test_deploy_gate_score_drift(self):
        """Benchmark score drift > 0.15 → REJECT."""
        from varis.m7_evolution.auto_retrain import evaluate_deploy_gate
        current = {"roc_auc": 0.85, "pr_auc": 0.82, "calibration_ece": 0.05}
        candidate = {"roc_auc": 0.87, "pr_auc": 0.84, "calibration_ece": 0.03}
        current_bm = {"BRCA1_p.Arg1699Trp": 0.85}
        candidate_bm = {"BRCA1_p.Arg1699Trp": 0.60}  # Drift of 0.25
        result = evaluate_deploy_gate(current, candidate, current_bm, candidate_bm)
        assert result["decision"] == "REJECT"
        assert "score drift" in result["reason"].lower()
```

**Step 2: Implement `evaluate_deploy_gate()` in auto_retrain.py**

Pure function: takes current metrics, candidate metrics, current benchmark scores, candidate benchmark scores. Returns `{"decision": "DEPLOY"/"REJECT", "reason": str, "metric_deltas": dict}`.

Hard gates: regression thresholds (0.005 AUC, 0.01 ECE), benchmark flip, score drift >0.15.
Soft gate: PR-AUC improves ≥0.01 or ECE improves ≥0.01.

Classification from score: >0.8 pathogenic, <0.2 benign, else uncertain.

**Step 3: Run tests, commit**

```
feat(m7): implement threshold-based deploy gate

Hard gates: no regression >0.005, no benchmark flip, no score drift >0.15.
Soft gate: PR-AUC or ECE must improve ≥0.01. Logs rejection reason.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 4: Concurrency Lock

**Files:**
- Create logic in `varis/m7_evolution/auto_retrain.py` (lock portion)
- Modify: `tests/test_m7_evolution.py`

**Step 1: Write failing tests**

```python
class TestConcurrencyLock:
    """Tests for filesystem-based concurrency lock."""

    def test_acquire_and_release(self, tmp_path):
        from varis.m7_evolution.auto_retrain import acquire_lock, release_lock
        lock_path = tmp_path / ".retrain.lock"
        assert acquire_lock(lock_path) is True
        assert lock_path.exists()
        release_lock(lock_path)
        assert not lock_path.exists()

    def test_lock_blocks_second_run(self, tmp_path):
        from varis.m7_evolution.auto_retrain import acquire_lock
        lock_path = tmp_path / ".retrain.lock"
        assert acquire_lock(lock_path) is True
        # Second acquire should fail (same PID is alive)
        assert acquire_lock(lock_path) is False

    def test_stale_lock_cleaned(self, tmp_path):
        from varis.m7_evolution.auto_retrain import acquire_lock
        lock_path = tmp_path / ".retrain.lock"
        # Write a lock with a dead PID
        lock_path.write_text(json.dumps({"pid": 99999999, "timestamp": "2026-01-01T00:00:00"}))
        # Should detect stale lock and acquire
        assert acquire_lock(lock_path) is True
```

**Step 2: Implement lock functions**

- `acquire_lock(lock_path)` → bool: write PID+timestamp, check if existing lock's PID is alive via `os.kill(pid, 0)`
- `release_lock(lock_path)` → None: remove lock file
- Use in retrain loop with `try/finally`

**Step 3: Run tests, commit**

```
feat(m7): implement filesystem concurrency lock

PID-based lock prevents concurrent retrain runs. Detects and
cleans stale locks from dead processes.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 5: Auto-Retrain Loop

**Files:**
- Replace: `varis/m7_evolution/auto_retrain.py` (complete the implementation)
- Modify: `tests/test_m7_evolution.py`

**Step 1: Write failing tests**

```python
from unittest.mock import patch, MagicMock

class TestRetrainLoop:
    """Tests for the full auto-retrain loop."""

    def test_retrain_loop_deploy(self, tmp_path):
        """Mocked: candidate improves → DEPLOY."""
        from varis.m7_evolution.auto_retrain import run_retrain_loop
        from varis.m7_evolution.evolution_log import init_evolution_log
        log_db = init_evolution_log(f"sqlite:///{tmp_path}/evo.db")
        # Mock the training pipeline to return good metrics
        with patch("varis.m7_evolution.auto_retrain.select_training_variants") as mock_select, \
             patch("varis.m7_evolution.auto_retrain.compute_features") as mock_compute, \
             patch("varis.m7_evolution.auto_retrain.train_and_evaluate") as mock_train, \
             patch("varis.m7_evolution.auto_retrain.run_benchmarks") as mock_bench:
            mock_select.return_value = MagicMock()
            mock_compute.return_value = (100, 50, 50, 0)
            mock_train.return_value = {
                "cv_results": {"roc_auc_mean": 0.90, "pr_auc_mean": 0.88},
            }
            mock_bench.return_value = {}
            result = run_retrain_loop(
                archive_root=tmp_path, log_db=log_db,
                lock_path=tmp_path / ".lock",
            )
        assert result["completed"] is True

    def test_retrain_loop_locked(self, tmp_path):
        """Lock held → exits without training."""
        from varis.m7_evolution.auto_retrain import run_retrain_loop, acquire_lock
        from varis.m7_evolution.evolution_log import init_evolution_log
        log_db = init_evolution_log(f"sqlite:///{tmp_path}/evo.db")
        lock_path = tmp_path / ".lock"
        acquire_lock(lock_path)
        result = run_retrain_loop(
            archive_root=tmp_path, log_db=log_db, lock_path=lock_path,
        )
        assert result["completed"] is False
        assert "lock" in result.get("reason", "").lower()
```

**Step 2: Implement `run_retrain_loop()`**

Full orchestration: lock → log start → train (reuse train.py functions) → benchmark → deploy gate → deploy/reject → cleanup → unlock → log complete.

Imports from `varis.m5_scoring.train`: `select_training_variants`, `compute_features`, `load_cached_records`, `train_and_evaluate`.

**Step 3: Run tests, commit**

```
feat(m7): implement full auto-retrain loop

Lock → select → compute → train → benchmark → gate → deploy/reject →
cleanup → unlock. Reuses m5_scoring/train.py pipeline.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 6: CLI Entry Point

**Files:**
- Create: `varis/m7_evolution/__main__.py`
- Modify: `varis/m7_evolution/__init__.py`

**Step 1: Implement __main__.py**

CLI with subcommands:
```
python -m varis.m7_evolution retrain
python -m varis.m7_evolution rollback --reason "text"
python -m varis.m7_evolution status
python -m varis.m7_evolution log --limit 20 --type DEPLOY
```

Uses argparse. Each subcommand calls the appropriate function from auto_retrain.py, model_archive.py, or evolution_log.py.

**Step 2: Update __init__.py** with proper imports.

**Step 3: Verify CLI runs**

```bash
python -m varis.m7_evolution status
python -m varis.m7_evolution log --limit 5
```

**Step 4: Commit**

```
feat(m7): implement CLI for retrain, rollback, status, and log

python -m varis.m7_evolution retrain|rollback|status|log
Scheduled via system cron for hands-off monthly retraining.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 7: API Endpoint for Evolution Log

**Files:**
- Modify: `varis/m6_platform/api/main.py`
- Modify: `tests/test_m6_api.py`

**Step 1: Add evolution-log endpoint to FastAPI**

```python
@app.get("/api/v1/evolution-log")
def get_evolution_log(limit: int = 50, event_type: str | None = None):
    """Retrieve evolution log entries."""
    from varis.m7_evolution.evolution_log import get_log, init_evolution_log
    log_db = init_evolution_log(DATABASE_URL)
    events = get_log(log_db, limit=limit, event_type=event_type)
    return {"events": events}
```

**Step 2: Add test**

```python
def test_evolution_log_endpoint(self, client):
    resp = client.get("/api/v1/evolution-log")
    assert resp.status_code == 200
    assert "events" in resp.json()
```

**Step 3: Commit**

```
feat(m7): add evolution-log API endpoint

GET /api/v1/evolution-log with optional limit and event_type filters.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 8: Bootstrap Current Model into Archive

**Files:**
- Create: `varis/m7_evolution/bootstrap.py` (one-time script)

**Step 1: Implement bootstrap script**

Takes the existing trained model at `data/models/` and archives it as the first production version:

```python
"""Bootstrap the existing trained model into the version archive.

Run once after initial training to set up the archive structure.
Usage: python -m varis.m7_evolution.bootstrap
"""
```

- Reads existing model files from `data/models/`
- Archives as `v2026.03` (current date-based version)
- Creates `current` symlink
- Logs DEPLOY event to evolution log
- Updates M5 `load_ensemble()` to load from `data/models/current/` if it exists, falling back to `data/models/` for backward compatibility

**Step 2: Run bootstrap, verify**

**Step 3: Commit**

```
feat(m7): add bootstrap script to archive initial trained model

One-time migration: moves existing model into version archive,
creates current symlink, logs initial deploy event.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 9: Final Verification and Cleanup

**Step 1: Run full test suite**

```
pytest tests/ -v --tb=short
```

**Step 2: Verify CLI commands work**

```
python -m varis.m7_evolution status
python -m varis.m7_evolution log --limit 5
```

**Step 3: Search for stale references**

Old event type constants, old function names from stubs.

**Step 4: Commit cleanup**

```
chore: Phase 6 cleanup — verify M7 tests and remove stale references

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```
