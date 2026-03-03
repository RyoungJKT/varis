# Phase 6 Design: M7 — Model Lifecycle Management

**Date:** 2026-03-04
**Status:** Approved
**Depends on:** M5 training pipeline, M6 API
**Scope:** Priorities 1-2 only (Evolution Log + Auto-Retrain Loop)

## Decisions Made

- **Priorities 1-2 only**: Evolution log + auto-retrain with governance. Tool scout and auto-integrator deferred.
- **SQLite table** in existing varis.db for evolution log
- **CLI + system cron** for scheduling (no in-process scheduler)
- **Relaxed deploy gate**: threshold-based, not "all must improve"
- **Atomic deploy** via symlink flip with file hash validation
- **Filesystem lock** for concurrency protection

## Model Version Archive

### Directory Structure

```
data/models/
├── current → archive/v2026.03/    (symlink to production version)
├── archive/
│   ├── v2026.03/
│   │   ├── catboost_model.cbm
│   │   ├── xgboost_model.json
│   │   ├── lightgbm_model.txt
│   │   ├── calibrator.pkl
│   │   ├── feature_columns.json
│   │   ├── training_metadata.json
│   │   └── version_metadata.json
│   ├── v2026.04/
│   └── ...
└── candidates/
    └── v2026.04_candidate/
```

### version_metadata.json

```json
{
  "version": "v2026.04",
  "status": "production",
  "deployed_at": "2026-04-01T02:00:00Z",
  "archived_at": null,
  "metrics": {
    "roc_auc": 0.87,
    "pr_auc": 0.84,
    "calibration_ece": 0.03,
    "brier_score": 0.12
  },
  "benchmark_scores": {
    "BRCA1_p.Arg1699Trp": 0.91,
    "BRCA1_p.Lys1183Arg": 0.08
  },
  "training_samples": 580,
  "clinvar_export_date": "2026-04-01",
  "file_hashes": {
    "catboost_model.cbm": "sha256:abc123...",
    "xgboost_model.json": "sha256:def456..."
  },
  "reproducibility": {
    "git_commit": "abc123",
    "python_version": "3.11.8",
    "catboost_version": "1.2.3",
    "xgboost_version": "2.0.1",
    "lightgbm_version": "4.1.0",
    "shap_version": "0.43.0",
    "feature_schema_version": "1.4.0",
    "random_seed": 42
  },
  "rejection_reason": null,
  "rollback_reason": null
}
```

### Retention Policy

- Keep last 5 candidates (reject or not)
- Keep 12 months of production versions
- `_cleanup_old_versions()` runs after each retrain cycle

## Evolution Log

### Database Table (SQLAlchemy, in existing varis.db)

```
evolution_log:
  id              — Integer, primary key
  event_type      — String: RETRAIN_START, RETRAIN_COMPLETE, DEPLOY,
                            REJECT, ROLLBACK, ERROR
  model_version   — String, nullable
  details         — Text (JSON: metrics, deltas, reasons, duration)
  timestamp       — DateTime, indexed
```

### API Endpoint

`GET /api/v1/evolution-log?limit=50&event_type=DEPLOY`

Returns list of events with timestamp, type, version, and details.

## Deploy Gate (Threshold-Based)

The strict "ALL metrics must improve" rule from CLAUDE.md is relaxed to be
practically usable while maintaining safety:

### Hard Gates (reject if any fail)

1. **No meaningful regression**: ROC-AUC drop ≤ 0.005 AND PR-AUC drop ≤ 0.005
   AND calibration ECE increase ≤ 0.01
2. **Benchmark stability**: No classification flip between likely_benign ↔
   likely_pathogenic on any benchmark variant
3. **Benchmark score drift**: No benchmark variant score changes by > 0.15

### Soft Gate (require at least one)

4. **Meaningful improvement**: PR-AUC improves by ≥ 0.01 OR calibration ECE
   improves by ≥ 0.01

### Decision Logic

```
if any hard gate fails → REJECT (log which gate and metric deltas)
if no soft gate passes → REJECT (log "no meaningful improvement")
else → DEPLOY
```

Movement within "uncertain" classification on benchmarks is allowed.
ROC-AUC is a guardrail; PR-AUC is the primary metric (better for imbalanced data).

## Auto-Retrain Loop

CLI command: `python -m varis.m7_evolution retrain`

### Flow

```
1. Acquire filesystem lock (data/models/.retrain.lock)
   - If lock exists and process alive → exit, log ERROR
   - If lock exists but stale → warn, proceed
2. Log RETRAIN_START
3. Download latest ClinVar (reuse train.py Phase A)
4. Select new variants — compare with previous manifest, find additions
5. Compute features for new variants (reuse train.py Phase B)
6. Train candidate model (reuse train.py Phase C)
   - Save to data/models/candidates/v{YYYY.MM}_candidate/
7. Run regression tests against benchmark variants
8. Load current production metrics from version_metadata.json
9. Compare candidate vs production metrics (deploy gate)
10. If DEPLOY:
    a. Write candidate to temp dir, validate file hashes
    b. Move current production to archive (update status)
    c. Atomic symlink flip: current → new version
    d. Update version_metadata.json status to "production"
    e. Log DEPLOY with metric deltas
11. If REJECT:
    a. Archive candidate as "rejected" with reason
    b. Log REJECT with metric deltas and which gate failed
12. Cleanup old candidates (keep last 5)
13. Release filesystem lock
14. Log RETRAIN_COMPLETE with duration
```

### Rollback

CLI: `python -m varis.m7_evolution rollback --reason "regression in BRCA1"`

1. Find most recent archived production version (status != "rejected")
2. Validate file hashes match stored metadata
3. Atomic symlink flip: current → previous version
4. Update metadata: old version → "rolled_back", previous → "production"
5. Log ROLLBACK with reason

## Concurrency Protection

Filesystem lock at `data/models/.retrain.lock`:
- Contains PID + timestamp
- On start: check if lock exists → check if PID is alive
  - Alive → exit with ERROR log
  - Dead → warn "stale lock", remove, proceed
- On exit (success or failure): remove lock
- Use `try/finally` to ensure lock removal

## Files

| File | Responsibility |
|------|---------------|
| `evolution_log.py` | Log events to DB, query events, get current model version |
| `auto_retrain.py` | Full retrain loop: lock → train → gate → deploy/reject → unlock |
| `model_archive.py` | Archive management: save/load/list versions, deploy (atomic), rollback, cleanup |
| `__init__.py` | Thin orchestrator |
| `__main__.py` | CLI: `retrain`, `rollback`, `status`, `log` subcommands |

## CLI

```
python -m varis.m7_evolution retrain              # Run full retrain loop
python -m varis.m7_evolution rollback --reason "text"  # Rollback to previous
python -m varis.m7_evolution status               # Show current model + recent history
python -m varis.m7_evolution log --limit 20       # Show evolution log
python -m varis.m7_evolution log --type DEPLOY    # Filter by event type
```

## Testing Strategy

| Test | What it verifies |
|------|-----------------|
| `test_log_event` | Event persists to DB with correct fields |
| `test_get_log_filtered` | Filter by event_type works |
| `test_archive_save_and_list` | Save version, list archived versions |
| `test_file_hashes` | SHA256 computed and stored correctly |
| `test_deploy_gate_passes` | Metrics improve → DEPLOY |
| `test_deploy_gate_regression` | Metric drops > threshold → REJECT |
| `test_deploy_gate_no_improvement` | No meaningful improvement → REJECT |
| `test_deploy_gate_benchmark_flip` | Classification flip → REJECT |
| `test_deploy_gate_score_drift` | Benchmark score drift > 0.15 → REJECT |
| `test_atomic_deploy` | Symlink flip, old version archived |
| `test_rollback` | Previous version restored, log updated |
| `test_concurrency_lock` | Second run exits when lock held |
| `test_stale_lock` | Dead PID lock detected and cleaned |
| `test_retention_cleanup` | Old candidates pruned to last 5 |
| `test_retrain_loop_mocked` | Full loop with mocked training |

## Naming

Internal module: `m7_evolution` (unchanged)
User-facing/docs: "Model Lifecycle Management" or "Scheduled Retraining"
API section: `/api/v1/evolution-log` (keep for URL stability)
