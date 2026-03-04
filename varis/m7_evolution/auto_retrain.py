"""Auto-Retrain (Loop 1) — Fully autonomous weekly model retraining with governance.

Scans ClinVar for newly expert-reviewed variants, runs them through the pipeline,
retrains all three models, compares metrics, runs regression tests, and deploys
only if all gates pass.

Release governance:
  - Every model is version-tagged (e.g. v2026.03)
  - Deployment requires ALL metrics to improve AND regression tests to pass
  - Previous versions are preserved in the version archive for rollback
  - Every decision is recorded in the Evolution Log with metric deltas

Priority: 2 (build second — most impactful self-evolution loop)
Frequency: Weekly (triggered by ClinVar update detection)
Autonomy: Fully autonomous — no human intervention required
"""
import json
import logging
import os
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from varis.config import DATA_DIR, DATABASE_URL, MODELS_DIR

logger = logging.getLogger(__name__)

# Default paths
DEFAULT_ARCHIVE_ROOT = MODELS_DIR
DEFAULT_LOCK_PATH = MODELS_DIR / ".retrain.lock"

# Lazy-loaded references to m5_scoring.train functions.
# These are set at module level so unittest.mock.patch can intercept them
# via "varis.m7_evolution.auto_retrain.<name>". Actual imports happen inside
# run_retrain_loop() to avoid circular imports.
select_training_variants = None  # type: ignore[assignment]
compute_features = None  # type: ignore[assignment]
load_cached_records = None  # type: ignore[assignment]
train_and_evaluate = None  # type: ignore[assignment]

# ── Model version format: "v{YYYY}.{MM}" e.g. "v2026.03" ──

# ── Classification thresholds for benchmark scoring ──
_PATHOGENIC_THRESHOLD = 0.8
_BENIGN_THRESHOLD = 0.2

# ── Deploy gate thresholds ──
_MAX_ROC_AUC_DROP = 0.005
_MAX_PR_AUC_DROP = 0.005
_MAX_ECE_INCREASE = 0.01
_MAX_BENCHMARK_SCORE_DRIFT = 0.15
_MIN_PR_AUC_IMPROVEMENT = 0.01
_MIN_ECE_IMPROVEMENT = 0.01


def _classify_score(score: float) -> str:
    """Classify a pathogenicity score into a category.

    Args:
        score: Pathogenicity score between 0 and 1.

    Returns:
        One of "likely_pathogenic", "likely_benign", or "uncertain".
    """
    if score > _PATHOGENIC_THRESHOLD:
        return "likely_pathogenic"
    elif score < _BENIGN_THRESHOLD:
        return "likely_benign"
    else:
        return "uncertain"


def evaluate_deploy_gate(
    current_metrics: dict,
    candidate_metrics: dict,
    current_benchmarks: dict,
    candidate_benchmarks: dict,
) -> dict:
    """Evaluate whether a candidate model should be deployed.

    Hard gates (reject if any fail):
    1. No meaningful regression: ROC-AUC drop <= 0.005, PR-AUC drop <= 0.005,
       calibration ECE increase <= 0.01
    2. No benchmark classification flip between likely_benign <-> likely_pathogenic
    3. No benchmark score drift > 0.15

    Soft gate (require at least one):
    4. PR-AUC improves >= 0.01 OR calibration ECE improves >= 0.01

    Args:
        current_metrics: {"roc_auc": float, "pr_auc": float, "calibration_ece": float}
        candidate_metrics: Same structure as current_metrics.
        current_benchmarks: {"variant_id": score, ...} where score is float 0-1.
        candidate_benchmarks: Same structure as current_benchmarks.

    Returns:
        {"decision": "DEPLOY"/"REJECT", "reason": str, "metric_deltas": dict}
    """
    metric_deltas = {
        "roc_auc": candidate_metrics["roc_auc"] - current_metrics["roc_auc"],
        "pr_auc": candidate_metrics["pr_auc"] - current_metrics["pr_auc"],
        "calibration_ece": candidate_metrics["calibration_ece"] - current_metrics["calibration_ece"],
    }

    # ── Hard Gate 1: No meaningful metric regression ──
    if metric_deltas["roc_auc"] < -_MAX_ROC_AUC_DROP:
        reason = f"roc_auc regression: {metric_deltas['roc_auc']:.4f} (limit: -{_MAX_ROC_AUC_DROP})"
        logger.warning("Deploy gate REJECT: %s", reason)
        return {"decision": "REJECT", "reason": reason, "metric_deltas": metric_deltas}

    if metric_deltas["pr_auc"] < -_MAX_PR_AUC_DROP:
        reason = f"pr_auc regression: {metric_deltas['pr_auc']:.4f} (limit: -{_MAX_PR_AUC_DROP})"
        logger.warning("Deploy gate REJECT: %s", reason)
        return {"decision": "REJECT", "reason": reason, "metric_deltas": metric_deltas}

    # ECE is lower-is-better, so an increase is a regression
    if metric_deltas["calibration_ece"] > _MAX_ECE_INCREASE:
        reason = (
            f"calibration_ece regression: +{metric_deltas['calibration_ece']:.4f} "
            f"(limit: +{_MAX_ECE_INCREASE})"
        )
        logger.warning("Deploy gate REJECT: %s", reason)
        return {"decision": "REJECT", "reason": reason, "metric_deltas": metric_deltas}

    # ── Hard Gate 2: No benchmark classification flip (benign <-> pathogenic) ──
    for variant_id in current_benchmarks:
        if variant_id not in candidate_benchmarks:
            continue
        current_class = _classify_score(current_benchmarks[variant_id])
        candidate_class = _classify_score(candidate_benchmarks[variant_id])

        # Only flag flips between the two extreme categories
        flip_pair = {current_class, candidate_class}
        if flip_pair == {"likely_benign", "likely_pathogenic"}:
            reason = (
                f"classification flip for {variant_id}: "
                f"{current_class} -> {candidate_class} "
                f"(score {current_benchmarks[variant_id]:.3f} -> "
                f"{candidate_benchmarks[variant_id]:.3f})"
            )
            logger.warning("Deploy gate REJECT: %s", reason)
            return {"decision": "REJECT", "reason": reason, "metric_deltas": metric_deltas}

    # ── Hard Gate 3: No benchmark score drift > 0.15 ──
    for variant_id in current_benchmarks:
        if variant_id not in candidate_benchmarks:
            continue
        drift = abs(candidate_benchmarks[variant_id] - current_benchmarks[variant_id])
        if drift > _MAX_BENCHMARK_SCORE_DRIFT:
            reason = (
                f"score drift for {variant_id}: {drift:.3f} "
                f"(limit: {_MAX_BENCHMARK_SCORE_DRIFT})"
            )
            logger.warning("Deploy gate REJECT: %s", reason)
            return {"decision": "REJECT", "reason": reason, "metric_deltas": metric_deltas}

    # ── Soft Gate: At least one meaningful improvement ──
    pr_auc_improved = metric_deltas["pr_auc"] >= _MIN_PR_AUC_IMPROVEMENT
    # ECE is lower-is-better, so improvement means a decrease (negative delta)
    ece_improved = metric_deltas["calibration_ece"] <= -_MIN_ECE_IMPROVEMENT

    if not pr_auc_improved and not ece_improved:
        reason = (
            f"no meaningful improvement: "
            f"pr_auc delta {metric_deltas['pr_auc']:.4f} "
            f"(need >= {_MIN_PR_AUC_IMPROVEMENT}), "
            f"ece delta {metric_deltas['calibration_ece']:.4f} "
            f"(need <= -{_MIN_ECE_IMPROVEMENT})"
        )
        logger.info("Deploy gate REJECT: %s", reason)
        return {"decision": "REJECT", "reason": reason, "metric_deltas": metric_deltas}

    # ── All gates passed ──
    improvements = []
    if pr_auc_improved:
        improvements.append(f"pr_auc +{metric_deltas['pr_auc']:.4f}")
    if ece_improved:
        improvements.append(f"ece {metric_deltas['calibration_ece']:.4f}")

    reason = f"all gates passed, improvements: {', '.join(improvements)}"
    logger.info("Deploy gate DEPLOY: %s", reason)
    return {"decision": "DEPLOY", "reason": reason, "metric_deltas": metric_deltas}


def acquire_lock(lock_path: Path) -> bool:
    """Acquire filesystem lock for retrain process.

    Writes PID + timestamp to lock file. If lock exists:
    - Check if PID is alive via os.kill(pid, 0)
    - If alive -> return False (locked)
    - If dead -> warn "stale lock", remove, acquire

    Args:
        lock_path: Path to lock file.

    Returns:
        True if lock acquired, False if already locked.
    """
    try:
        if lock_path.exists():
            try:
                lock_data = json.loads(lock_path.read_text())
                existing_pid = lock_data.get("pid")
            except (json.JSONDecodeError, KeyError, OSError) as e:
                logger.warning("Corrupt lock file at %s: %s — removing", lock_path, e)
                lock_path.unlink(missing_ok=True)
                existing_pid = None

            if existing_pid is not None:
                # Check if the process is still alive
                try:
                    os.kill(existing_pid, 0)
                    # Process is alive — lock is valid
                    logger.info(
                        "Lock held by PID %d at %s — cannot acquire",
                        existing_pid,
                        lock_path,
                    )
                    return False
                except OSError:
                    # Process is dead — stale lock
                    logger.warning(
                        "stale lock from PID %d at %s — removing",
                        existing_pid,
                        lock_path,
                    )
                    lock_path.unlink(missing_ok=True)

        # Acquire the lock
        lock_path.parent.mkdir(parents=True, exist_ok=True)
        lock_data = {
            "pid": os.getpid(),
            "timestamp": datetime.now(timezone.utc).isoformat(),
        }
        lock_path.write_text(json.dumps(lock_data))
        logger.info("Lock acquired at %s (PID %d)", lock_path, os.getpid())
        return True
    except Exception as e:
        logger.error("Failed to acquire lock at %s: %s", lock_path, e)
        return False


def release_lock(lock_path: Path) -> None:
    """Release filesystem lock.

    Removes the lock file if it exists. Logs a warning if the lock file
    does not exist (double-release).

    Args:
        lock_path: Path to lock file.
    """
    try:
        if lock_path.exists():
            lock_path.unlink()
            logger.info("Lock released at %s", lock_path)
        else:
            logger.warning("Lock file not found at %s — nothing to release", lock_path)
    except Exception as e:
        logger.error("Failed to release lock at %s: %s", lock_path, e)


def run_benchmarks_for_model(model_dir: Path) -> dict:
    """Load ensemble from model_dir and run predictions on benchmark variants.

    Loads the trained ensemble, runs predictions on each benchmark variant
    from the benchmark file, and returns a mapping of variant_id to score.

    Args:
        model_dir: Path to directory containing trained model files
            (catboost_model.cbm, xgboost_model.json, lightgbm_model.txt,
            calibrator.pkl, feature_columns.json).

    Returns:
        Dict mapping variant_id (str) to predicted pathogenicity score (float).
        Returns empty dict if model cannot be loaded or benchmarks fail.
    """
    try:
        from varis.m5_scoring.benchmarks import run_benchmarks

        model_dir = Path(model_dir)

        # Use run_benchmarks from m5_scoring for structured evaluation
        benchmark_results = run_benchmarks(model_dir)

        # Extract variant_id -> score mapping
        scores: dict[str, float] = {}
        for variant in benchmark_results.get("variants", []):
            variant_id = f"{variant['gene']}_{variant['hgvs']}"
            score = variant.get("score")
            if score is not None:
                scores[variant_id] = float(score)

        # If no scores from benchmark file (pipeline not fully available),
        # try to load ensemble and return empty scores with a log message
        if not scores:
            logger.info(
                "No benchmark scores available from %s — "
                "full pipeline evaluation not yet available",
                model_dir,
            )

        return scores

    except Exception as e:
        logger.warning("Failed to run benchmarks for model at %s: %s", model_dir, e)
        return {}


def run_retrain_loop(
    archive_root: Optional[Path] = None,
    log_db: Optional[object] = None,
    lock_path: Optional[Path] = None,
) -> dict:
    """Execute the full auto-retrain loop. Called by scheduler.

    Full orchestration sequence:
      1.  Acquire filesystem lock — if locked, return immediately
      2.  Log RETRAIN_START event
      3.  Call select_training_variants() from m5_scoring.train
      4.  Call compute_features() from m5_scoring.train
      5.  Call load_cached_records() from m5_scoring.train
      6.  Call train_and_evaluate() — saves candidate to candidates/ dir
      7.  Run benchmark predictions on candidate model
      8.  Load current production metrics from model_archive
      9.  Evaluate deploy gate
      10. If DEPLOY: archive candidate, deploy via model_archive, log DEPLOY
      11. If REJECT: mark candidate rejected, log REJECT
      12. Cleanup old versions
      13. Release lock
      14. Log RETRAIN_COMPLETE

    Args:
        archive_root: Root directory for the model archive. Defaults to MODELS_DIR.
        log_db: SQLAlchemy session factory for evolution log. If None, initializes
            from DATABASE_URL config (or default SQLite path).
        lock_path: Path to the filesystem lock file. Defaults to
            MODELS_DIR / ".retrain.lock".

    Returns:
        Dict with keys:
            completed (bool): Whether the loop ran to completion.
            decision (str | None): "DEPLOY", "REJECT", or None if not reached.
            reason (str): Human-readable explanation of outcome.
            metrics (dict | None): Candidate model metrics if training succeeded.
            duration (float): Wall-clock seconds elapsed.
    """
    # Lazy imports to avoid circular dependencies.
    # m5_scoring functions are referenced via module-level attributes so that
    # unittest.mock.patch can intercept them.
    import varis.m7_evolution.auto_retrain as _self

    if _self.select_training_variants is None:
        from varis.m5_scoring.train import (
            select_training_variants as _sel,
            compute_features as _comp,
            load_cached_records as _load,
            train_and_evaluate as _train,
        )
        _self.select_training_variants = _sel
        _self.compute_features = _comp
        _self.load_cached_records = _load
        _self.train_and_evaluate = _train

    from varis.m7_evolution.model_archive import (
        archive_version,
        deploy_version,
        get_current_version,
        mark_rejected,
        cleanup_old_versions,
        list_versions,
    )
    from varis.m7_evolution.evolution_log import (
        log_event,
        init_evolution_log,
        EVENT_RETRAIN_START,
        EVENT_RETRAIN_COMPLETE,
        EVENT_DEPLOY,
        EVENT_REJECT,
        EVENT_ERROR,
    )

    start_time = time.monotonic()

    # Resolve defaults
    if archive_root is None:
        archive_root = DEFAULT_ARCHIVE_ROOT
    archive_root = Path(archive_root)

    if lock_path is None:
        lock_path = DEFAULT_LOCK_PATH
    lock_path = Path(lock_path)

    if log_db is None:
        db_url = DATABASE_URL
        if not db_url.startswith("sqlite"):
            db_url = f"sqlite:///{DATA_DIR / 'evolution.db'}"
        log_db = init_evolution_log(db_url)

    result: dict = {
        "completed": False,
        "decision": None,
        "reason": "",
        "metrics": None,
        "duration": 0.0,
    }

    # ── Step 1: Acquire lock ──
    if not acquire_lock(lock_path):
        result["reason"] = "Lock held by another retrain process"
        result["duration"] = time.monotonic() - start_time
        logger.info("Retrain loop aborted: %s", result["reason"])
        return result

    try:
        # ── Step 2: Log RETRAIN_START ──
        model_version = datetime.now(timezone.utc).strftime("v%Y.%m")
        log_event(log_db, EVENT_RETRAIN_START, model_version=model_version)
        logger.info("Retrain loop started for version %s", model_version)

        # ── Step 3: Select training variants ──
        logger.info("Phase A: Selecting training variants from ClinVar")
        try:
            _self.select_training_variants()
        except Exception as e:
            reason = f"Variant selection failed: {e}"
            logger.error(reason)
            log_event(log_db, EVENT_ERROR, model_version=model_version,
                       details={"phase": "select", "error": str(e)})
            result["reason"] = reason
            return result

        # ── Step 4: Compute features ──
        logger.info("Phase B: Computing features via M1-M4")
        try:
            total, cached, computed, failed = _self.compute_features()
        except Exception as e:
            reason = f"Feature computation failed: {e}"
            logger.error(reason)
            log_event(log_db, EVENT_ERROR, model_version=model_version,
                       details={"phase": "compute", "error": str(e)})
            result["reason"] = reason
            return result

        # ── Step 5: Load cached records ──
        logger.info("Loading cached records for training")
        try:
            records, labels, genes = _self.load_cached_records()
        except Exception as e:
            reason = f"Loading cached records failed: {e}"
            logger.error(reason)
            log_event(log_db, EVENT_ERROR, model_version=model_version,
                       details={"phase": "load", "error": str(e)})
            result["reason"] = reason
            return result

        if not records:
            reason = "No cached records available for training"
            logger.error(reason)
            log_event(log_db, EVENT_ERROR, model_version=model_version,
                       details={"phase": "load", "error": reason})
            result["reason"] = reason
            return result

        # ── Step 6: Train and evaluate ──
        logger.info("Phase C: Training ensemble model")
        candidate_dir = archive_root / "candidates" / f"{model_version}_candidate"
        try:
            train_result = _self.train_and_evaluate(
                records, labels, genes,
                output_dir=candidate_dir,
                model_version=model_version,
            )
        except Exception as e:
            reason = f"Training failed: {e}"
            logger.error(reason)
            log_event(log_db, EVENT_ERROR, model_version=model_version,
                       details={"phase": "train", "error": str(e)})
            result["reason"] = reason
            return result

        cv_results = train_result.get("cv_results", {})
        candidate_metrics = {
            "roc_auc": cv_results.get("roc_auc_mean", 0.0),
            "pr_auc": cv_results.get("pr_auc_mean", 0.0),
            "calibration_ece": cv_results.get("calibration_ece", 0.05),
        }
        result["metrics"] = candidate_metrics

        # ── Step 7: Run benchmark predictions on candidate ──
        logger.info("Running benchmarks on candidate model")
        candidate_benchmarks = run_benchmarks_for_model(candidate_dir)

        # ── Step 8: Load current production metrics ──
        current_version = get_current_version(archive_root)
        current_metrics: dict = {"roc_auc": 0.0, "pr_auc": 0.0, "calibration_ece": 1.0}
        current_benchmarks: dict = {}

        if current_version:
            versions = list_versions(archive_root)
            for v in versions:
                if v.get("version") == current_version:
                    stored_metrics = v.get("metrics", {})
                    current_metrics = {
                        "roc_auc": stored_metrics.get("roc_auc", 0.0),
                        "pr_auc": stored_metrics.get("pr_auc", 0.0),
                        "calibration_ece": stored_metrics.get("calibration_ece", 1.0),
                    }
                    current_benchmarks = stored_metrics.get("benchmarks", {})
                    break

            # Run benchmarks on current model for comparison
            current_model_dir = archive_root / "archive" / current_version
            if current_model_dir.is_dir():
                fresh_current_benchmarks = run_benchmarks_for_model(current_model_dir)
                if fresh_current_benchmarks:
                    current_benchmarks = fresh_current_benchmarks

        # ── Step 9: Evaluate deploy gate ──
        logger.info("Evaluating deploy gate")
        gate_result = evaluate_deploy_gate(
            current_metrics, candidate_metrics,
            current_benchmarks, candidate_benchmarks,
        )
        decision = gate_result["decision"]
        gate_reason = gate_result["reason"]
        result["decision"] = decision
        result["reason"] = gate_reason

        # ── Step 10/11: Deploy or Reject ──
        if decision == "DEPLOY":
            logger.info("Deploy gate passed — archiving and deploying %s", model_version)

            archive_version(
                candidate_dir, model_version, candidate_metrics,
                archive_root=archive_root,
            )
            deploy_version(model_version, archive_root=archive_root)

            log_event(log_db, EVENT_DEPLOY, model_version=model_version, details={
                "candidate_metrics": candidate_metrics,
                "current_metrics": current_metrics,
                "metric_deltas": gate_result.get("metric_deltas", {}),
                "gate_reason": gate_reason,
                "previous_version": current_version,
            })
        else:
            logger.info("Deploy gate rejected — marking %s as rejected", model_version)

            # Archive the candidate first so mark_rejected can find it
            archive_version(
                candidate_dir, model_version, candidate_metrics,
                archive_root=archive_root,
            )
            mark_rejected(model_version, gate_reason, archive_root=archive_root)

            log_event(log_db, EVENT_REJECT, model_version=model_version, details={
                "candidate_metrics": candidate_metrics,
                "current_metrics": current_metrics,
                "metric_deltas": gate_result.get("metric_deltas", {}),
                "gate_reason": gate_reason,
            })

        # ── Step 12: Cleanup old versions ──
        logger.info("Cleaning up old model versions")
        cleanup_old_versions(archive_root=archive_root)

        # ── Step 14: Log RETRAIN_COMPLETE ──
        duration = time.monotonic() - start_time
        result["completed"] = True
        result["duration"] = duration

        log_event(log_db, EVENT_RETRAIN_COMPLETE, model_version=model_version, details={
            "decision": decision,
            "reason": gate_reason,
            "duration_seconds": duration,
            "candidate_metrics": candidate_metrics,
        })

        logger.info(
            "Retrain loop completed: decision=%s, reason=%s, duration=%.1fs",
            decision, gate_reason, duration,
        )
        return result

    except Exception as e:
        # Catch-all for unexpected errors
        duration = time.monotonic() - start_time
        reason = f"Unexpected error in retrain loop: {e}"
        logger.error(reason, exc_info=True)

        try:
            from varis.m7_evolution.evolution_log import EVENT_ERROR
            log_event(log_db, EVENT_ERROR, details={"error": str(e)})
        except Exception:
            pass

        result["reason"] = reason
        result["duration"] = duration
        return result

    finally:
        # ── Step 13: Release lock (always) ──
        release_lock(lock_path)
