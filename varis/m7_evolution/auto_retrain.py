"""Auto-Retrain (Loop 1) — Fully autonomous monthly model retraining with governance.

Scans ClinVar for newly expert-reviewed variants, runs them through the pipeline,
retrains all three models, compares metrics, runs regression tests, and deploys
only if all gates pass.

Release governance:
  - Every model is version-tagged (e.g. v2026.03)
  - Deployment requires ALL metrics to improve AND regression tests to pass
  - Previous versions are preserved in the version archive for rollback
  - Every decision is recorded in the Evolution Log with metric deltas

Priority: 2 (build second — most impactful self-evolution loop)
Frequency: Monthly (triggered by ClinVar update detection)
Autonomy: Fully autonomous — no human intervention required
"""
import json
import logging
import os
import signal
from datetime import datetime, timezone
from pathlib import Path

from varis.m7_evolution.evolution_log import log_event, EVENT_AUTO_RETRAIN

logger = logging.getLogger(__name__)

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


def check_for_new_clinvar_data() -> list[dict] | None:
    """Scan ClinVar for newly expert-reviewed variants since last retrain.

    Returns list of new variant records, or None if no new data.
    """
    pass


def retrain_and_evaluate(new_variants: list[dict]) -> dict:
    """Retrain ensemble on expanded dataset and compare with current model.

    Process:
    1. Run new variants through M1-M4 pipeline for feature extraction
    2. Add to training dataset
    3. Retrain CatBoost + XGBoost + LightGBM with Optuna tuning
    4. Evaluate on held-out test set
    5. Compare all metrics with current model

    Returns:
        Dict with old_metrics, new_metrics, and deploy_decision.
    """
    pass


def run_regression_tests(candidate_model, benchmark_path: str = "tests/benchmark_variants.json") -> dict:
    """Run candidate model against fixed benchmark variant set.

    The benchmark set contains well-characterised variants (pathogenic, benign,
    and known edge cases) that must NEVER be included in training data.

    Any unexpected reclassification of a benchmark variant triggers rejection.

    Returns:
        Dict with passed (bool), total_variants, reclassified (list),
        and details for logging.
    """
    pass


def deploy_if_better(old_metrics: dict, new_metrics: dict, regression_results: dict) -> bool:
    """Deploy new model only if ALL key metrics improve AND regression tests pass.

    Metrics checked: ROC-AUC, PR-AUC, precision, recall.
    If any metric decreases OR any regression test fails, keep current model.

    On deploy:
      - Tag new model with version stamp (v{YYYY}.{MM})
      - Archive current production model to version history
      - Log full metric comparison to Evolution Log

    On reject:
      - Archive candidate model as rejected
      - Log rejection reason with metric comparison

    Returns:
        True if deployed, False if kept current model.
    """
    pass


def rollback_to_previous(reason: str) -> bool:
    """Restore the previous production model from the version archive.

    Called when a deployed model produces unexpected results in production.
    Records the rollback event and reason in the Evolution Log.

    Returns:
        True if rollback succeeded, False if no previous version available.
    """
    pass


def get_model_version_history() -> list[dict]:
    """Return list of all archived model versions with metadata.

    Each entry contains: version tag, deploy date, metrics at deploy time,
    training variant count, and status (production/archived/rejected/rolled-back).
    """
    pass


def run_retrain_loop():
    """Execute the full auto-retrain loop. Called by scheduler.

    Full sequence:
    1. Check for new ClinVar data
    2. If new data: retrain and evaluate candidate model
    3. Run regression tests on candidate
    4. Deploy only if metrics improve AND regression passes
    5. Log everything to Evolution Log with metric deltas
    """
    pass
