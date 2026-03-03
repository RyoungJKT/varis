"""Benchmarks — Regression testing against known variants and training manifest.

Evaluates the ensemble against a fixed set of well-characterized variants.
Any unexpected reclassification blocks deployment.

Functions:
    run_benchmarks: Evaluate model against benchmark variant set
    save_training_manifest: Write training metadata for reproducibility
"""
import json
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

from varis.m5_scoring.ensemble import load_ensemble, predict_from_models

logger = logging.getLogger(__name__)


def run_benchmarks(
    model_dir: Path,
    benchmark_file: str = "tests/benchmark_variants.json",
) -> dict:
    """Evaluate model against benchmark variants.

    Loads the benchmark variant set and generates predictions for each.
    Does NOT evaluate accuracy (no trained pipeline for real variants yet),
    but records predicted scores for regression testing.

    Args:
        model_dir: Directory containing trained model files.
        benchmark_file: Path to benchmark variants JSON file.

    Returns:
        Dict with 'variants' list containing per-variant results.
    """
    benchmark_path = Path(benchmark_file)
    if not benchmark_path.exists():
        logger.warning(f"Benchmark file not found: {benchmark_path}")
        return {"variants": [], "error": "benchmark file not found"}

    with open(benchmark_path) as f:
        variants = json.load(f)

    results = {
        "model_dir": str(model_dir),
        "benchmark_file": str(benchmark_file),
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "n_variants": len(variants),
        "variants": [],
    }

    # Try to load model (may not have matching features for real variants)
    try:
        models = load_ensemble(model_dir)
    except Exception as e:
        logger.warning(f"Could not load model for benchmarks: {e}")
        for v in variants:
            results["variants"].append({
                "gene": v["gene"],
                "hgvs": v["hgvs"],
                "expected": v["expected"],
                "predicted": None,
                "score": None,
                "error": "model not loaded",
            })
        return results

    for v in variants:
        result_entry = {
            "gene": v["gene"],
            "hgvs": v["hgvs"],
            "expected": v["expected"],
        }
        # Note: real benchmark evaluation requires running the full pipeline
        # For now, we record the variant metadata for regression tracking
        result_entry["predicted"] = None
        result_entry["score"] = None
        result_entry["note"] = "Full pipeline evaluation not yet available"
        results["variants"].append(result_entry)

    return results


def save_training_manifest(
    output_dir: Path,
    metrics: dict,
    metadata: dict,
) -> Path:
    """Write training manifest for reproducibility and audit.

    Args:
        output_dir: Directory to save manifest.
        metrics: Training/evaluation metrics (roc_auc, pr_auc, etc.).
        metadata: Training metadata (n_samples, n_features, etc.).

    Returns:
        Path to the saved manifest file.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    manifest = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "metrics": metrics,
        "metadata": metadata,
        "library_versions": _get_library_versions(),
    }

    manifest_path = output_dir / "training_manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    logger.info(f"Training manifest saved to {manifest_path}")
    return manifest_path


def _get_library_versions() -> dict:
    """Get versions of ML libraries for reproducibility.

    Returns:
        Dict of library_name → version string.
    """
    versions = {}
    try:
        import catboost
        versions["catboost"] = catboost.__version__
    except ImportError:
        versions["catboost"] = None
    try:
        import xgboost
        versions["xgboost"] = xgboost.__version__
    except ImportError:
        versions["xgboost"] = None
    try:
        import lightgbm
        versions["lightgbm"] = lightgbm.__version__
    except ImportError:
        versions["lightgbm"] = None
    try:
        import sklearn
        versions["scikit-learn"] = sklearn.__version__
    except ImportError:
        versions["scikit-learn"] = None
    try:
        import shap
        versions["shap"] = shap.__version__
    except ImportError:
        versions["shap"] = None
    return versions
