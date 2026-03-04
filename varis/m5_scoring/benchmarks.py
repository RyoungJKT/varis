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

# Cache directory for benchmark variant records (separate from training cache)
BENCHMARK_CACHE_DIR = Path("data/benchmarks")


def run_benchmarks(
    model_dir: Path,
    benchmark_file: str = "tests/benchmark_variants.json",
) -> dict:
    """Evaluate model against benchmark variants.

    For each benchmark variant:
    1. Check benchmark cache for existing VariantRecord
    2. If not cached, run full M1-M4 pipeline and cache the result
    3. Extract features and score with the model
    4. Compare predicted classification to expected

    Args:
        model_dir: Directory containing trained model files.
        benchmark_file: Path to benchmark variants JSON file.

    Returns:
        Dict with 'variants' list containing per-variant results,
        'n_correct' count, and 'accuracy' float.
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
        "n_correct": 0,
        "accuracy": 0.0,
    }

    # Load model
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
                "correct": None,
                "error": "model not loaded",
            })
        return results

    for v in variants:
        result_entry = _evaluate_benchmark_variant(
            v["gene"], v["hgvs"], v["expected"], models
        )
        results["variants"].append(result_entry)
        if result_entry.get("correct"):
            results["n_correct"] += 1

    n_scored = sum(1 for v in results["variants"] if v.get("score") is not None)
    if n_scored > 0:
        results["accuracy"] = results["n_correct"] / n_scored

    return results


def _evaluate_benchmark_variant(
    gene: str,
    hgvs: str,
    expected: str,
    models: dict,
) -> dict:
    """Evaluate a single benchmark variant against the model.

    Args:
        gene: Gene symbol.
        hgvs: HGVS protein notation.
        expected: Expected classification.
        models: Loaded ensemble models dict.

    Returns:
        Dict with gene, hgvs, expected, predicted, score, correct, error.
    """
    from varis.m5_scoring.feature_extractor import extract_features

    variant_id = f"{gene}_{hgvs}"
    entry = {
        "gene": gene,
        "hgvs": hgvs,
        "expected": expected,
        "predicted": None,
        "score": None,
        "correct": None,
        "error": None,
    }

    # Step 1: Get VariantRecord (from cache or pipeline)
    try:
        record = _get_benchmark_record(gene, hgvs)
    except Exception as e:
        entry["error"] = f"pipeline failed: {e}"
        logger.warning(f"Benchmark {variant_id}: pipeline failed: {e}")
        return entry

    # Step 2: Extract features
    try:
        features = extract_features(record)
    except Exception as e:
        entry["error"] = f"feature extraction failed: {e}"
        logger.warning(f"Benchmark {variant_id}: feature extraction failed: {e}")
        return entry

    # Step 3: Score with model
    try:
        scores = predict_from_models(models, features)
        entry["score"] = scores["score_ensemble"]
        entry["predicted"] = scores["classification"]
        entry["correct"] = entry["predicted"] == expected
        logger.info(
            f"Benchmark {variant_id}: score={scores['score_ensemble']:.3f} "
            f"predicted={scores['classification']} expected={expected} "
            f"correct={entry['correct']}"
        )
    except Exception as e:
        entry["error"] = f"prediction failed: {e}"
        logger.warning(f"Benchmark {variant_id}: prediction failed: {e}")

    return entry


def _get_benchmark_record(gene: str, hgvs: str):
    """Get VariantRecord for a benchmark variant, using cache or pipeline.

    Args:
        gene: Gene symbol.
        hgvs: HGVS protein notation.

    Returns:
        VariantRecord with M1-M4 features populated.
    """
    from varis.models.variant_record import VariantRecord

    variant_id = f"{gene}_{hgvs}"
    cache_path = BENCHMARK_CACHE_DIR / f"{variant_id}.json"

    # Check cache
    if cache_path.exists():
        logger.info(f"Benchmark cache hit: {variant_id}")
        with open(cache_path) as f:
            record_dict = json.load(f)
        return VariantRecord.from_dict(record_dict)

    # Run pipeline
    logger.info(f"Benchmark cache miss: {variant_id} — running pipeline")
    from varis.pipeline import investigate

    record = investigate(gene, hgvs)

    # Cache for future runs
    BENCHMARK_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "w") as f:
        json.dump(record.to_dict(), f, default=str)
    logger.info(f"Benchmark cached: {variant_id}")

    return record


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
