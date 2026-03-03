"""SHAP Explainer — Per-prediction feature importance for clinical interpretability.

This is what makes Varis different from AlphaMissense. Every prediction comes
with a ranked list of features showing exactly what drove the decision.
A genetic counselor can say: "The variant is predicted pathogenic primarily
because the DDG is +3.7 kcal/mol and the position is 100% conserved."

Uses TreeSHAP for each base model and averages SHAP values across models.
"""
import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import shap

from varis.config import MODELS_DIR
from varis.models.variant_record import VariantRecord

logger = logging.getLogger(__name__)


def explain_from_models(
    models: dict,
    features_dict: dict,
    top_n: int = 10,
) -> list[dict]:
    """Generate SHAP explanation for a single prediction.

    Computes TreeSHAP for each base model, averages SHAP values,
    and returns the top N features sorted by |SHAP|.

    Args:
        models: Dict from ensemble.load_ensemble().
        features_dict: Feature dict (column_name → value).
        top_n: Number of top features to return.

    Returns:
        List of dicts with keys: feature, value, shap.
        Sorted by |shap| descending.
    """
    columns = models["feature_columns"]
    values = [features_dict.get(col) for col in columns]
    X = pd.DataFrame([values], columns=columns)

    shap_arrays = []

    # CatBoost TreeSHAP
    try:
        cb_explainer = shap.TreeExplainer(models["catboost"])
        cb_shap = cb_explainer.shap_values(X)
        if isinstance(cb_shap, list):
            cb_shap = cb_shap[1]  # Class 1 (pathogenic)
        shap_arrays.append(np.array(cb_shap).flatten())
    except Exception as e:
        logger.warning(f"CatBoost SHAP failed: {e}")

    # XGBoost TreeSHAP
    try:
        xgb_explainer = shap.TreeExplainer(models["xgboost"])
        xgb_shap = xgb_explainer.shap_values(X)
        if isinstance(xgb_shap, list):
            xgb_shap = xgb_shap[1]
        shap_arrays.append(np.array(xgb_shap).flatten())
    except Exception as e:
        logger.warning(f"XGBoost SHAP failed: {e}")

    # LightGBM TreeSHAP
    try:
        lgbm_explainer = shap.TreeExplainer(models["lightgbm"])
        lgbm_shap = lgbm_explainer.shap_values(X)
        if isinstance(lgbm_shap, list):
            lgbm_shap = lgbm_shap[1]
        shap_arrays.append(np.array(lgbm_shap).flatten())
    except Exception as e:
        logger.warning(f"LightGBM SHAP failed: {e}")

    if not shap_arrays:
        logger.warning("All SHAP computations failed")
        return []

    # Average SHAP values across models
    avg_shap = np.mean(shap_arrays, axis=0)

    # Build result list
    results = []
    for i, col in enumerate(columns):
        results.append({
            "feature": col,
            "value": features_dict.get(col),
            "shap": float(avg_shap[i]),
        })

    # Sort by |shap| descending, take top N
    results.sort(key=lambda x: abs(x["shap"]), reverse=True)
    return results[:top_n]


def explain_prediction(
    variant_record: VariantRecord,
    features: dict,
    model_dir: Optional[Path] = None,
) -> VariantRecord:
    """Generate SHAP explanation and populate VariantRecord.

    Args:
        variant_record: Must have ensemble scores from predict().
        features: Same feature dict used for prediction.
        model_dir: Directory with model files. Defaults to MODELS_DIR.

    Returns:
        VariantRecord with shap_top_features populated.
    """
    if model_dir is None:
        model_dir = MODELS_DIR

    from varis.m5_scoring.ensemble import load_ensemble
    models = load_ensemble(model_dir)
    top_features = explain_from_models(models, features)
    variant_record.shap_top_features = top_features
    return variant_record


def compute_global_importance(
    models: dict,
    X: pd.DataFrame,
) -> dict:
    """Compute global feature importance via mean |SHAP| across samples.

    Args:
        models: Dict from ensemble.load_ensemble().
        X: Feature matrix (DataFrame) to compute importance over.

    Returns:
        Dict of feature_name → mean absolute SHAP value.
    """
    columns = models["feature_columns"]

    all_shap = []

    # Compute SHAP for each model across all samples
    for model_key in ["catboost", "xgboost", "lightgbm"]:
        try:
            explainer = shap.TreeExplainer(models[model_key])
            shap_values = explainer.shap_values(X[columns])
            if isinstance(shap_values, list):
                shap_values = shap_values[1]
            all_shap.append(np.abs(np.array(shap_values)))
        except Exception as e:
            logger.warning(f"{model_key} global SHAP failed: {e}")

    if not all_shap:
        return {col: 0.0 for col in columns}

    # Average across models, then mean across samples
    avg_abs_shap = np.mean(all_shap, axis=0)  # (n_samples, n_features)
    mean_importance = np.mean(avg_abs_shap, axis=0)  # (n_features,)

    return {col: float(mean_importance[i]) for i, col in enumerate(columns)}
