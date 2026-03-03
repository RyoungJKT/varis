"""Ensemble — Three-model gradient-boosted tree ensemble.

CatBoost (primary): Best on mixed feature types, native categorical handling.
XGBoost (secondary): Industry standard, robust SHAP integration.
LightGBM (secondary): Fastest training, often highest raw accuracy on tabular data.

Final prediction = calibrated average of all three model predictions.
Model disagreement is flagged — signals a variant that warrants closer review.
"""
import json
import logging
import pickle
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from catboost import CatBoostClassifier
from lightgbm import LGBMClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from xgboost import XGBClassifier

from varis.config import (
    CATBOOST_PARAMS,
    LIGHTGBM_PARAMS,
    MODELS_DIR,
    XGBOOST_PARAMS,
)
from varis.models.variant_record import VariantRecord

logger = logging.getLogger(__name__)


def _classify(score: float) -> str:
    """Convert ensemble score to classification.

    Args:
        score: Calibrated ensemble score in [0, 1].

    Returns:
        'likely_pathogenic' (>0.8), 'uncertain' (0.2-0.8), or 'likely_benign' (<0.2).
    """
    if score > 0.8:
        return "likely_pathogenic"
    elif score < 0.2:
        return "likely_benign"
    return "uncertain"


def _compute_model_agreement(scores: dict[str, float]) -> str:
    """Determine agreement level between models.

    Args:
        scores: Dict of model_name → predicted probability.

    Returns:
        'high' if max-min spread < 0.1, 'moderate' if < 0.2, 'low' otherwise.
    """
    vals = list(scores.values())
    spread = max(vals) - min(vals)
    if spread < 0.1:
        return "high"
    elif spread < 0.2:
        return "moderate"
    return "low"


def _train_catboost(X: pd.DataFrame, y: pd.Series, params: dict) -> CatBoostClassifier:
    """Train CatBoost classifier with given parameters.

    Args:
        X: Training features.
        y: Training labels (0=benign, 1=pathogenic).
        params: CatBoost parameters.

    Returns:
        Trained CatBoostClassifier.
    """
    model = CatBoostClassifier(**params)
    model.fit(X, y, verbose=0)
    return model


def _train_xgboost(X: pd.DataFrame, y: pd.Series, params: dict) -> XGBClassifier:
    """Train XGBoost classifier with given parameters.

    Args:
        X: Training features.
        y: Training labels (0=benign, 1=pathogenic).
        params: XGBoost parameters.

    Returns:
        Trained XGBClassifier.
    """
    model = XGBClassifier(**params)
    model.fit(X, y, verbose=False)
    return model


def _train_lightgbm(X: pd.DataFrame, y: pd.Series, params: dict) -> LGBMClassifier:
    """Train LightGBM classifier with given parameters.

    Args:
        X: Training features.
        y: Training labels (0=benign, 1=pathogenic).
        params: LightGBM parameters.

    Returns:
        Trained LGBMClassifier.
    """
    model = LGBMClassifier(**params)
    model.fit(X, y)
    return model


def train_ensemble(
    X: pd.DataFrame,
    y: pd.Series,
    output_dir: Path = MODELS_DIR,
) -> dict:
    """Train all three models, fit calibrator, save all artifacts.

    Args:
        X: Training features (DataFrame).
        y: Training labels (Series, 0=benign, 1=pathogenic).
        output_dir: Where to save trained model weights.

    Returns:
        Dict of training metrics for each model and the ensemble.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Split: 80% train, 20% calibration
    X_train, X_cal, y_train, y_cal = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y,
    )

    # Train individual models
    logger.info("Training CatBoost...")
    cb_model = _train_catboost(X_train, y_train, CATBOOST_PARAMS)

    logger.info("Training XGBoost...")
    xgb_model = _train_xgboost(X_train, y_train, XGBOOST_PARAMS)

    logger.info("Training LightGBM...")
    lgbm_model = _train_lightgbm(X_train, y_train, LIGHTGBM_PARAMS)

    # Get raw probabilities on calibration set
    cb_probs = cb_model.predict_proba(X_cal)[:, 1]
    xgb_probs = xgb_model.predict_proba(X_cal)[:, 1]
    lgbm_probs = lgbm_model.predict_proba(X_cal)[:, 1]
    raw_avg = (cb_probs + xgb_probs + lgbm_probs) / 3.0

    # Fit Platt calibrator (logistic regression on raw ensemble average)
    calibrator = LogisticRegression(random_state=42)
    calibrator.fit(raw_avg.reshape(-1, 1), y_cal)

    # Save models
    cb_model.save_model(str(output_dir / "catboost_model.cbm"))
    xgb_model.save_model(str(output_dir / "xgboost_model.json"))
    lgbm_model.booster_.save_model(str(output_dir / "lightgbm_model.txt"))

    with open(output_dir / "calibrator.pkl", "wb") as f:
        pickle.dump(calibrator, f)

    with open(output_dir / "feature_columns.json", "w") as f:
        json.dump(list(X.columns), f)

    # Training metadata
    metadata = {
        "n_samples": len(X),
        "n_features": len(X.columns),
        "feature_columns": list(X.columns),
        "class_distribution": {
            "benign": int((y == 0).sum()),
            "pathogenic": int((y == 1).sum()),
        },
    }
    with open(output_dir / "training_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    logger.info(f"Ensemble saved to {output_dir}")
    return metadata


def load_ensemble(model_dir: Path = MODELS_DIR) -> dict:
    """Load trained model weights from disk.

    Args:
        model_dir: Directory containing model files.

    Returns:
        Dict with keys: catboost, xgboost, lightgbm, calibrator, feature_columns.
    """
    model_dir = Path(model_dir)

    cb_model = CatBoostClassifier()
    cb_model.load_model(str(model_dir / "catboost_model.cbm"))

    xgb_model = XGBClassifier()
    xgb_model.load_model(str(model_dir / "xgboost_model.json"))

    import lightgbm as lgb
    lgbm_booster = lgb.Booster(model_file=str(model_dir / "lightgbm_model.txt"))

    with open(model_dir / "calibrator.pkl", "rb") as f:
        calibrator = pickle.load(f)

    with open(model_dir / "feature_columns.json") as f:
        feature_columns = json.load(f)

    return {
        "catboost": cb_model,
        "xgboost": xgb_model,
        "lightgbm": lgbm_booster,
        "calibrator": calibrator,
        "feature_columns": feature_columns,
    }


def predict_from_models(models: dict, features_dict: dict) -> dict:
    """Run all three models and compute ensemble prediction.

    Args:
        models: Dict from load_ensemble().
        features_dict: Feature dict (column_name → value).

    Returns:
        Dict with score_ensemble, score_catboost, score_xgboost,
        score_lightgbm, classification, model_agreement,
        confidence_lower, confidence_upper.
    """
    columns = models["feature_columns"]
    values = [features_dict.get(col) for col in columns]
    X = pd.DataFrame([values], columns=columns)

    # Individual model predictions
    cb_prob = float(models["catboost"].predict_proba(X)[0, 1])
    xgb_prob = float(models["xgboost"].predict_proba(X)[0, 1])

    # LightGBM booster uses raw predict (returns probability directly)
    lgbm_raw = models["lightgbm"].predict(X)
    lgbm_prob = float(lgbm_raw[0])

    # Raw ensemble average
    raw_avg = (cb_prob + xgb_prob + lgbm_prob) / 3.0

    # Calibrate
    calibrated = float(
        models["calibrator"].predict_proba(np.array([[raw_avg]]))[0, 1]
    )
    calibrated = max(0.0, min(1.0, calibrated))

    # Individual scores dict for agreement
    individual = {
        "catboost": cb_prob,
        "xgboost": xgb_prob,
        "lightgbm": lgbm_prob,
    }

    # Confidence interval from model spread
    vals = list(individual.values())
    spread = max(vals) - min(vals)
    confidence_lower = max(0.0, calibrated - spread / 2)
    confidence_upper = min(1.0, calibrated + spread / 2)

    return {
        "score_ensemble": calibrated,
        "score_catboost": cb_prob,
        "score_xgboost": xgb_prob,
        "score_lightgbm": lgbm_prob,
        "classification": _classify(calibrated),
        "model_agreement": _compute_model_agreement(individual),
        "confidence_lower": confidence_lower,
        "confidence_upper": confidence_upper,
    }


def predict(
    variant_record: VariantRecord,
    features: dict,
    model_dir: Optional[Path] = None,
) -> VariantRecord:
    """Convenience wrapper: load models, predict, populate VariantRecord.

    Args:
        variant_record: The current investigation record.
        features: Feature dict from feature_extractor.
        model_dir: Directory with model files. Defaults to MODELS_DIR.

    Returns:
        VariantRecord with ensemble scores populated.
    """
    if model_dir is None:
        model_dir = MODELS_DIR

    models = load_ensemble(model_dir)
    scores = predict_from_models(models, features)

    variant_record.score_ensemble = scores["score_ensemble"]
    variant_record.score_catboost = scores["score_catboost"]
    variant_record.score_xgboost = scores["score_xgboost"]
    variant_record.score_lightgbm = scores["score_lightgbm"]
    variant_record.classification = scores["classification"]
    variant_record.model_agreement = scores["model_agreement"]
    variant_record.confidence_lower = scores["confidence_lower"]
    variant_record.confidence_upper = scores["confidence_upper"]
    variant_record.features_used = len(
        [v for v in features.values() if v is not None]
    )

    return variant_record
