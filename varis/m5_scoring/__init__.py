"""M5: ML Scoring Engine — Trains and runs the pathogenicity ensemble.

Three gradient-boosted tree models (CatBoost, XGBoost, LightGBM) trained on
structural features extracted by the pipeline. Ensemble average is the final score.
All three models handle missing features natively — the ensemble trains on
whatever features are available.

Why tree-based instead of deep learning: Varis's argument against AlphaMissense
is that black-box predictions aren't clinically trustworthy. Using deep learning
would undermine that narrative.

Depends on: Any subset of M3 + M4 features (works with 3-15 features)
Populates: All scoring.* and evidence.* fields
"""
import logging

from varis.config import MODELS_DIR

logger = logging.getLogger(__name__)


def run(variant_record):
    """Execute M5: predict pathogenicity using the ML ensemble.

    Pipeline: extract features → predict with ensemble → SHAP explain →
    map evidence tags. Each step is independent — if one fails, the
    others still run.

    Args:
        variant_record: VariantRecord with M3/M4 features populated.

    Returns:
        VariantRecord with scoring and evidence fields populated.
    """
    from varis.m5_scoring.ensemble import load_ensemble, predict_from_models
    from varis.m5_scoring.evidence_mapper import map_evidence_tags
    from varis.m5_scoring.feature_extractor import extract_features
    from varis.m5_scoring.shap_explainer import explain_from_models

    # Extract features
    try:
        features = extract_features(variant_record)
    except Exception as e:
        logger.warning(f"M5 feature extraction failed: {e}")
        variant_record.mark_module_failed("M5")
        return variant_record

    # Load models and predict
    try:
        models = load_ensemble(MODELS_DIR)
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
    except Exception as e:
        logger.warning(f"M5 ensemble prediction failed: {e}")
        variant_record.mark_module_failed("M5")
        return variant_record

    # SHAP explanation
    try:
        top_features = explain_from_models(models, features)
        variant_record.shap_top_features = top_features
    except Exception as e:
        logger.warning(f"SHAP explanation failed: {e}")
        variant_record.mark_module_failed("M5.shap")

    # Evidence tags
    try:
        variant_record = map_evidence_tags(variant_record)
    except Exception as e:
        logger.warning(f"Evidence mapping failed: {e}")
        variant_record.mark_module_failed("M5.evidence")

    variant_record.mark_module_completed("M5")
    return variant_record
