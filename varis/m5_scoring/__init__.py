"""M5: ML Scoring Engine — Trains and runs the pathogenicity ensemble.

Three gradient-boosted tree models (CatBoost, XGBoost, LightGBM) trained on
structural features extracted by the pipeline. Ensemble average is the final score.
All three models handle missing features natively — the ensemble trains on
whatever features are available.

Why tree-based instead of deep learning: Varis's argument against AlphaMissense
is that black-box predictions aren't clinically trustworthy. Using deep learning
would undermine that narrative.

Depends on: Any subset of M3 + M4 features (works with 3-15 features)
Populates: All scoring.* and acmg.* fields
"""
import logging
logger = logging.getLogger(__name__)

def run(variant_record):
    """Execute M5: predict pathogenicity using the ML ensemble."""
    from varis.m5_scoring.feature_extractor import extract_features
    from varis.m5_scoring.ensemble import predict
    from varis.m5_scoring.shap_explainer import explain_prediction
    from varis.m5_scoring.acmg_mapper import map_acmg_codes

    try:
        features = extract_features(variant_record)
        variant_record = predict(variant_record, features)
    except Exception as e:
        logger.warning(f"M5 ensemble prediction failed: {e}")
        variant_record.mark_module_failed("M5.ensemble")
        return variant_record

    try:
        variant_record = explain_prediction(variant_record, features)
    except Exception as e:
        logger.warning(f"SHAP explanation failed: {e}")
        variant_record.mark_module_failed("M5.shap")

    try:
        variant_record = map_acmg_codes(variant_record)
    except Exception as e:
        logger.warning(f"ACMG mapping failed: {e}")
        variant_record.mark_module_failed("M5.acmg")

    variant_record.mark_module_completed("M5")
    return variant_record
