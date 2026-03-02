"""SHAP Explainer — Per-prediction feature importance for clinical interpretability.

This is what makes Varis different from AlphaMissense. Every prediction comes
with a ranked list of features showing exactly what drove the decision.
A genetic counselor can say: "The variant is predicted pathogenic primarily
because the ΔΔG is +3.7 kcal/mol and the position is 100% conserved."
"""
import logging
from varis.models.variant_record import VariantRecord
logger = logging.getLogger(__name__)

def explain_prediction(variant_record: VariantRecord, features: dict) -> VariantRecord:
    """Generate SHAP explanation for this prediction.

    Uses the primary CatBoost model for SHAP values. Populates
    shap_top_features with ranked feature importance for this variant.

    Args:
        variant_record: Must have ensemble scores from predict().
        features: Same feature dict used for prediction.

    Returns:
        VariantRecord with shap_top_features populated.
    """
    pass

def generate_shap_summary(feature_matrix, model) -> dict:
    """Generate global SHAP summary across all training data.

    Used for the benchmark report and VarisDB global feature importance display.
    """
    pass
