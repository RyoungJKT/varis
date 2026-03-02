"""Ensemble — Three-model gradient-boosted tree ensemble.

CatBoost (primary): Best on mixed feature types, native categorical handling.
XGBoost (secondary): Industry standard, robust SHAP integration.
LightGBM (secondary): Fastest training, often highest raw accuracy on tabular data.

Final prediction = average of all three model predictions.
Model disagreement is flagged — signals a variant that warrants closer review.
"""
import logging
from pathlib import Path
from varis.models.variant_record import VariantRecord
from varis.config import ENSEMBLE_MODELS, MODELS_DIR
logger = logging.getLogger(__name__)

def predict(variant_record: VariantRecord, features: dict) -> VariantRecord:
    """Run all three models and compute ensemble prediction.

    Args:
        variant_record: The current investigation record.
        features: Feature dict from feature_extractor.

    Returns:
        VariantRecord with ensemble score, individual scores, CI, classification.
    """
    pass

def train_ensemble(feature_matrix, labels, output_dir: Path = MODELS_DIR) -> dict:
    """Train all three models with Optuna hyperparameter tuning.

    Args:
        feature_matrix: Training features (DataFrame).
        labels: Training labels (Series, 0=benign, 1=pathogenic).
        output_dir: Where to save trained model weights.

    Returns:
        Dict of training metrics for each model and the ensemble.
    """
    pass

def _train_catboost(X, y, params: dict) -> object:
    """Train CatBoost classifier with given parameters."""
    pass

def _train_xgboost(X, y, params: dict) -> object:
    """Train XGBoost classifier with given parameters."""
    pass

def _train_lightgbm(X, y, params: dict) -> object:
    """Train LightGBM classifier with given parameters."""
    pass

def _compute_model_agreement(scores: dict[str, float]) -> str:
    """Determine agreement level between models.

    Returns 'high' if max-min spread < 0.1, 'moderate' if < 0.2, 'low' otherwise.
    """
    pass

def _classify(score: float) -> str:
    """Convert ensemble score to classification.

    Returns 'likely_pathogenic' (>0.8), 'uncertain' (0.2-0.8), or 'likely_benign' (<0.2).
    """
    pass

def load_ensemble(model_dir: Path = MODELS_DIR) -> dict:
    """Load trained model weights from disk."""
    pass
