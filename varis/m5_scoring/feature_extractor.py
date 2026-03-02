"""Feature Extractor — Builds the feature vector from VariantRecord for ML input.

Extracts the ~15 structural features from the VariantRecord and formats them
for the gradient-boosted tree ensemble. Missing features remain None —
CatBoost, XGBoost, and LightGBM all handle this natively.
"""
import logging
from varis.models.variant_record import VariantRecord
logger = logging.getLogger(__name__)

def extract_features(variant_record: VariantRecord) -> dict:
    """Extract ML feature dict from VariantRecord.

    Returns:
        Dict of feature_name -> value (or None for missing features).
        The ensemble handles None natively.
    """
    return variant_record.get_structural_features()

def build_feature_matrix(records: list[VariantRecord]) -> tuple:
    """Build feature matrix and label vector from multiple VariantRecords.

    Used during training. Extracts features from each record and
    assembles into a DataFrame suitable for model training.

    Args:
        records: List of VariantRecords with known classifications.

    Returns:
        Tuple of (feature_dataframe, label_series).
    """
    pass
