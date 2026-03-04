"""Feature Extractor — Builds the feature vector from VariantRecord for ML input.

Extracts structural features from the VariantRecord, applies one-hot encoding
for categoricals, converts domain_name to binary in_domain flag, and provides
strict column alignment for model input.

Functions:
    extract_features: VariantRecord → ordered dict of features
    build_feature_vector: dict + columns → ordered list for model input
    simulate_missingness: VariantRecord + groups → VariantRecord with dropped features
"""
import copy
import logging
import random
from collections import OrderedDict
from typing import Optional

from varis.models.variant_record import NullReason, VariantRecord

logger = logging.getLogger(__name__)

# One-hot categories and their possible values
ONE_HOT_CATEGORIES = {
    "secondary_structure": ["helix", "sheet", "coil"],
    "burial_category": ["core", "surface"],
    "charge_change": ["+ve → neutral", "neutral → +ve", "+ve → -ve",
                      "-ve → +ve", "-ve → neutral", "neutral → -ve",
                      "no_change"],
    "mutation_site_confidence_bucket": ["very_high", "confident", "low", "very_low"],
}

# Feature groups for missingness simulation
FEATURE_GROUPS = {
    "ddg": ["ddg_evoef2", "ddg_foldx", "ddg_pyrosetta"],
    "sasa": ["solvent_accessibility_relative", "burial_category"],
    "dssp": ["secondary_structure_name"],
    "conservation": ["conservation_score"],
    "contacts": ["contacts_wt", "hbonds_wt", "packing_density"],
    "domain": ["domain_name"],
}


def extract_features(variant_record: VariantRecord) -> dict:
    """Extract ML feature dict from VariantRecord with one-hot encoding.

    Gets raw features via record.get_ml_features(), then:
    - One-hot encodes categorical columns
    - Converts domain_name to binary in_domain flag
    - Removes raw categorical columns
    - Returns deterministically ordered dict

    Args:
        variant_record: The VariantRecord to extract features from.

    Returns:
        Ordered dict of feature_name → value (or None for missing).
    """
    raw = variant_record.get_ml_features()

    features = OrderedDict()

    # Numeric features (pass through directly)
    numeric_keys = [
        "ddg_evoef2", "ddg_foldx", "ddg_pyrosetta",
        "solvent_accessibility_relative",
        "contacts_wt", "hbonds_wt", "packing_density",
        "conservation_score",
        "mutation_site_plddt",
        "gnomad_frequency", "alphamissense_score",
    ]
    for key in numeric_keys:
        features[key] = raw.get(key)

    # Binary: in_domain from domain_name
    domain_name = raw.get("domain_name")
    if domain_name is None:
        features["in_domain"] = None
    else:
        features["in_domain"] = 1 if domain_name else 0

    # Availability flags
    availability_keys = [
        "ddg_available", "sasa_available", "dssp_available",
        "conservation_available", "domain_available", "contacts_available",
    ]
    for key in availability_keys:
        features[key] = raw.get(key)

    # One-hot encode categoricals
    for cat_name, categories in ONE_HOT_CATEGORIES.items():
        raw_value = raw.get(cat_name)
        for cat_value in categories:
            col_name = f"{cat_name}_{cat_value}"
            if raw_value is None:
                features[col_name] = None
            elif raw_value == cat_value:
                features[col_name] = 1
            else:
                features[col_name] = 0

    return features


def build_feature_vector(features: dict, columns: list[str]) -> list:
    """Convert feature dict to ordered list matching column order.

    Args:
        features: Feature dict from extract_features().
        columns: Ordered list of column names (from training).

    Returns:
        List of values in column order. None values are preserved.

    Raises:
        ValueError: If a required column is missing from features.
    """
    vector = []
    for col in columns:
        if col not in features:
            raise ValueError(f"Missing feature column: '{col}'")
        vector.append(features[col])
    return vector


def simulate_missingness(
    variant_record: VariantRecord,
    groups: Optional[list[str]] = None,
    rate: float = 0.15,
) -> VariantRecord:
    """Simulate missing feature blocks for training robustness.

    Randomly drops entire feature groups and sets their availability
    flags to False. This teaches the model to handle real-world
    missingness without learning artifacts.

    Args:
        variant_record: The record to modify (will be deep-copied).
        groups: Specific groups to drop. If None, randomly selects
                based on rate.
        rate: Probability of dropping each group (when groups is None).

    Returns:
        Modified VariantRecord with dropped features.
    """
    record = copy.deepcopy(variant_record)

    if groups is None:
        groups = [g for g in FEATURE_GROUPS if random.random() < rate]

    for group in groups:
        if group not in FEATURE_GROUPS:
            logger.warning(f"Unknown feature group: {group}")
            continue

        # Set all features in this group to None
        for field_name in FEATURE_GROUPS[group]:
            record.set_with_reason(field_name, None, NullReason.INTENTIONALLY_SKIPPED)

        # Set availability flag to False
        record.set_feature_status(group, False, NullReason.INTENTIONALLY_SKIPPED)

    return record
