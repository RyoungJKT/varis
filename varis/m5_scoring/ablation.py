"""Ablation — Feature group removal to prove each module adds value.

Drops entire feature groups (e.g., all conservation features), retrains
the ensemble, and compares metrics. If removing a feature group doesn't
hurt performance, that module may not be contributing.

Functions:
    drop_feature_group: Remove a feature group from DataFrame
    run_ablation: Full ablation study across specified groups
"""
import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold

from varis.m5_scoring.ensemble import train_ensemble, load_ensemble, predict_from_models

logger = logging.getLogger(__name__)

# Feature group → (feature columns, availability flag column)
ABLATION_GROUPS = {
    "conservation": {
        "features": ["conservation_score"],
        "flag": "conservation_available",
    },
    "ddg": {
        "features": ["ddg_evoef2", "ddg_foldx", "ddg_pyrosetta"],
        "flag": "ddg_available",
    },
    "sasa": {
        "features": ["solvent_accessibility_relative", "burial_category_core", "burial_category_surface"],
        "flag": "sasa_available",
    },
    "dssp": {
        "features": ["secondary_structure_helix", "secondary_structure_sheet", "secondary_structure_coil"],
        "flag": "dssp_available",
    },
    "contacts": {
        "features": ["contacts_wt", "hbonds_wt", "packing_density"],
        "flag": "contacts_available",
    },
    "domain": {
        "features": ["in_domain"],
        "flag": "domain_available",
    },
}


def drop_feature_group(X: pd.DataFrame, group_name: str) -> pd.DataFrame:
    """Remove a feature group by setting columns to NaN/False.

    Args:
        X: Feature DataFrame.
        group_name: Name of group to drop (key in ABLATION_GROUPS).

    Returns:
        Modified copy of X with group features set to NaN
        and availability flag set to False.
    """
    result = X.copy()
    group = ABLATION_GROUPS.get(group_name, {})

    for col in group.get("features", []):
        if col in result.columns:
            result[col] = np.nan

    flag_col = group.get("flag")
    if flag_col and flag_col in result.columns:
        result[flag_col] = False

    return result


def run_ablation(
    X: pd.DataFrame,
    y: pd.Series,
    genes: pd.Series,
    groups: Optional[list[str]] = None,
    output_dir: Optional[Path] = None,
) -> dict:
    """Run ablation study: remove each feature group, retrain, evaluate.

    Args:
        X: Feature DataFrame.
        y: Labels (0=benign, 1=pathogenic).
        genes: Gene labels for stratified splitting.
        groups: Groups to ablate. Defaults to all in ABLATION_GROUPS.
        output_dir: Directory for saving ablation models.

    Returns:
        Dict of {f"without_{group}": {"roc_auc": float, ...}}.
    """
    if groups is None:
        groups = list(ABLATION_GROUPS.keys())
    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    results = {}

    for group_name in groups:
        logger.info(f"Ablation: removing {group_name}")
        X_ablated = drop_feature_group(X, group_name)

        # Train on ablated data and evaluate via cross-validation
        try:
            aucs = []
            sgkf = StratifiedGroupKFold(n_splits=3, shuffle=True, random_state=42)

            for fold_idx, (train_idx, test_idx) in enumerate(sgkf.split(X_ablated, y, genes)):
                X_train = X_ablated.iloc[train_idx]
                y_train = y.iloc[train_idx]
                X_test = X_ablated.iloc[test_idx]
                y_test = y.iloc[test_idx]

                fold_dir = output_dir / f"ablation_{group_name}_fold{fold_idx}" if output_dir else Path(f"/tmp/ablation_{group_name}_fold{fold_idx}")
                fold_dir.mkdir(parents=True, exist_ok=True)

                train_ensemble(X_train, y_train, output_dir=fold_dir)
                models = load_ensemble(fold_dir)

                preds = []
                for _, row in X_test.iterrows():
                    scores = predict_from_models(models, row.to_dict())
                    preds.append(scores["score_ensemble"])

                auc = roc_auc_score(y_test, preds)
                aucs.append(auc)

            results[f"without_{group_name}"] = {
                "roc_auc": float(np.mean(aucs)),
                "roc_auc_std": float(np.std(aucs)),
                "folds": len(aucs),
            }
            logger.info(f"  {group_name} ablated: AUC={np.mean(aucs):.3f} ± {np.std(aucs):.3f}")

        except Exception as e:
            logger.warning(f"Ablation failed for {group_name}: {e}")
            results[f"without_{group_name}"] = {
                "roc_auc": None,
                "error": str(e),
            }

    return results
