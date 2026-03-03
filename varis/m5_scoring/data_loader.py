"""Data Loader — ClinVar variant loading and gene-stratified splitting.

Provides gene-stratified cross-validation (StratifiedGroupKFold) to ensure
no gene appears in both train and test — the honest evaluation strategy.

Functions:
    gene_stratified_split: DataFrame → generator of (train_idx, test_idx)
    load_clinvar_variants: filepath → DataFrame of labeled variants
    build_training_dataset: list[VariantRecord] → DataFrame feature matrix
"""
import logging
from pathlib import Path
from typing import Generator, Optional

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold

from varis.m5_scoring.feature_extractor import extract_features
from varis.models.variant_record import VariantRecord

logger = logging.getLogger(__name__)


def gene_stratified_split(
    df: pd.DataFrame,
    n_splits: int = 5,
) -> Generator[tuple[np.ndarray, np.ndarray], None, None]:
    """Gene-stratified cross-validation splits.

    Uses StratifiedGroupKFold: groups by gene, stratified on label.
    No gene appears in both train and test within a fold.

    Args:
        df: DataFrame with 'gene' and 'label' columns.
        n_splits: Number of CV folds.

    Yields:
        Tuples of (train_indices, test_indices).
    """
    sgkf = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=42)
    groups = df["gene"]
    labels = df["label"]

    for train_idx, test_idx in sgkf.split(df, labels, groups):
        yield train_idx, test_idx


def load_clinvar_variants(filepath: str) -> pd.DataFrame:
    """Parse ClinVar variant_summary.txt into labeled DataFrame.

    Filters to missense variants with Pathogenic or Benign classifications.

    Args:
        filepath: Path to ClinVar variant_summary.txt file.

    Returns:
        DataFrame with columns: gene, variant, label (1=pathogenic, 0=benign).
    """
    filepath = Path(filepath)
    if not filepath.exists():
        logger.warning(f"ClinVar file not found: {filepath}")
        return pd.DataFrame(columns=["gene", "variant", "label"])

    try:
        df = pd.read_csv(filepath, sep="\t", low_memory=False)
    except Exception as e:
        logger.warning(f"Failed to read ClinVar file: {e}")
        return pd.DataFrame(columns=["gene", "variant", "label"])

    # Filter to missense single nucleotide variants
    if "Type" in df.columns:
        df = df[df["Type"] == "single nucleotide variant"]

    # Filter to clear classifications
    if "ClinicalSignificance" in df.columns:
        path_mask = df["ClinicalSignificance"].str.contains(
            "Pathogenic", case=False, na=False,
        ) & ~df["ClinicalSignificance"].str.contains(
            "Conflicting", case=False, na=False,
        )
        benign_mask = df["ClinicalSignificance"].str.contains(
            "Benign", case=False, na=False,
        ) & ~df["ClinicalSignificance"].str.contains(
            "Conflicting", case=False, na=False,
        )
        df = df[path_mask | benign_mask].copy()
        df["label"] = path_mask[df.index].astype(int)

    # Extract gene and variant columns
    gene_col = "GeneSymbol" if "GeneSymbol" in df.columns else "gene"
    result = pd.DataFrame({
        "gene": df[gene_col] if gene_col in df.columns else None,
        "variant": df.get("Name", df.index.astype(str)),
        "label": df["label"],
    })

    return result.dropna(subset=["gene", "label"]).reset_index(drop=True)


def build_training_dataset(
    records: list[VariantRecord],
    labels: Optional[list[int]] = None,
) -> pd.DataFrame:
    """Build feature matrix from VariantRecords.

    Args:
        records: List of VariantRecords with features computed.
        labels: Optional list of labels (0=benign, 1=pathogenic).
            If None, labels are not included.

    Returns:
        DataFrame with features and optional 'label' column.
    """
    rows = []
    for record in records:
        features = extract_features(record)
        rows.append(features)

    df = pd.DataFrame(rows)

    if labels is not None:
        df["label"] = labels

    return df
