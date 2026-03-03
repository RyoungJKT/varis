# Phase 4: M5 ML Scoring Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Implement the ML ensemble scoring pipeline — feature extraction, training, prediction, SHAP explanations, and evidence tagging — so M6 (UI) has real scores and explanations to display.

**Architecture:** Monolithic file-per-concern: `feature_extractor.py` (one-hot, strict column alignment) → `ensemble.py` (train 3 models, calibrate, predict) → `shap_explainer.py` (per-model TreeSHAP averaged) → `evidence_mapper.py` (multi-signal evidence tags). Training uses StratifiedGroupKFold with missingness simulation. Models saved to `data/models/`.

**Tech Stack:** CatBoost, XGBoost, LightGBM (gradient-boosted trees), scikit-learn (calibration, metrics, splits), SHAP (TreeSHAP), pandas/numpy (feature matrix), httpx (ClinVar download)

---

### Task 1: Schema Update — Rename ACMG to Evidence Tags (v1.4.0)

**Files:**
- Modify: `varis/models/variant_record.py`
- Modify: `tests/conftest.py`
- Modify: `tests/test_m5_scoring.py`

**Step 1: Write failing tests**

Replace `tests/test_m5_scoring.py` entirely:

```python
"""Tests for M5: ML Scoring Engine."""
import pytest
import numpy as np
from unittest.mock import MagicMock, patch
from varis.models.variant_record import (
    VariantRecord, create_variant_record, RECORD_SCHEMA_VERSION,
)


class TestSchemaV140:
    """Verify schema v1.4.0 — evidence tags replace ACMG codes."""

    def test_schema_version(self):
        assert RECORD_SCHEMA_VERSION == "1.4.0"

    def test_evidence_tag_fields_exist(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert hasattr(record, "evidence_tags")
        assert hasattr(record, "evidence_computational_support")
        assert hasattr(record, "evidence_rarity")
        assert hasattr(record, "evidence_energetics")
        assert hasattr(record, "evidence_domain_context")

    def test_old_acmg_fields_removed(self):
        record = create_variant_record("BRCA1", "p.Arg1699Trp")
        assert not hasattr(record, "acmg_codes")
        assert not hasattr(record, "acmg_pp3")
        assert not hasattr(record, "acmg_pp2")
        assert not hasattr(record, "acmg_ps3_proxy")
        assert not hasattr(record, "acmg_pm1")
```

**Step 2: Update VariantRecord**

1. Bump `RECORD_SCHEMA_VERSION = "1.4.0"`
2. Replace the ACMG EVIDENCE section:
```python
    # =========================================================================
    # EVIDENCE TAGS — Computational evidence inspired by ACMG criteria (M5)
    # NOTE: These are computational suggestions, NOT clinical adjudications.
    # Each tag requires multiple independent signals to prevent circular reasoning.
    # =========================================================================
    evidence_tags: Optional[list[str]] = None
    evidence_computational_support: Optional[bool] = None  # PP3-like: multi-signal
    evidence_rarity: Optional[bool] = None                  # PM2-like: gnomAD rarity
    evidence_energetics: Optional[bool] = None              # DDG-based, NOT PS3
    evidence_domain_context: Optional[bool] = None          # PM1-like: critical domain
    acmg_pm5: Optional[bool] = None                         # Keep: cross-ref to known pathogenic
```

3. Update `conftest.py` `fully_populated_record`: rename `acmg_codes` → `evidence_tags`, values to `["computational_support", "rarity_evidence", "energetics_support"]`, rename individual boolean fields.

4. Fix any other test references to old ACMG field names.

**Step 3: Run tests**

Run: `pytest tests/ -v --tb=short`

**Step 4: Commit**

```
refactor(schema): rename ACMG codes to evidence tags (v1.4.0)

Evidence tags are computational suggestions, not ACMG adjudications.
Rename acmg_codes → evidence_tags, acmg_pp3 → evidence_computational_support,
acmg_pp2 → evidence_rarity, acmg_ps3_proxy → evidence_energetics,
acmg_pm1 → evidence_domain_context.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 2: Feature Extractor — One-Hot Encoding and Strict Column Alignment

**Files:**
- Replace: `varis/m5_scoring/feature_extractor.py`
- Modify: `tests/test_m5_scoring.py`

**Step 1: Write failing tests**

Add to `tests/test_m5_scoring.py`:

```python
class TestFeatureExtractor:
    """Tests for feature_extractor.py — one-hot, strict alignment."""

    def test_extract_from_full_record(self, fully_populated_record):
        """Full record extracts all features as ordered dict."""
        from varis.m5_scoring.feature_extractor import extract_features
        features = extract_features(fully_populated_record)
        assert isinstance(features, dict)
        assert len(features) > 15  # Expanded by one-hot

    def test_extract_from_partial_record(self, m1_completed_record):
        """Missing features are None, not absent."""
        from varis.m5_scoring.feature_extractor import extract_features
        features = extract_features(m1_completed_record)
        assert features["ddg_foldx"] is None
        assert features.get("ddg_available") is not None  # Should be False or None

    def test_one_hot_encoding(self, fully_populated_record):
        """Categoricals expanded to binary columns."""
        from varis.m5_scoring.feature_extractor import extract_features
        features = extract_features(fully_populated_record)
        # secondary_structure_name="helix" → one-hot columns
        assert "secondary_structure_helix" in features
        assert "secondary_structure_sheet" in features
        assert "secondary_structure_coil" in features
        assert features["secondary_structure_helix"] == 1
        assert features["secondary_structure_sheet"] == 0

    def test_in_domain_flag(self, fully_populated_record):
        """domain_name converted to binary in_domain flag."""
        from varis.m5_scoring.feature_extractor import extract_features
        features = extract_features(fully_populated_record)
        assert "in_domain" in features
        assert features["in_domain"] == 1  # domain_name is "BRCT"
        assert "domain_name" not in features  # Raw name removed

    def test_feature_extraction_allows_null(self, m1_completed_record):
        """None values pass through (real missingness)."""
        from varis.m5_scoring.feature_extractor import extract_features
        features = extract_features(m1_completed_record)
        # Structural features should be None
        assert features["solvent_accessibility_relative"] is None

    def test_extract_with_column_order(self, fully_populated_record):
        """Features extracted in deterministic order."""
        from varis.m5_scoring.feature_extractor import extract_features
        f1 = extract_features(fully_populated_record)
        f2 = extract_features(fully_populated_record)
        assert list(f1.keys()) == list(f2.keys())

    def test_build_feature_vector(self, fully_populated_record):
        """Converts feature dict to ordered numpy array for model input."""
        from varis.m5_scoring.feature_extractor import extract_features, build_feature_vector
        features = extract_features(fully_populated_record)
        columns = list(features.keys())
        vector = build_feature_vector(features, columns)
        assert len(vector) == len(columns)

    def test_build_feature_vector_fails_on_missing_column(self, fully_populated_record):
        """Raises ValueError if a required column is missing."""
        from varis.m5_scoring.feature_extractor import extract_features, build_feature_vector
        features = extract_features(fully_populated_record)
        bad_columns = list(features.keys()) + ["nonexistent_column"]
        with pytest.raises(ValueError, match="Missing feature"):
            build_feature_vector(features, bad_columns)

    def test_missingness_simulation(self, fully_populated_record):
        """Drops feature blocks, sets availability flags."""
        from varis.m5_scoring.feature_extractor import simulate_missingness
        import copy
        record = copy.deepcopy(fully_populated_record)
        # Force a specific group to be dropped
        modified = simulate_missingness(record, groups=["sasa"])
        features = modified.get_ml_features()
        assert features["sasa_available"] is False
        assert modified.solvent_accessibility_relative is None
```

**Step 2: Implement feature_extractor.py**

Key implementation:
- `extract_features(variant_record)` → `dict` of feature_name: value
  - Gets raw features via `record.get_ml_features()`
  - One-hot encodes: secondary_structure_name, burial_category, charge_change, mutation_site_confidence_bucket
  - Converts domain_name → in_domain (1/0/None)
  - Removes raw categorical columns
  - Returns ordered dict
- `build_feature_vector(features, columns)` → list
  - Takes feature dict + column order, returns values in that order
  - Raises ValueError if column name not in features dict
  - None values stay as None (float("nan") for numpy)
- `simulate_missingness(variant_record, groups=None, rate=0.15)` → VariantRecord
  - If groups is None, randomly pick from ["ddg", "sasa", "dssp", "conservation", "contacts", "domain"]
  - Set all features in that group to None
  - Set availability flag to False

Define `FEATURE_GROUPS` mapping:
```python
FEATURE_GROUPS = {
    "ddg": ["ddg_foldx", "ddg_pyrosetta"],
    "sasa": ["solvent_accessibility_relative", "burial_category"],
    "dssp": ["secondary_structure_name"],
    "conservation": ["conservation_score"],
    "contacts": ["contacts_wt", "hbonds_wt", "packing_density"],
    "domain": ["domain_name"],
}
```

**Step 3: Run tests, commit**

```
feat(m5): implement feature extractor with one-hot encoding

One-hot for categoricals. Binary in_domain flag. Strict column
alignment with fail-fast on missing columns. Missingness simulation.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 3: Ensemble — Train, Calibrate, Predict

**Files:**
- Replace: `varis/m5_scoring/ensemble.py`
- Modify: `tests/test_m5_scoring.py`

**Step 1: Write failing tests**

Add to `tests/test_m5_scoring.py`:

```python
class TestEnsemble:
    """Tests for ensemble.py — train, calibrate, predict."""

    def test_classify_thresholds(self):
        """Classification thresholds: >0.8 pathogenic, <0.2 benign."""
        from varis.m5_scoring.ensemble import _classify
        assert _classify(0.9) == "likely_pathogenic"
        assert _classify(0.5) == "uncertain"
        assert _classify(0.1) == "likely_benign"
        assert _classify(0.8) == "uncertain"  # boundary: not >0.8
        assert _classify(0.81) == "likely_pathogenic"
        assert _classify(0.2) == "uncertain"  # boundary: not <0.2
        assert _classify(0.19) == "likely_benign"

    def test_model_agreement(self):
        """Spread thresholds for agreement."""
        from varis.m5_scoring.ensemble import _compute_model_agreement
        assert _compute_model_agreement({"a": 0.8, "b": 0.85, "c": 0.82}) == "high"
        assert _compute_model_agreement({"a": 0.7, "b": 0.85, "c": 0.82}) == "moderate"
        assert _compute_model_agreement({"a": 0.5, "b": 0.85, "c": 0.82}) == "low"

    def test_train_and_predict_roundtrip(self, tmp_path):
        """Train on synthetic data, save, load, predict."""
        from varis.m5_scoring.ensemble import train_ensemble, load_ensemble, predict_from_models
        import pandas as pd
        # Small synthetic dataset: 50 samples, 5 features
        np.random.seed(42)
        n = 50
        X = pd.DataFrame({
            "feat1": np.random.randn(n),
            "feat2": np.random.randn(n),
            "feat3": np.random.randn(n),
            "feat4": np.random.randn(n),
            "feat5": np.random.randn(n),
        })
        y = pd.Series((X["feat1"] > 0).astype(int))  # Simple decision boundary
        train_ensemble(X, y, output_dir=tmp_path)
        # Verify files saved
        assert (tmp_path / "catboost_model.cbm").exists()
        assert (tmp_path / "xgboost_model.json").exists()
        assert (tmp_path / "lightgbm_model.txt").exists()
        assert (tmp_path / "calibrator.pkl").exists()
        assert (tmp_path / "feature_columns.json").exists()
        # Load and predict
        models = load_ensemble(tmp_path)
        scores = predict_from_models(models, X.iloc[0].to_dict())
        assert "score_ensemble" in scores
        assert 0.0 <= scores["score_ensemble"] <= 1.0
        assert "score_catboost" in scores
        assert "classification" in scores
        assert "model_agreement" in scores

    def test_calibration_range(self, tmp_path):
        """Calibrated score in [0, 1]."""
        from varis.m5_scoring.ensemble import train_ensemble, load_ensemble, predict_from_models
        import pandas as pd
        np.random.seed(42)
        n = 100
        X = pd.DataFrame({"f1": np.random.randn(n), "f2": np.random.randn(n)})
        y = pd.Series((X["f1"] > 0).astype(int))
        train_ensemble(X, y, output_dir=tmp_path)
        models = load_ensemble(tmp_path)
        for i in range(10):
            scores = predict_from_models(models, X.iloc[i].to_dict())
            assert 0.0 <= scores["score_ensemble"] <= 1.0

    def test_save_load_roundtrip(self, tmp_path):
        """Save and load produces same predictions."""
        from varis.m5_scoring.ensemble import train_ensemble, load_ensemble, predict_from_models
        import pandas as pd
        np.random.seed(42)
        n = 50
        X = pd.DataFrame({"f1": np.random.randn(n), "f2": np.random.randn(n)})
        y = pd.Series((X["f1"] > 0).astype(int))
        train_ensemble(X, y, output_dir=tmp_path)
        models = load_ensemble(tmp_path)
        sample = X.iloc[0].to_dict()
        s1 = predict_from_models(models, sample)
        s2 = predict_from_models(models, sample)
        assert s1["score_ensemble"] == pytest.approx(s2["score_ensemble"])
```

**Step 2: Implement ensemble.py**

Key functions:
- `train_ensemble(X, y, output_dir)` — trains 3 models with config params, fits Platt calibrator on held-out fold, saves all artifacts
- `load_ensemble(model_dir)` → dict with loaded models + calibrator + feature_columns
- `predict_from_models(models, features_dict)` → dict with score_ensemble, score_catboost, score_xgboost, score_lightgbm, classification, model_agreement, confidence_lower, confidence_upper
- `predict(variant_record, features)` → VariantRecord (convenience wrapper that loads models and calls predict_from_models)
- `_classify(score)` → str
- `_compute_model_agreement(scores)` → str
- `_train_catboost(X, y, params)`, `_train_xgboost(X, y, params)`, `_train_lightgbm(X, y, params)`

Training saves:
- Individual model files (cbm, json, txt)
- `calibrator.pkl` (sklearn CalibratedClassifierCV or LogisticRegression)
- `feature_columns.json` (strict column order)
- `training_metadata.json` (manifest)

Calibration: Use `sklearn.linear_model.LogisticRegression` fitted on raw_score → label from a held-out 20% of training data.

**Step 3: Run tests, commit**

```
feat(m5): implement ensemble training, calibration, and prediction

CatBoost + XGBoost + LightGBM with Platt calibration. Classification
thresholds, model agreement, save/load roundtrip.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 4: SHAP Explainer

**Files:**
- Replace: `varis/m5_scoring/shap_explainer.py`
- Modify: `tests/test_m5_scoring.py`

**Step 1: Write failing tests**

```python
class TestSHAPExplainer:
    """Tests for shap_explainer.py — per-variant and global explanations."""

    def test_shap_top_features_format(self, tmp_path):
        """Returns list of {feature, value, shap} dicts."""
        from varis.m5_scoring.shap_explainer import explain_from_models
        from varis.m5_scoring.ensemble import train_ensemble, load_ensemble
        import pandas as pd
        np.random.seed(42)
        n = 50
        X = pd.DataFrame({"f1": np.random.randn(n), "f2": np.random.randn(n)})
        y = pd.Series((X["f1"] > 0).astype(int))
        train_ensemble(X, y, output_dir=tmp_path)
        models = load_ensemble(tmp_path)
        top_features = explain_from_models(models, X.iloc[0].to_dict(), top_n=10)
        assert isinstance(top_features, list)
        assert len(top_features) <= 10
        for item in top_features:
            assert "feature" in item
            assert "value" in item
            assert "shap" in item

    def test_shap_sorted_by_abs_value(self, tmp_path):
        """SHAP features sorted by |shap| descending."""
        from varis.m5_scoring.shap_explainer import explain_from_models
        from varis.m5_scoring.ensemble import train_ensemble, load_ensemble
        import pandas as pd
        np.random.seed(42)
        n = 50
        X = pd.DataFrame({"f1": np.random.randn(n), "f2": np.random.randn(n)})
        y = pd.Series((X["f1"] > 0).astype(int))
        train_ensemble(X, y, output_dir=tmp_path)
        models = load_ensemble(tmp_path)
        top_features = explain_from_models(models, X.iloc[0].to_dict())
        shap_abs = [abs(f["shap"]) for f in top_features]
        assert shap_abs == sorted(shap_abs, reverse=True)

    def test_global_importance(self, tmp_path):
        """Global importance is a dict of feature→mean_abs_shap."""
        from varis.m5_scoring.shap_explainer import compute_global_importance
        from varis.m5_scoring.ensemble import train_ensemble, load_ensemble
        import pandas as pd
        np.random.seed(42)
        n = 50
        X = pd.DataFrame({"f1": np.random.randn(n), "f2": np.random.randn(n)})
        y = pd.Series((X["f1"] > 0).astype(int))
        train_ensemble(X, y, output_dir=tmp_path)
        models = load_ensemble(tmp_path)
        importance = compute_global_importance(models, X)
        assert isinstance(importance, dict)
        assert "f1" in importance
        assert all(v >= 0 for v in importance.values())
```

**Step 2: Implement shap_explainer.py**

Key functions:
- `explain_from_models(models, features_dict, top_n=10)` → list[dict]
  - Compute TreeSHAP for each base model
  - Average SHAP values across models
  - Sort by |SHAP|, return top N
- `explain_prediction(variant_record, features)` → VariantRecord
  - Loads models, calls explain_from_models, sets shap_top_features
- `compute_global_importance(models, X)` → dict
  - Mean |SHAP| across all samples per feature

TreeSHAP approach:
- CatBoost: `shap.TreeExplainer(catboost_model)`
- XGBoost: `shap.TreeExplainer(xgboost_model)`
- LightGBM: `shap.TreeExplainer(lightgbm_model)`
- Average the SHAP value arrays

**Step 3: Run tests, commit**

```
feat(m5): implement SHAP explainer with per-model averaging

Per-variant TreeSHAP averaged across 3 models. Top 10 features
sorted by |SHAP|. Global importance precomputed at training.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 5: Evidence Mapper (Replaces ACMG Mapper)

**Files:**
- Create: `varis/m5_scoring/evidence_mapper.py`
- Remove: `varis/m5_scoring/acmg_mapper.py`
- Modify: `tests/test_m5_scoring.py`

**Step 1: Write failing tests**

```python
class TestEvidenceMapper:
    """Tests for evidence_mapper.py — multi-signal evidence tags."""

    def test_computational_support_multi_signal(self, fully_populated_record):
        """Requires conservation > 0.9 AND (buried OR high score)."""
        from varis.m5_scoring.evidence_mapper import map_evidence_tags
        fully_populated_record.conservation_score = 0.95
        fully_populated_record.burial_category = "core"
        fully_populated_record.score_ensemble = 0.85
        result = map_evidence_tags(fully_populated_record)
        assert result.evidence_computational_support is True
        assert "computational_support" in result.evidence_tags

    def test_computational_support_needs_multiple_signals(self, fully_populated_record):
        """Score alone is not enough — needs conservation too."""
        from varis.m5_scoring.evidence_mapper import map_evidence_tags
        fully_populated_record.conservation_score = 0.5  # Low conservation
        fully_populated_record.score_ensemble = 0.95  # High score
        fully_populated_record.burial_category = "surface"
        result = map_evidence_tags(fully_populated_record)
        assert result.evidence_computational_support is False

    def test_rarity_evidence(self, fully_populated_record):
        """gnomAD < 0.0001 triggers rarity."""
        from varis.m5_scoring.evidence_mapper import map_evidence_tags
        fully_populated_record.gnomad_frequency = 0.00003
        result = map_evidence_tags(fully_populated_record)
        assert result.evidence_rarity is True

    def test_rarity_absent_from_gnomad(self, fully_populated_record):
        """Absent from gnomAD (None) triggers rarity."""
        from varis.m5_scoring.evidence_mapper import map_evidence_tags
        fully_populated_record.gnomad_frequency = None
        result = map_evidence_tags(fully_populated_record)
        assert result.evidence_rarity is True

    def test_energetics_support(self, fully_populated_record):
        """DDG > 2.0 triggers energetics."""
        from varis.m5_scoring.evidence_mapper import map_evidence_tags
        fully_populated_record.ddg_mean = 3.5
        result = map_evidence_tags(fully_populated_record)
        assert result.evidence_energetics is True

    def test_energetics_no_ddg(self, fully_populated_record):
        """No DDG available → energetics is False, not error."""
        from varis.m5_scoring.evidence_mapper import map_evidence_tags
        fully_populated_record.ddg_mean = None
        fully_populated_record.ddg_available = False
        result = map_evidence_tags(fully_populated_record)
        assert result.evidence_energetics is False

    def test_domain_context(self, fully_populated_record):
        """Critical domain triggers domain context."""
        from varis.m5_scoring.evidence_mapper import map_evidence_tags
        fully_populated_record.domain_criticality = "critical"
        result = map_evidence_tags(fully_populated_record)
        assert result.evidence_domain_context is True

    def test_domain_context_peripheral(self, fully_populated_record):
        """Peripheral domain does NOT trigger domain context."""
        from varis.m5_scoring.evidence_mapper import map_evidence_tags
        fully_populated_record.domain_criticality = "peripheral"
        result = map_evidence_tags(fully_populated_record)
        assert result.evidence_domain_context is False
```

**Step 2: Implement evidence_mapper.py**

```python
"""Evidence Mapper — Maps computational evidence to suggested tags.

Evidence tags are inspired by ACMG criteria but are NOT clinical
adjudications. Each tag requires multiple independent signals to
prevent circular reasoning ("the model proves itself").

Tags:
  computational_support (PP3-like): conservation > 0.9 AND (buried core OR high score)
  rarity_evidence (PM2-like): gnomAD < 0.0001 or absent
  energetics_support: DDG > 2.0 kcal/mol
  domain_context (PM1-like): critical domain
"""
```

Also remove old `acmg_mapper.py` via `git rm`.

**Step 3: Run tests, commit**

```
feat(m5): implement evidence mapper with multi-signal requirements

Replaces ACMG mapper. computational_support requires conservation AND
burial. rarity_evidence uses gnomAD threshold. energetics_support from
DDG. domain_context from criticality. Removes acmg_mapper.py.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 6: Data Loader — ClinVar Download and Gene-Stratified Split

**Files:**
- Replace: `varis/m5_scoring/data_loader.py` (if exists, or create)
- Modify: `tests/test_m5_scoring.py`

**Step 1: Write failing tests**

```python
class TestDataLoader:
    """Tests for data_loader.py — ClinVar parsing and splitting."""

    def test_gene_stratified_split(self):
        """No gene appears in both train and test."""
        from varis.m5_scoring.data_loader import gene_stratified_split
        import pandas as pd
        df = pd.DataFrame({
            "gene": ["BRCA1"] * 20 + ["TP53"] * 20 + ["CFTR"] * 20,
            "label": [1] * 10 + [0] * 10 + [1] * 10 + [0] * 10 + [1] * 10 + [0] * 10,
            "feat": range(60),
        })
        for train_idx, test_idx in gene_stratified_split(df, n_splits=3):
            train_genes = set(df.iloc[train_idx]["gene"])
            test_genes = set(df.iloc[test_idx]["gene"])
            assert train_genes.isdisjoint(test_genes), "Gene in both train and test!"

    def test_split_preserves_class_ratio(self):
        """Stratified: similar pathogenic/benign ratio in each fold."""
        from varis.m5_scoring.data_loader import gene_stratified_split
        import pandas as pd
        df = pd.DataFrame({
            "gene": ["BRCA1"] * 20 + ["TP53"] * 20 + ["CFTR"] * 20 + ["MLH1"] * 20 + ["MSH2"] * 20,
            "label": ([1] * 10 + [0] * 10) * 5,
            "feat": range(100),
        })
        ratios = []
        for train_idx, test_idx in gene_stratified_split(df, n_splits=5):
            test_labels = df.iloc[test_idx]["label"]
            ratio = test_labels.mean()
            ratios.append(ratio)
        # Each fold should have roughly 50% pathogenic (within tolerance)
        for r in ratios:
            assert 0.2 < r < 0.8, f"Extreme class imbalance in fold: {r}"
```

**Step 2: Implement data_loader.py**

Key functions:
- `gene_stratified_split(df, n_splits=5)` → generator of (train_idx, test_idx)
  - Uses `sklearn.model_selection.StratifiedGroupKFold`
  - Groups = df["gene"], stratify on df["label"]
- `load_clinvar_variants(filepath)` → pd.DataFrame
  - Parses ClinVar variant_summary.txt
  - Filters to missense, Pathogenic/Benign labels
  - Returns DataFrame with gene, variant, label columns
- `build_training_dataset(records: list[VariantRecord])` → pd.DataFrame
  - Extracts features from each VariantRecord
  - Builds feature matrix with labels

**Step 3: Run tests, commit**

```
feat(m5): implement data loader with gene-stratified split

StratifiedGroupKFold ensures no gene leaks between train/test.
ClinVar variant parser and feature matrix builder.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 7: Ablation Infrastructure

**Files:**
- Create: `varis/m5_scoring/ablation.py`
- Modify: `tests/test_m5_scoring.py`

**Step 1: Write failing tests**

```python
class TestAblation:
    """Tests for ablation.py — feature group removal."""

    def test_ablation_removes_group(self):
        """Removing a feature group sets those columns to NaN."""
        from varis.m5_scoring.ablation import drop_feature_group
        import pandas as pd
        X = pd.DataFrame({
            "conservation_score": [0.9, 0.8],
            "ddg_foldx": [1.0, 2.0],
            "sasa_available": [True, True],
            "conservation_available": [True, True],
        })
        result = drop_feature_group(X, "conservation")
        assert result["conservation_score"].isna().all()
        assert result["conservation_available"].eq(False).all()
        # Other features untouched
        assert not result["ddg_foldx"].isna().any()

    def test_ablation_runs_and_returns_metrics(self, tmp_path):
        """Run ablation produces metrics dict."""
        from varis.m5_scoring.ablation import run_ablation
        import pandas as pd
        np.random.seed(42)
        n = 60
        X = pd.DataFrame({
            "f1": np.random.randn(n),
            "f2": np.random.randn(n),
            "conservation_score": np.random.randn(n),
            "conservation_available": [True] * n,
        })
        y = pd.Series((X["f1"] > 0).astype(int))
        genes = pd.Series(["A"] * 20 + ["B"] * 20 + ["C"] * 20)
        results = run_ablation(X, y, genes, groups=["conservation"], output_dir=tmp_path)
        assert "without_conservation" in results
        assert "roc_auc" in results["without_conservation"]
```

**Step 2: Implement ablation.py**

Key functions:
- `drop_feature_group(X, group_name)` → DataFrame with group columns set to NaN/False
- `run_ablation(X, y, genes, groups, output_dir)` → dict of {group_name: metrics}
  - For each group: drop features, retrain ensemble, evaluate, report metrics

**Step 3: Run tests, commit**

```
feat(m5): implement ablation infrastructure

Drop feature groups, retrain, compare metrics. Proves each module
adds value. Includes leakage ablation for external predictors.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 8: M5 Orchestrator and Integration Test

**Files:**
- Modify: `varis/m5_scoring/__init__.py`
- Modify: `tests/test_m5_scoring.py`

**Step 1: Write failing tests**

```python
class TestM5Orchestrator:
    """Tests for M5 orchestration."""

    def test_m5_predict_with_mock_models(self, fully_populated_record, tmp_path):
        """Full M5 pipeline with trained models."""
        from varis.m5_scoring import run
        from varis.m5_scoring.ensemble import train_ensemble
        import pandas as pd
        # Train tiny models
        np.random.seed(42)
        n = 50
        X = pd.DataFrame({col: np.random.randn(n) for col in [
            "ddg_foldx", "ddg_pyrosetta", "solvent_accessibility_relative",
            "contacts_wt", "hbonds_wt", "packing_density",
            "conservation_score", "mutation_site_plddt",
            "gnomad_frequency", "alphamissense_score",
            "in_domain", "ddg_available", "sasa_available",
            "dssp_available", "conservation_available",
            "domain_available", "contacts_available",
            "secondary_structure_helix", "secondary_structure_sheet",
            "secondary_structure_coil", "burial_core", "burial_surface",
        ]})
        y = pd.Series((X["conservation_score"] > 0).astype(int))
        train_ensemble(X, y, output_dir=tmp_path)

        with patch("varis.m5_scoring.ensemble.MODELS_DIR", tmp_path):
            result = run(fully_populated_record)

        assert "M5" in result.modules_completed
        assert result.score_ensemble is not None
        assert 0.0 <= result.score_ensemble <= 1.0
        assert result.classification in ("likely_pathogenic", "uncertain", "likely_benign")
        assert result.evidence_tags is not None

    def test_m5_no_features(self, empty_record):
        """No features → M5 fails gracefully."""
        from varis.m5_scoring import run
        result = run(empty_record)
        assert "M5" in (result.modules_failed or [])
```

**Step 2: Update M5 orchestrator**

Replace `varis/m5_scoring/__init__.py` to use new module names:
- Import from `evidence_mapper` (not `acmg_mapper`)
- Import `extract_features` and `build_feature_vector` from `feature_extractor`
- Load models via `ensemble.load_ensemble()`
- Call `predict_from_models`, `explain_from_models`, `map_evidence_tags`

**Step 3: Run full test suite**

```
pytest tests/ -v --tb=short
```

**Step 4: Commit**

```
feat(m5): implement M5 orchestrator with full pipeline

Extract features → predict with ensemble → SHAP explain →
map evidence tags. Integration test with trained models.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 9: Benchmark Framework and Training Manifest

**Files:**
- Replace: `varis/m5_scoring/benchmarks.py`
- Modify: `tests/benchmark_variants.json`
- Modify: `tests/test_m5_scoring.py`

**Step 1: Expand benchmark variants**

Update `tests/benchmark_variants.json` to have at least a template with 5 variants (the full 50+ will be added as data becomes available):

```json
[
  {"gene": "BRCA1", "hgvs": "p.Arg1699Trp", "expected": "likely_pathogenic"},
  {"gene": "BRCA1", "hgvs": "p.Lys1183Arg", "expected": "likely_benign"},
  {"gene": "TP53", "hgvs": "p.Arg175His", "expected": "likely_pathogenic"},
  {"gene": "BRCA2", "hgvs": "p.Asp2723His", "expected": "likely_pathogenic"},
  {"gene": "CFTR", "hgvs": "p.Gly551Asp", "expected": "likely_pathogenic"}
]
```

**Step 2: Implement benchmarks.py**

Key functions:
- `run_benchmarks(model_dir, benchmark_file)` → dict of results per variant
- `save_training_manifest(output_dir, metrics, metadata)` → writes training_metadata.json

**Step 3: Run tests, commit**

```
feat(m5): implement benchmark framework and training manifest

Benchmark regression tests against known variants. Training
manifest with version, metrics, library versions, preprocessing.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 10: Final Cleanup and Verification

**Step 1: Search for stale references**

- `acmg_codes`, `acmg_pp3`, `acmg_pp2`, `acmg_ps3_proxy`, `acmg_pm1` (except `acmg_pm5`)
- `from varis.m5_scoring.acmg_mapper`
- Old `map_acmg_codes` references

**Step 2: Verify imports**

```
python -c "from varis.m5_scoring import run"
```

**Step 3: Run full test suite**

```
pytest tests/ -v --tb=short
```

**Step 4: Commit cleanup**

```
chore: Phase 4 cleanup — remove stale ACMG references

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```
