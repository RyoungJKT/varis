# Phase 4 Design: M5 — ML Ensemble Scoring

**Date:** 2026-03-03
**Status:** Approved
**Depends on:** M1 (ingestion), M3 (structural features), M4 (conservation)

## Decisions Made

- **ClinVar download** for training data, train on small curated subset (~500 variants, ~20 genes)
- **Config defaults** for model hyperparameters (no Optuna — deferred to M7)
- **Monolithic architecture** — one file per concern, shared feature matrix
- **Gene-stratified split** (StratifiedGroupKFold) as headline evaluation
- **Platt calibration** on separate calibration fold
- **Missingness simulation** during training (drop 10-20% of feature blocks)
- **Ablation infrastructure** to prove each module adds value
- **Evidence tags** (not ACMG codes) — honest computational evidence labeling
- **One-hot encoding** for small categoricals, binary `in_domain` flag
- **Per-model TreeSHAP** averaged, explaining pre-calibration score

## Module Architecture

### Files

| File | Responsibility |
|------|---------------|
| `data_loader.py` | Download ClinVar variants, filter missense, build feature matrix, gene-stratified split |
| `feature_extractor.py` | Extract features from VariantRecord into ordered array. One-hot categoricals. Strict alignment with `feature_columns.json`. Missingness simulation. |
| `ensemble.py` | Train CatBoost + XGBoost + LightGBM. Calibrate (Platt). Predict. Classify. Model agreement. Save/load. |
| `shap_explainer.py` | Per-model TreeSHAP averaged. Top 10 features per variant. Global importance precomputed. |
| `evidence_mapper.py` | Map computational evidence tags (replaces acmg_mapper.py). Multi-signal requirements. |
| `ablation.py` | Drop feature groups, retrain, compare metrics. Leakage ablation (without external predictors). |
| `__init__.py` | Orchestrate: extract → predict → explain → map evidence. |

### Data Flow — Training

```
ClinVar variant_summary.txt
  → filter missense with Pathogenic/Benign labels
  → run M1-M4 pipeline per variant (or load from cache)
  → feature_extractor: build feature matrix (one-hot, missingness sim)
  → StratifiedGroupKFold by gene (5 folds)
  → train 3 models per fold with config defaults
  → calibrate on held-out calibration fold (Platt scaling)
  → evaluate: ROC-AUC, PR-AUC, precision/recall, calibration ECE
  → save models + manifest to data/models/v2026.03/
```

### Data Flow — Inference

```
VariantRecord (after M1-M4)
  → feature_extractor: extract ordered features (strict column alignment)
    - Fail fast if column NAME missing (bug)
    - Allow None values (real missingness, models handle natively)
  → ensemble.predict: 3 individual scores + calibrated ensemble
  → shap_explainer: per-model TreeSHAP, average, top 10 features
  → evidence_mapper: multi-signal evidence tags
  → populate VariantRecord
```

## Feature Engineering

### Categorical Handling

One-hot encode (shared matrix across all 3 models):

| Feature | Values | Encoding |
|---------|--------|----------|
| `secondary_structure_name` | helix, sheet, coil | 3 binary columns |
| `burial_category` | core, surface | 2 binary columns |
| `charge_change` | +ve→neutral, etc. | N binary columns |
| `mutation_site_confidence_bucket` | high, medium, low | 3 binary columns |

Binary flag (replaces high-cardinality categorical):
- `domain_name` → `in_domain` (True/False). Domain criticality is more useful than the name.

### Feature Columns (strict order)

`feature_columns.json` saved at training time. Contains:
- All numeric features (ddg_foldx, solvent_accessibility_relative, etc.)
- One-hot column names (secondary_structure_helix, secondary_structure_sheet, etc.)
- Availability flag columns (ddg_available, sasa_available, etc.)

At inference, the extractor builds the vector in exactly this order.
Missing column name → `ValueError` (fail fast — this is a bug).
Missing value (None) → allowed, models handle natively.

### Missingness Simulation

During training, for 10-20% of samples:
1. Pick a random feature group (ddg, sasa, dssp, conservation, contacts, domain)
2. Set all features in that group to None
3. Set corresponding `*_available` flag to False
4. This teaches the model to degrade gracefully when tools fail

## Model Training

### Split Strategy

**StratifiedGroupKFold** (5 folds):
- Groups = gene symbols (no gene appears in both train and test)
- Stratified by label (maintains pathogenic/benign ratio per fold)
- One fold reserved for calibration (not used for evaluation)

### Models

All trained with config defaults from `varis/config.py`:
- CatBoost: 1000 iterations, lr=0.05, depth=6
- XGBoost: 1000 estimators, lr=0.05, depth=6
- LightGBM: 1000 estimators, lr=0.05, depth=6

### Calibration

Platt scaling (logistic regression on raw ensemble score) fitted on the
dedicated calibration fold. Maps raw model output to calibrated probability.
Saved as `calibrator.pkl`.

### Ensemble Score

```
raw_score = mean(catboost_pred, xgboost_pred, lightgbm_pred)
score_ensemble = calibrator.predict_proba(raw_score)
```

### Classification

```
score > 0.8  → "likely_pathogenic"
score < 0.2  → "likely_benign"
otherwise    → "uncertain"
```

Labeled as model classification, not clinical classification.

### Model Agreement

```
spread = max(individual_scores) - min(individual_scores)
spread < 0.1  → "high"
spread < 0.2  → "moderate"
otherwise     → "low"
```

## SHAP Explanations

### Per-Variant (on demand)

1. Compute TreeSHAP for each base model separately (fast, tree-native)
2. Average SHAP values across 3 models
3. SHAP explains the **pre-calibration** score
4. Return top 10 features sorted by |SHAP value|
5. Cache with VariantRecord (`shap_top_features`)

Format:
```json
[
  {"feature": "conservation_score", "value": 0.98, "shap": 0.15},
  {"feature": "solvent_accessibility_relative", "value": 0.03, "shap": 0.12},
  ...
]
```

### Global Importance (precomputed at training)

Mean |SHAP| across all training samples per feature. Saved in manifest
as `global_shap_importance.json`. Used by M6 for summary visualizations.

### Calibration Transparency

UI shows: "SHAP values explain the raw model score. The displayed probability
is calibrated separately." Calibration curve available in the manifest.

## Evidence Tags (Replaces ACMG Mapper)

Evidence tags are **computational suggestions inspired by ACMG criteria**,
not ACMG adjudications. Each tag requires multiple independent signals
to prevent circular reasoning.

| Tag | Criteria | ACMG Analog | Note |
|-----|----------|-------------|------|
| `computational_support` | High conservation (>0.9) AND (buried core OR high model score) — multiple independent signals required | PP3-like | Not just "model says pathogenic" |
| `rarity_evidence` | gnomAD frequency < 0.0001 or absent | PM2-like | Not PP2 (that's gene constraint) |
| `energetics_support` | ddg_mean > 2.0 kcal/mol (when available) | Functional-impact proxy | NOT PS3 — computational, not functional assay |
| `domain_context` | In critical domain (domain_criticality=="critical") with low benign variation in that domain | PM1-like | Requires domain-specific benign rate check |

### Schema Changes

Rename in VariantRecord:
- `acmg_codes` → `evidence_tags` (list[str])
- `acmg_pm1` → `evidence_domain_context` (bool)
- `acmg_pm5` → keep as-is (cross-ref to known pathogenic at position)
- `acmg_pp3` → `evidence_computational_support` (bool)
- `acmg_pp2` → `evidence_rarity` (bool)
- `acmg_ps3_proxy` → `evidence_energetics` (bool)

## Ablation Studies

### Feature Group Ablations

Drop each feature group, retrain, report metric change:
- Without structural features (M3): drop all SASA/DSSP/contacts/DDG
- Without conservation (M4): drop conservation_score, position_entropy, etc.
- Without external predictors: drop AlphaMissense, gnomAD frequency
- Without availability flags: drop all *_available columns

### Leakage Ablation

- Train without AlphaMissense (likely trained on ClinVar)
- Train without conservation scores from ConSurf (may overlap)
- Report: "Some external predictors may share training labels; ablations show robustness"

## Training Manifest (`training_metadata.json`)

```json
{
  "model_version": "v2026.03",
  "feature_schema_version": "1.3.0",
  "dataset": {
    "source": "clinvar",
    "export_date": "2026-03-01",
    "export_file_hash": "sha256:abc123...",
    "total_variants": 523,
    "pathogenic": 312,
    "benign": 211,
    "genes": 22
  },
  "split_strategy": "stratified_group_kfold",
  "num_folds": 5,
  "calibration_method": "platt_scaling",
  "metrics": {
    "roc_auc": {"mean": 0.87, "std": 0.03},
    "pr_auc": {"mean": 0.84, "std": 0.04},
    "calibration_ece": 0.03
  },
  "feature_columns": ["ddg_foldx", "solvent_accessibility_relative", ...],
  "preprocessing": {
    "categorical_encoding": "one_hot",
    "missingness_strategy": "native_handling",
    "missingness_simulation_rate": 0.15,
    "missingness_block_groups": ["ddg", "sasa", "dssp", "conservation", "contacts", "domain"]
  },
  "library_versions": {
    "catboost": "1.2.x",
    "xgboost": "2.0.x",
    "lightgbm": "4.0.x",
    "shap": "0.43.x",
    "scikit-learn": "1.3.x"
  },
  "ablation_results": {
    "without_structural": {"roc_auc": 0.78},
    "without_conservation": {"roc_auc": 0.82},
    "without_external": {"roc_auc": 0.85}
  },
  "git_commit": "abc123",
  "training_timestamp": "2026-03-03T12:00:00Z"
}
```

## Model Artifacts (`data/models/v2026.03/`)

```
catboost_model.cbm
xgboost_model.json
lightgbm_model.txt
calibrator.pkl
feature_columns.json
categorical_mappings.json
training_metadata.json
global_shap_importance.json
```

## Testing Strategy

- Default: offline, no real training (mock models)
- Integration tests use small synthetic datasets
- Benchmark regression tests against `tests/benchmark_variants.json`

| Test | What it verifies |
|------|-----------------|
| `test_feature_extraction_order` | Features match feature_columns.json exactly |
| `test_feature_extraction_fails_on_mismatch` | Raises ValueError if column name missing |
| `test_feature_extraction_allows_null_values` | None values pass through (not a bug) |
| `test_one_hot_encoding` | Categoricals correctly expanded |
| `test_missingness_simulation` | Drops feature blocks, sets availability flags |
| `test_ensemble_predict_shape` | Returns 3 scores + ensemble |
| `test_calibration_range` | Calibrated score in [0, 1] |
| `test_classification_thresholds` | >0.8 pathogenic, <0.2 benign |
| `test_model_agreement` | Spread thresholds correct |
| `test_shap_per_model_averaged` | SHAP from 3 models, averaged |
| `test_shap_top_features_format` | List of {feature, value, shap} dicts |
| `test_shap_global_importance` | Precomputed, sums to ~1.0 |
| `test_evidence_computational_support` | Multi-signal requirement |
| `test_evidence_rarity` | gnomAD threshold |
| `test_evidence_energetics` | DDG threshold |
| `test_evidence_domain_context` | Critical domain check |
| `test_save_load_roundtrip` | Save models, load, predict same result |
| `test_gene_stratified_split` | No gene in both train and test |
| `test_ablation_runs` | Feature group removal changes metrics |
| `test_m5_integration` | Full pipeline on fixture data |

## Dependencies

Uses existing optional ML deps from pyproject.toml:
- catboost, xgboost, lightgbm (ensemble models)
- scikit-learn (calibration, metrics, StratifiedGroupKFold)
- shap (TreeSHAP explanations)
- No new dependencies needed
