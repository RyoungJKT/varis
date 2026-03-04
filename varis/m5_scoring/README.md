# M5 Scoring — ML Ensemble for Variant Classification

See CLAUDE.md for architecture overview and coding standards.

## Current Model: v2026.03

Trained on 2026-03-04. Deployed to `data/models/`.

### Performance (Gene-Stratified 5-Fold CV)

| Fold | Held-Out Genes | ROC-AUC | PR-AUC |
|------|---------------|---------|--------|
| 1 | CFTR, PALB2, RB1, TP53 | 0.884 | 0.853 |
| 2 | BRCA1, BRCA2, MSH2, RAD51C | 0.827 | 0.816 |
| 3 | CDH1, MEN1, MSH6, RET | 0.860 | 0.857 |
| 4 | APC, MLH1, PTEN, VHL | 0.860 | 0.867 |
| 5 | ATM, CHEK2, PMS2 | 0.799 | 0.770 |
| **Mean** | | **0.846 +/- 0.030** | **0.833 +/- 0.036** |

Entire genes are held out from training in each fold — this is the honest headline number.

## Training Data Selection

The model trains on **confidently labeled variants only** — never on VUS:

- **Pathogenic / Likely pathogenic** from ClinVar (excluding conflicting interpretations)
- **Benign / Likely benign** from ClinVar (excluding conflicting interpretations)
- **VUS are explicitly excluded** — they are what we are trying to predict
- Minimum 1-star ClinVar review status required (filters out unreviewed submissions)

### Actual Training Set (v2026.03)

- **535 variants** selected from ClinVar (271 pathogenic, 264 benign)
- **19 genes** included: BRCA1, BRCA2, TP53, CFTR, MSH2, MLH1, MSH6, PMS2, PTEN, RB1, APC, VHL, MEN1, RET, CDH1, PALB2, ATM, CHEK2, RAD51C
- **1 gene skipped**: RAD51D (insufficient pathogenic variants in ClinVar)
- Up to **15 per class per gene** (15 pathogenic + 15 benign = 30 per gene max)
- Minimum 5 per class per gene required
- 5 benchmark variants (BRCA1 p.Arg1699Trp, etc.) held out — never trained on

### How Often (M7 Auto-Retrain)

- **Weekly cadence** — model versions use format `v{YYYY}.{MM}` (e.g., `v2026.03`)
- Downloads fresh ClinVar variant_summary.txt.gz (cached for 7 days)
- Recomputes features via M1–M4 pipeline for any new variants
- Trains candidate model with gene-stratified cross-validation
- Evaluates candidate against deploy gates (no metric regressions, no benchmark classification flips)
- Only deploys if metrics improve and no regressions are detected

### Key Design Insight

VUS are the **output**, not the input. The model learns the boundary between known pathogenic and known benign variants using structural, conservation, and population features. It then classifies uncertain variants based on where they fall relative to that boundary.

Simulated missingness during training (15% of samples get random feature groups dropped) teaches the model to handle real-world pipeline failures without learning artifacts from systematic missingness.

## Pipeline

Three-phase training pipeline (see `train.py`):

1. **Phase A (select):** Download ClinVar, select balanced missense variants per gene
2. **Phase B (compute):** Run M1–M4 on each variant, cache VariantRecords (resumable)
3. **Phase C (train):** Build feature matrix, gene-stratified CV, train final ensemble

```bash
python -m varis.m5_scoring select          # Phase A only
python -m varis.m5_scoring compute         # Phase B only (resumable)
python -m varis.m5_scoring train-only      # Phase C only (from cached features)
python -m varis.m5_scoring train           # All phases
python -m varis.m5_scoring status          # Show progress
```

### Phase B Details

For each of the 535 variants, Phase B runs the full M1–M4 pipeline:
- **M1**: ClinVar lookup, UniProt protein fetch, gnomAD frequencies, AlphaMissense score
- **M2**: AlphaFold structure download (ESMFold fallback), pLDDT extraction
- **M3**: EvoEF2 DDG, FreeSASA, DSSP, contact map, H-bonds, InterPro domains
- **M4**: UniProt orthologs, Clustal Omega MSA, conservation scoring (ConSurf fallback)

Results are cached as JSON VariantRecords. Phase B is resumable — if interrupted, it picks up from the last completed variant.

## Ensemble

Three gradient-boosted tree models, each trained independently:

- **CatBoost** — handles categoricals natively, ordered boosting
- **XGBoost** — strong regularization, histogram-based
- **LightGBM** — leaf-wise growth, fast training

Final score = mean of three model scores, calibrated via logistic regression.
SHAP values computed per-prediction using TreeSHAP for interpretability.

## Feature Groups

The ensemble uses up to 16 features organized in groups:

| Group | Features | Source |
|-------|----------|--------|
| DDG (stability) | ddg_evoef2, ddg_foldx, ddg_pyrosetta, ddg_mean | M3 |
| Structure | solvent_accessibility_relative, secondary_structure, contacts_wt, hbonds_wt, packing_density | M3 |
| Domain | domain_criticality | M3 (InterPro) |
| Confidence | plddt_at_residue, plddt_mean | M2 |
| Population | gnomad_frequency | M1 |
| Conservation | conservation_score | M4 |
| External | alphamissense_score | M1 |

Ablation testing (`ablation.py`) measures per-group importance by dropping each group and retraining.
