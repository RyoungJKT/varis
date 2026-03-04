# Varis — Structural Evidence for Every Variant

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://python.org)
[![CI](https://github.com/RussellGenetics/varis/actions/workflows/ci.yml/badge.svg)](https://github.com/RussellGenetics/varis/actions)

**An AI-powered platform that doesn't just classify mutations — it investigates them.**

Varis is a structural investigator for genetic variants. While tools like AlphaMissense give clinicians a single pathogenicity score, Varis investigates *why* a mutation is damaging — orchestrating a multi-tool bioinformatics pipeline to calculate protein destabilization, map functional domains, and generate the structural evidence clinicians need to diagnose rare disease.

A [Russell Genetics](https://russellgenetics.org) product.

## Project Status (v2026.03)

| Component | Status | Details |
|-----------|--------|---------|
| M1: Data Ingestion | Complete | ClinVar, UniProt, gnomAD, AlphaMissense clients |
| M2: Structure Engine | Complete | AlphaFold + ESMFold fallback |
| M3: Structural Analysis | Complete | EvoEF2 DDG, FreeSASA, DSSP, contacts, InterPro domains |
| M4: Conservation | Complete | UniProt orthologs + Clustal Omega + ConSurf fallback |
| M5: ML Scoring | Complete | Ensemble trained on 535 ClinVar variants (ROC-AUC 0.846) |
| M6: Platform | Complete | FastAPI + PostgreSQL + React frontend + PDF reports |
| M7: Self-Evolution | Partial | Auto-retrain loop done; tool scout & auto-integrator stubbed |
| CI/CD | Complete | GitHub Actions (pytest, ruff, mypy), Docker, docker-compose |
| Tests | 269 passing | Across all modules |

## Quick Start

```bash
# Install
pip install -e ".[all]"

# Investigate a variant
python -m varis BRCA1 p.Arg1699Trp

# Run validation suite
python -m varis --validate
```

## Architecture

Varis is built as 7 independent modules connected by a shared Variant Record:

| Module | What It Does |
|--------|-------------|
| M1: Data Ingestion | Parses variants, retrieves data from ClinVar/gnomAD/UniProt/AlphaFold |
| M2: Structure Engine | Obtains and prepares 3D protein structures |
| M3: Structural Analysis | Extracts structural features (ΔΔG, SASA, DSSP, contacts, domains) |
| M4: Conservation | Calculates evolutionary conservation (independent of M2/M3) |
| M5: ML Scoring | Trains/runs interpretable ensemble (CatBoost + XGBoost + LightGBM) |
| M6: Platform (VarisDB) | Database, API, frontend, 3D viewer, PDF reports |
| M7: Self-Evolution | Auto-retrain, tool discovery, evolution log |

**Key design principle:** If any module fails, the pipeline continues with whatever data it has. The ML ensemble handles missing features natively.

### Structural Features Computed

For each variant, the pipeline extracts up to 16 features:

| Feature | Source | Description |
|---------|--------|-------------|
| ddg_evoef2 | EvoEF2 | Stability change (kcal/mol), positive = destabilizing |
| ddg_foldx | FoldX (optional) | Stability change from FoldX |
| ddg_pyrosetta | PyRosetta (optional) | Stability change from PyRosetta |
| ddg_mean | M3 orchestrator | Mean of available DDG methods |
| solvent_accessibility_relative | FreeSASA | Relative SASA (0=buried, 1=exposed) |
| burial_category | FreeSASA | core/surface classification |
| secondary_structure | DSSP | helix/sheet/coil at mutation site |
| contacts_wt | BioPython | Heavy-atom contacts within 4.5A |
| hbonds_wt | BioPython | Hydrogen bonds at mutation site |
| packing_density | BioPython | Local packing (contacts / target atoms) |
| domain_criticality | InterPro | critical/important/peripheral |
| plddt_at_residue | AlphaFold/ESMFold | Structure confidence at mutation site |
| plddt_mean | AlphaFold/ESMFold | Mean structure confidence |
| gnomad_frequency | gnomAD | Population allele frequency |
| conservation_score | Clustal/ConSurf | Evolutionary conservation (0-1) |
| alphamissense_score | AlphaMissense | External pathogenicity score |

## ML Model

### Training Data

- **535 variants** from ClinVar (271 pathogenic, 264 benign)
- **19 genes**: BRCA1, BRCA2, TP53, CFTR, MSH2, MLH1, MSH6, PMS2, PTEN, RB1, APC, VHL, MEN1, RET, CDH1, PALB2, ATM, CHEK2, RAD51C
- RAD51D skipped (insufficient pathogenic variants)
- VUS explicitly excluded — they are the prediction target
- 5 benchmark variants held out, never trained on
- Simulated 15% feature missingness during training for robustness

### Performance (Gene-Stratified Cross-Validation)

| Fold | Held-Out Genes | ROC-AUC | PR-AUC |
|------|---------------|---------|--------|
| 1 | CFTR, PALB2, RB1, TP53 | 0.884 | 0.853 |
| 2 | BRCA1, BRCA2, MSH2, RAD51C | 0.827 | 0.816 |
| 3 | CDH1, MEN1, MSH6, RET | 0.860 | 0.857 |
| 4 | APC, MLH1, PTEN, VHL | 0.860 | 0.867 |
| 5 | ATM, CHEK2, PMS2 | 0.799 | 0.770 |
| **Mean** | | **0.846 +/- 0.030** | **0.833 +/- 0.036** |

These are honest numbers — entire genes are held out from training to prevent leakage.

### Ensemble Architecture

Three gradient-boosted tree models, each trained independently:
- **CatBoost** — handles categoricals natively, ordered boosting
- **XGBoost** — strong regularization, histogram-based
- **LightGBM** — leaf-wise growth, fast training

Final score = mean of three model scores. SHAP values computed per-prediction using TreeSHAP for full interpretability.

## Prerequisites & Licensing

Varis's own code is open source under the **MIT license**. However, some structural
analysis tools require separate free academic licenses that **cannot be redistributed**
with Varis:

| Tool | License | Status | Required? |
|------|---------|--------|-----------|
| EvoEF2 | MIT | Primary DDG tool, compiled from source | Recommended |
| FoldX | CRG Barcelona (free academic) | Optional fallback DDG | No |
| PyRosetta | RosettaCommons (free academic) | Optional fallback DDG | No |

**If optional tools are not installed**, Varis's fallback architecture will skip the
affected analysis steps and continue with the remaining tools. The ML ensemble handles
missing features natively. Results will have reduced feature depth but remain valid.

Other dependencies (FreeSASA, DSSP, HMMER, BioPython, etc.) are fully open source
and installed automatically via `pip install`.

### Installing EvoEF2

EvoEF2 is the primary stability prediction tool (MIT licensed). To compile from source:

```bash
git clone https://github.com/tommyhuangthu/EvoEF2.git /tmp/EvoEF2
cd /tmp/EvoEF2 && g++ -O3 -ffast-math -o EvoEF2 src/*.cpp
cp EvoEF2 ~/bin/
cp -r library ~/bin/library    # Required: parameter files must be next to binary
export EVOEF2_BINARY=~/bin/EvoEF2
```

## ACMG Evidence Codes

Varis suggests ACMG evidence codes (PP3, PM1, PM5, PS3-proxy, PP2) to **support
professional review**. It does **not** replace ACMG adjudication by qualified variant
curation teams. PS3-proxy is computational structural evidence, not wet-lab functional
data, and should be weighted accordingly. All reports include this disclaimer.

## Deployment

Varis uses a split deployment architecture: the React frontend is served from a global CDN, and the FastAPI backend runs as a separate service with its own database.

### Frontend — Vercel (free)

The React + Tailwind + Mol* frontend at `varis/m6_platform/frontend/` deploys as a static site.

```bash
# Build the frontend
cd varis/m6_platform/frontend
npm install && npm run build

# Deploy to Vercel
npx vercel --prod
```

Set the environment variable in Vercel:
```
VITE_API_URL=https://api.varis.russellgenetics.org
```

This serves the app on a global CDN with automatic HTTPS. Free tier is sufficient.

### Backend — Railway (~$7–25/mo)

The FastAPI API + PostgreSQL database deploy via Docker.

```bash
# Deploy with Railway CLI
railway login
railway init
railway up
```

Or use the included `docker-compose.yml` for any Docker-compatible host:
```bash
docker-compose up -d
```

This starts three services:
- **api** — FastAPI on port 8000 (with EvoEF2 compiled in the image)
- **db** — PostgreSQL 16 on port 5432
- **frontend** — Dev server on port 5173 (local dev only; production uses Vercel)

Environment variables for the backend:
```bash
DATABASE_URL=postgresql+asyncpg://varis:varis@db:5432/varis
CORS_ORIGINS=https://varis.russellgenetics.org
NCBI_API_KEY=           # Optional, increases rate limit
CLINVAR_API_KEY=        # Required for ClinVar submission
EVOEF2_BINARY=EvoEF2   # Compiled in Docker image
```

### Weekly Auto-Retrain

M7 retrains the ML ensemble weekly against fresh ClinVar data. Set up a cron job or GitHub Action on the backend server:

```bash
# Manual trigger
python -m varis.m7_evolution retrain

# Or via cron (every Monday at 3am UTC)
0 3 * * 1 cd /app && python -m varis.m7_evolution retrain >> /app/data/logs/retrain.log 2>&1
```

Data storage on the backend:
```
data/
├── variant_summary.txt.gz    # ClinVar dump (~1GB, cached 7 days)
├── structures/               # PDB files (~5MB each)
├── conservation/             # MSA cache per protein
├── training/variants/        # Cached VariantRecord JSONs
├── models/                   # Trained ensemble (catboost/xgboost/lightgbm .pkl)
└── varis.db                  # SQLite (dev) — PostgreSQL in production
```

### Architecture Summary

```
Vercel (free)                          Railway ($7-25/mo)
┌──────────────────────┐               ┌──────────────────────────────┐
│  React + Tailwind    │               │  FastAPI                     │
│  Mol* 3D viewer      │──── API ────▶ │  ├── /api/v1/investigate/*  │
│  SHAP waterfall      │    calls      │  ├── /api/v1/variants/*     │
│  Reliability strip   │               │  └── /api/v1/jobs/*         │
└──────────────────────┘               │                              │
                                       │  PostgreSQL                  │
                                       │  EvoEF2 binary               │
                                       │  ML models (data/models/)    │
                                       │  Structure cache             │
                                       └──────────────────────────────┘
```

## Open Science

Everything is open: the pipeline, the models, the database, the training data, and the evolution log.

- **VarisDB**: Free, searchable database at varisdb.russellgenetics.org
- **ML Model**: Weights on HuggingFace
- **Dataset**: 200K+ variants × 15 structural features on Zenodo
- **ClinVar**: Automated evidence submissions
- **Notebooks**: 5 educational Google Colab tutorials

## License

MIT — because the science should be free. Note: Some structural analysis
dependencies (FoldX, PyRosetta) require separate free academic licenses.
See [Prerequisites & Licensing](#prerequisites--licensing) above.
