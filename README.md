# Varis — Structural Evidence for Every Variant

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://python.org)

**An AI-powered platform that doesn't just classify mutations — it investigates them.**

Varis is a structural investigator for genetic variants. While tools like AlphaMissense give clinicians a single pathogenicity score, Varis investigates *why* a mutation is damaging — orchestrating a 41-tool bioinformatics pipeline to calculate protein destabilization, map functional domains, and generate the structural evidence clinicians need to diagnose rare disease.

A [Russell Genetics](https://russellgenetics.org) product.

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

## Prerequisites & Licensing

Varis's own code is open source under the **MIT license**. However, some structural
analysis tools require separate free academic licenses that **cannot be redistributed**
with Varis:

| Tool | License Required From | Cost | Link |
|------|----------------------|------|------|
| FoldX | CRG Barcelona | Free (academic) | https://foldxsuite.crg.eu/ |
| PyRosetta | RosettaCommons | Free (academic) | https://www.rosettacommons.org/software/license-and-download |

**If these tools are not installed**, Varis's fallback architecture will skip the
affected analysis steps and continue with the remaining tools. The ML ensemble handles
missing features natively. Results will have reduced feature depth but remain valid.

Other dependencies (FreeSASA, DSSP, HMMER, BioPython, etc.) are fully open source
and installed automatically via `pip install`.

## ACMG Evidence Codes

Varis suggests ACMG evidence codes (PP3, PM1, PM5, PS3-proxy, PP2) to **support
professional review**. It does **not** replace ACMG adjudication by qualified variant
curation teams. PS3-proxy is computational structural evidence, not wet-lab functional
data, and should be weighted accordingly. All reports include this disclaimer.

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
