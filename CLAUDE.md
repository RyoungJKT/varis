# CLAUDE.md — Varis Project Instructions

## What This Project Is

Varis is a structural investigator for genetic variants. It determines whether a patient's
genetic mutation causes disease by investigating the structural evidence — not by classifying
with a black-box score. It is a product of Russell Genetics.

Repository: https://github.com/RussellGenetics/varis

## The One Rule That Overrides Everything

**If a module fails, set its fields to `null` in the Variant Record and continue.**
Never crash. Never raise unhandled exceptions that stop the pipeline.
The system degrades gracefully — fewer features, but always functional.

## Architecture

### The Variant Record

Every module communicates through ONE shared data structure: `VariantRecord` (defined in
`varis/models/variant_record.py`). Each module reads from it and writes its results back.

- Module writes its results to specific fields in the VariantRecord
- If a module fails, its fields remain `None`
- No module imports from or depends on another module's internals
- Each module is a standalone unit with defined inputs and outputs

### The 7 Modules

| Module | Directory | What It Does | Depends On |
|--------|-----------|-------------|------------|
| M1 | `m1_ingestion` | Parses variants, retrieves data from public databases | Nothing (entry point) |
| M2 | `m2_structure` | Obtains and prepares 3D protein structure | M1 (needs protein ID) |
| M3 | `m3_structural_analysis` | Extracts structural features from 3D structure | M2 (needs structure) |
| M4 | `m4_conservation` | Calculates evolutionary conservation | M1 (needs sequence, NOT M2/M3) |
| M5 | `m5_scoring` | Trains/runs ML pathogenicity ensemble | Any subset of M3 + M4 features |
| M6 | `m6_platform` | VarisDB: database, API, frontend, reports | M5 (needs scores to display) |
| M7 | `m7_evolution` | Auto-retrain, tool discovery, evolution log | M5 + M6 |

**Critical independence:** M3 and M4 are completely independent. If all structural analysis
fails (M3), conservation (M4) still works. If conservation fails, structural analysis still works.
The ML model (M5) takes whatever features it gets.

### Fallback Philosophy

Every critical capability has a primary tool, a fallback, and a worst case:

| Capability | Primary | Fallback 1 | Fallback 2 | Worst Case |
|-----------|---------|-----------|-----------|------------|
| 3D structure | AlphaFold DB | ESMFold | Swiss-Model | Skip structural features |
| ΔΔG stability | EvoEF2 | PyRosetta (optional) | DDGun | Skip ΔΔG, use other features |
| Secondary structure | DSSP | mdtraj DSSP | BioPython | Assign from pLDDT |
| Solvent accessibility | FreeSASA | DSSP SASA | BioPython | Estimate from depth |
| Conservation | BLAST + Clustal | PSI-BLAST + MAFFT | ConSurf API | gnomAD frequency proxy |
| Domain ID | HMMER + Pfam | InterProScan API | UniProt annotations | Skip |
| 3D visualization | Mol* | NGL Viewer | 3Dmol.js | Static images |
| Search | Elasticsearch | PostgreSQL FTS | SQLite FTS5 | Python search |
| Frontend | React + Tailwind | Next.js | Plain HTML | Streamlit |
| ML model | CatBoost+XGBoost+LightGBM | Any two models | Any single model | Logistic Regression |

## Coding Standards

### Python
- Python 3.11+
- Type hints on ALL function signatures
- Docstrings on ALL public functions (Google style)
- Every function that can fail must catch exceptions and return None/default, never crash
- Use `logging` module, never `print()` for status messages
- Import from the module's own directory, never reach into another module's internals

### Error Handling Pattern

```python
import logging

logger = logging.getLogger(__name__)

def some_tool_wrapper(variant_record: dict) -> dict:
    """Runs SomeTool and adds results to variant record.

    Args:
        variant_record: The shared VariantRecord dict.

    Returns:
        The variant_record with new fields populated (or None on failure).
    """
    try:
        # ... do the work ...
        variant_record["some_field"] = result
    except Exception as e:
        logger.warning(f"SomeTool failed for {variant_record.get('variant_id', 'unknown')}: {e}")
        variant_record["some_field"] = None
    return variant_record
```

### File Organization
- One tool wrapper per file
- Each file has: imports, logger, main function(s), helper functions
- Test file mirrors source file: `foldx_wrapper.py` → `test_m3_structural.py::test_foldx_wrapper`

### Naming Conventions
- Files: `snake_case.py`
- Classes: `PascalCase`
- Functions: `snake_case`
- Constants: `UPPER_SNAKE_CASE`
- Variant Record fields: `snake_case`

## The Variant Record Schema

The canonical schema is in `varis/models/variant_record.py`. Key field groups:

- `input.*` — What the user provided (gene, variant notation)
- `normalization.*` — Coordinate mapping, canonical transcript, position validation
- `identifiers.*` — Cross-referenced IDs (ClinVar, UniProt, etc.)
- `protein.*` — Protein sequence and metadata
- `structure.*` — 3D structure data and paths
- `structural_features.*` — FoldX, DSSP, FreeSASA, etc. results
- `feature_availability.*` — Explicit flags: ddg_available, sasa_available, etc.
- `conservation.*` — BLAST, alignment, conservation scores
- `scoring.*` — ML ensemble predictions, SHAP, confidence
- `metadata.*` — Timestamps, pipeline version, tool versions, schema version
- `null_reasons` — Dict mapping field names to NullReason codes

### Critical Rules

1. **Schema versioning**: Every record has `record_schema_version` (semver). Increment
   on any field change. Never mix schema versions in training data.

2. **Null reason codes**: NEVER set a field to None without recording why in
   `null_reasons`. Use `record.set_with_reason(field, value, reason)` instead of
   direct assignment. Valid reasons: `not_attempted`, `tool_missing`, `tool_crashed`,
   `license_unavailable`, `low_confidence_structure`, `rate_limited`, `timed_out`,
   `no_data_available`, `intentionally_skipped`, `coordinate_mapping_failed`,
   `upstream_dependency_failed`, `isoform_mismatch`, `validation_failed`.

3. **Feature availability flags**: When a structural feature is computed (or fails),
   call `record.set_feature_status("ddg", True/False, reason)`. These flags are
   included in ML features via `record.get_ml_features()` to prevent the model
   from learning artifacts from missingness.

4. **Coordinate validation**: After mapping positions between databases, ALWAYS verify
   the reference amino acid matches. A mismatch means wrong structural analysis.

## Build Phases

Building in this order. Each phase produces a working system.

1. **Phase 1 (Weeks 1-3):** M1 — CLI tool that retrieves all data from public databases
2. **Phase 2 (Weeks 3-6):** M2 + M3 — Structural analysis, easy tools first
3. **Phase 3 (Weeks 5-7):** M4 — Conservation (independent of Phase 2)
4. **Phase 4 (Weeks 6-9):** M5 — ML ensemble training and SHAP
5. **Phase 5 (Weeks 8-11):** M6 — VarisDB web platform
6. **Phase 6 (Weeks 10-14):** M7 — Self-evolution and launch

## The 3-Day Rule

Never spend more than 3 days stuck on one tool. Switch to the fallback, document
the failure, and keep moving. A working system with 8 features is infinitely better
than a broken system that was supposed to have 15.

## M5 Evaluation Rules (Phase 4)

When building the ML ensemble, follow these evaluation rules strictly:

- **Lead with gene-stratified split**: Hold out entire genes from training. If BRCA1
  is in the test set, no BRCA1 variants appear in training. This is the honest headline
  number — expect 5-10 points lower than random split.
- **Also implement**: Time-split validation (train on pre-cutoff, test on post-cutoff)
  and random stratified k-fold (for comparison only, never the headline).
- **Metrics**: ROC-AUC, PR-AUC, precision/recall at threshold, calibration plots.
- **Simulated missingness**: During training, randomly drop entire feature blocks for
  10-20% of samples. Set feature_available_* flags to False for dropped features.
  This prevents the model from learning shortcuts from systematic missingness.
- **Benchmark regression tests**: Maintain a fixed set of 50+ well-characterized
  variants in `tests/benchmark_variants.json`. Re-evaluate with every candidate model.
  Any unexpected reclassification blocks deployment.

## M7 Governance Rules (Phase 6)

When building auto-retrain (Loop 1), implement these governance safeguards:

- **Model versioning**: Tag every model with date stamp (v2026.MM). Archive previous
  versions. Current + candidate + 6 months of history always available.
- **Deploy gate**: New model deploys ONLY if ALL metrics improve AND regression tests
  pass. If any metric drops, keep current model, archive candidate as rejected, log why.
- **Rollback protocol**: Previous model version restorable immediately from archive.
  Record rollback in Evolution Log with reason.
- **Build governance BEFORE enabling auto-retrain.** The version archive and regression
  test infrastructure is what makes monthly retraining safe.

## Testing

- Every module has a test file in `tests/`
- Use `pytest` as the test runner
- The validation variant `BRCA1 p.Arg1699Trp` is the canonical test case
- Tests should work offline where possible (mock API responses)
- `tests/conftest.py` has shared fixtures including a sample VariantRecord

## Do NOT

- Import from one module into another (use the VariantRecord as the interface)
- Let any single tool failure crash the pipeline
- Use deep learning for the ML model (tree-based models with SHAP for interpretability)
- Hardcode API keys (use environment variables via .env)
- Skip type hints or docstrings
- Use `print()` instead of `logging`
- Import `requests` — use `httpx` everywhere (sync: `httpx.Client()`, async: `httpx.AsyncClient()`)

## Licensing

### Code License
Varis's own code is MIT license.

### ΔΔG Tool Licensing & Data Strategy
- **EvoEF2** (MIT license): Primary ΔΔG tool. Fully open — no restrictions on computed data.
- **PyRosetta** (RosettaCommons): Optional secondary ΔΔG tool. Free for non-commercial use,
  non-redistributable. Cannot be bundled.
- **FoldX** (CRG Barcelona): Not used. License too restrictive for open data publishing.

### VarisDB Data License
- **VarisDB computed data**: CC BY 4.0 (fully open) for all EvoEF2-computed values.
- **PyRosetta-derived values**: If included, flagged separately in the database with
  `ddg_source="pyrosetta"` so users know the provenance. These values carry the same
  CC BY 4.0 license (the non-commercial restriction is on running the software, not
  on publishing computed results — confirmed by PyRosetta license terms).
- When both EvoEF2 and PyRosetta values are available, agreement between them increases
  confidence. The `ddg_mean` field averages available methods.

### Non-Redistributable Dependencies
The README, Dockerfile, and any installation docs must clearly state that PyRosetta
requires its own license. The fallback architecture means the pipeline works without
it (EvoEF2 covers ΔΔG, reduced features if neither is available).

## Infrastructure Rules

### HTTP Client: Standardize on `httpx`

Use `httpx` for ALL outbound HTTP calls. Do not use `requests`.
- FastAPI handlers and async code: `httpx.AsyncClient()`
- CLI scripts and synchronous pipeline code: `httpx.Client()`
- This eliminates duplicate mocking patterns, inconsistent timeout/retry behavior,
  and blocking calls in async handlers

The same principle applies to all infrastructure: one library per capability.
Scientific tool fallbacks (FoldX vs PyRosetta) are intentional alternatives where
only one runs at a time. Infrastructure duplication (two HTTP clients, two PDF
generators loaded simultaneously) is just mess. If a fallback tool exists, gate it
behind a runtime check — don't import both at module level.

### API Architecture: Separate Read from Compute

Two very different endpoints — don't let them blur:
- `GET /api/v1/investigate/{variant_id}` — **reads pre-computed results** from the
  database + computes SHAP. Fast. This is what the React frontend calls.
- `POST /api/v1/variants/investigate` — **triggers the full M1-M5 pipeline**.
  Expensive (minutes). Returns a job ID immediately. Pipeline runs in a background
  worker (Celery/RQ). Never run this synchronously in a request handler.

Cache AlphaFold structure downloads — the same PDB file should not be re-downloaded
for every variant in the same protein.

### Pydantic Response Models

Define Pydantic models for all API responses. FastAPI's value is automatic validation
and OpenAPI docs — use it:
- `InvestigationResponse` — the 4-section JSON payload (structure, features, prediction, explanation)
- `VariantSummary` — lightweight record for search results
- `SearchResponse` — paginated list of VariantSummary
- `JobStatus` — for async investigation submissions

Use Pydantic for request validation too. Variant IDs must be validated (alphanumeric +
limited punctuation only) to prevent path traversal or injection.

### Error Handling in API Layer

The pipeline error pattern (catch, log, set null, continue) does not apply to the API.
API errors need structured responses:
- Use FastAPI exception handlers for consistent error JSON
- Return appropriate HTTP status codes (404 for unknown variant, 422 for invalid input,
  503 if database is down)
- Log errors as structured JSON (not print statements) for debugging

### Security Defaults

- CORS: Do NOT ship `allow_origins=["*"]` in production. Default to localhost only.
  Add a config variable `CORS_ORIGINS` for deployment.
- Rate limiting: Add rate limits on `/investigate` and `/variants/investigate`
  endpoints — these are computationally expensive.
- Input sanitization: Validate all identifiers (variant IDs, gene names) against
  strict patterns before passing to database queries.

### Testing the API Layer

- Use FastAPI's `TestClient` for endpoint tests
- Database tests: use a separate test database or SQLite in-memory, with transaction
  rollback between tests
- Mock external API clients at the client boundary (e.g., mock `ClinVarClient.fetch`,
  not the HTTP library directly) for clean test isolation
- Network client classes should accept an `httpx.Client` parameter for dependency
  injection in tests

## M6 Platform (VarisDB) — UI Architecture Rules

The VarisDB frontend must be **honest** — every visual element reflects real computed data.

### Backend to Frontend Contract

FastAPI serves `/api/v1/investigate/{variant_id}` returning four JSON sections:
- `structure`: source, chain, residue_index, plddt_at_residue, coordinate_mapping_confidence
- `features`: list of { name, value, units, evidence_tag, available }
- `prediction`: score, classification, confidence bounds, model_agreement, individual_scores
- `explanation`: list of { feature, value, shap } sorted by |shap| — from ACTUAL trained models

The `explanation` section is populated from the VariantRecord's `scoring.shap_top_features`
field, which is computed by the M5 SHAP explainer using the trained ensemble models.

**CRITICAL**: SHAP values are computed server-side by the real CatBoost/XGBoost/LightGBM
models. The React frontend is a pure visualization layer. No ML computation in the browser.

### Frontend Layout (React + Tailwind + Mol*)

1. **Left panel**: Mol* 3D viewer (pdbe-molstar). Color by pLDDT, highlight mutation residue.
   - coordinate_mapping_confidence "exact"/"high" -> highlight residue normally
   - "low" -> yellow banner: "Residue mapping uncertain"
   - "failed" -> red banner, do NOT highlight any residue
2. **Right panel**: Reliability Strip (top) + SHAP Waterfall (animated, Recharts) + ACMG Evidence Panel
3. **Reliability Strip**: pLDDT badge (green/yellow/red), feature availability indicators
   with null_reason on hover, ensemble agreement flag

### Rules
- Never render a SHAP chart from mock/placeholder data
- Never highlight a residue when coordinate_mapping_confidence is "failed"
- Always show the Reliability Strip — it prevents overclaiming
- ACMG codes are labeled "Suggested evidence" not "Assigned evidence"
- The investigation endpoint reads directly from the VariantRecord — same schema,
  same null_reasons, same feature_available flags
