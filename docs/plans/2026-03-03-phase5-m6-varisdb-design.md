# Phase 5 Design: M6 — VarisDB Web Platform

**Date:** 2026-03-03
**Status:** Approved
**Depends on:** M1-M5 (all pipeline modules complete)
**Scope:** Layers 1-4 (API + React + SHAP waterfall + Mol* 3D viewer)

## Decisions Made

- **Layers 1-4** in scope (API, React search/detail, SHAP waterfall, Mol* viewer)
- **SQLAlchemy + DATABASE_URL** — SQLite for dev/tests, PostgreSQL for production
- **ThreadPoolExecutor** with DB-backed job tracking (no Celery/Redis)
- **Vite + React + Tailwind** frontend
- **Recharts** for SHAP waterfall
- **pdbe-molstar** for 3D structure viewer
- Consistent endpoint naming: `investigations` as the resource noun
- Idempotent investigations: return cached result unless `force=true`

## Backend: FastAPI + SQLAlchemy

### API Endpoints

| Method | Path | What it does |
|--------|------|-------------|
| `POST` | `/api/v1/investigations` | Submit new investigation (async, returns job_id). Idempotent: returns existing result if variant already investigated (unless `force=true`). |
| `GET` | `/api/v1/investigations/{variant_id}` | Read pre-computed investigation result (fast, for frontend). Returns 4-section payload. |
| `GET` | `/api/v1/jobs/{job_id}` | Poll job status (queued/running/succeeded/failed + current_step) |
| `GET` | `/api/v1/variants?q={query}` | Search variants by gene, HGVS, ClinVar ID, rsID |
| `GET` | `/api/v1/variants/stats` | Database statistics (variant count, genes, model version) |
| `GET` | `/health` | Health check |

### Database Schema (SQLAlchemy ORM)

Two tables, dialect-agnostic (no PostgreSQL-specific types):

**`variants` table:**
```
id              — Integer, primary key, autoincrement
variant_id      — String, unique, indexed (e.g., "BRCA1_p.Arg1699Trp")
gene            — String, indexed
clinvar_id      — String, nullable, indexed
classification  — String, nullable
model_version   — String, nullable, indexed
pipeline_version — String, nullable
record_json     — Text (full VariantRecord as JSON, including null_reasons)
created_at      — DateTime, indexed
updated_at      — DateTime
```

**`jobs` table:**
```
id              — Integer, primary key, autoincrement
job_id          — String, unique, indexed (UUID)
variant_id      — String, indexed
status          — String (queued/running/succeeded/failed)
current_step    — String, nullable (e.g., "M3: Structural Analysis")
error_message   — Text, nullable
created_at      — DateTime
started_at      — DateTime, nullable
completed_at    — DateTime, nullable
```

JSON stored as TEXT (works in both SQLite and PostgreSQL). Indexed columns:
`variant_id`, `gene`, `clinvar_id`, `model_version`, `created_at`.

`DATABASE_URL` config: defaults to `sqlite:///data/varis.db`, overridden by
environment variable for PostgreSQL in production.

### Background Worker

`ThreadPoolExecutor(max_workers=2)` with DB-backed job tracking:

1. `POST /investigations` creates a `jobs` row with status `queued`
2. Executor thread updates status to `running` with `current_step`
3. Runs M1→M2→M3→M4→M5 pipeline, updating `current_step` after each module
4. On success: saves VariantRecord to `variants` table, updates job to `succeeded`
5. On failure: updates job to `failed` with `error_message`
6. Job timeout: configurable (default 10 minutes), marks as `failed: timeout`
7. On server startup: scan for `running` jobs, mark as `failed: worker_restart`

**Idempotency:** If `variant_id` already exists in `variants` table and `force` is
not set, return the existing investigation immediately without running the pipeline.

### Pydantic Response Models

5-section investigation response:

```python
class StructureSection(BaseModel):
    source: str | None              # "alphafold" / "esmfold"
    pdb_path: str | None
    residue_index: int | None
    plddt_at_residue: float | None
    plddt_mean: float | None
    coordinate_mapping_confidence: str | None
    mutation_site_confidence_bucket: str | None

class FeatureItem(BaseModel):
    name: str
    value: float | str | None
    units: str | None
    evidence_tag: str | None
    available: bool

class PredictionSection(BaseModel):
    score: float | None
    classification: str | None
    confidence_lower: float | None
    confidence_upper: float | None
    model_agreement: str | None
    individual_scores: dict | None

class ExplanationItem(BaseModel):
    feature: str
    value: float | None
    shap: float

class ProvenanceSection(BaseModel):
    data_sources: list[str]         # ["ClinVar", "gnomAD", "UniProt", "AlphaFold"]
    model_version: str | None
    pipeline_version: str | None
    investigation_timestamp: str | None
    modules_completed: list[str]

class InvestigationResponse(BaseModel):
    variant_id: str
    gene: str
    hgvs: str
    structure: StructureSection
    features: list[FeatureItem]
    prediction: PredictionSection
    explanation: list[ExplanationItem]
    provenance: ProvenanceSection
```

`explanation` items come pre-sorted by |SHAP| descending from M5. Frontend renders
without re-sorting.

### Input Validation

Two-layer validation:
1. **Security layer:** length cap (256 chars), allowed characters (alphanumeric, `.`, `_`, `-`, `(`, `)`)
2. **Domain layer:** parse HGVS / recognize ClinVar ID / rsID format
3. Friendly error: "Unsupported variant format. Expected: GENE p.Xxx000Yyy"

### Security

- **CORS:** localhost only by default. `CORS_ORIGINS` env var for deployment.
- **Rate limiting:** in-memory counter on `POST /investigations` only. Clear error:
  "Rate limit exceeded; retry in 60s." Phase 5 placeholder, Redis later.

### Error Handling

Structured JSON errors via FastAPI exception handlers:
- 404: variant not found
- 422: invalid input format (with helpful message)
- 429: rate limit exceeded
- 503: database unavailable

## Frontend: Vite + React + Tailwind

### Pages

1. **Search page** (`/`) — search bar, recent investigations, database stats
2. **Investigation page** (`/variant/{variant_id}`) — main 4-panel view
3. **Job status page** (`/jobs/{job_id}`) — step-by-step progress indicator

### Investigation Page Layout

```
┌──────────────────────────┬──────────────────────────────────┐
│                          │  Reliability Strip               │
│   Mol* 3D Viewer         │  [pLDDT: 92.4] [SASA ✓] [DSSP  │
│   (pdbe-molstar)         │  ✓] [DDG ✗ tool_missing] ...    │
│                          ├──────────────────────────────────┤
│   Color by pLDDT         │  SHAP Waterfall (Recharts)       │
│   Highlight mutation     │  ┃▓▓▓▓▓▓▓ conservation  +0.24  │
│   residue                │  ┃▓▓▓▓▓ burial          +0.19  │
│                          │  ┃▓▓▓ contacts           +0.12  │
│   Confidence banners:    │  → Score: 0.91 [0.85-0.96]     │
│   low pLDDT → badge      ├──────────────────────────────────┤
│   failed mapping → red   │  Evidence Tags                   │
│                          │  • Computational support (PP3)   │
│                          │  • Rarity evidence (PM2)         │
│                          │  "Suggested computational        │
│                          │   evidence"                      │
│                          ├──────────────────────────────────┤
│                          │  Provenance                      │
│                          │  Sources: ClinVar, gnomAD, ...   │
│                          │  Model: v2026.03                 │
└──────────────────────────┴──────────────────────────────────┘
```

### React Components

| Component | Responsibility |
|-----------|---------------|
| `App.jsx` | React Router, global layout |
| `SearchPage.jsx` | Search bar (gene/HGVS input), results list, database stats |
| `InvestigationPage.jsx` | Fetches `/investigations/{id}`, renders left/right panels |
| `JobStatusPage.jsx` | Polls `/jobs/{id}`, shows M1→M5 step progress |
| `MolstarViewer.jsx` | pdbe-molstar wrapper, pLDDT coloring, confidence-aware residue highlight |
| `ReliabilityStrip.jsx` | pLDDT badge, feature availability icons with null_reason tooltips |
| `ShapWaterfall.jsx` | Recharts horizontal bar chart, pre-sorted, animated |
| `EvidencePanel.jsx` | Evidence tags with descriptions, "Suggested computational evidence" |
| `PredictionBadge.jsx` | Score circle, classification label, confidence interval, model agreement |
| `ProvenanceFooter.jsx` | Data sources, model version, timestamp |

### Mol* Rules

- Color by pLDDT (green ≥90, yellow ≥70, red <70)
- Highlight mutation residue based on `coordinate_mapping_confidence`:
  - `"exact"` / `"high"` → highlight normally (red sphere)
  - `"low"` → yellow banner: "Residue mapping uncertain"
  - `"failed"` → red banner: "Cannot locate residue in structure", NO highlight
- If `mutation_site_confidence_bucket="low"` (even with exact mapping): show
  "Low structure confidence at this site" badge

### Job Progress UI

When polling a running job, show step indicator:

```
M1 Ingestion     ✓ Complete
M2 Structure     ✓ Complete
M3 Analysis      ● Running...
M4 Conservation  ○ Pending
M5 Scoring       ○ Pending
```

No ETA, just the current step. Updates on each poll.

## Testing Strategy

### Backend Tests

- FastAPI `TestClient` for all endpoints
- SQLite in-memory (`sqlite:///:memory:`) with transaction rollback
- Mock pipeline for async job tests
- Validate Pydantic response schemas
- Edge case: job persistence on restart (scan running → mark failed)
- Idempotency: same variant returns cached result

### Frontend Tests

- React Testing Library for component rendering
- Mock fetch for API integration
- Mol*: test mount, props, confidence banners (not WebGL rendering)
- SHAP waterfall: test data rendering, not animation

## Dependencies

**Backend (existing in pyproject.toml `[platform]`):**
- fastapi, uvicorn, sqlalchemy (asyncpg not needed for SQLite)

**Frontend (new, in `varis/m6_platform/frontend/package.json`):**
- react, react-dom, react-router-dom
- tailwindcss, @tailwindcss/vite
- recharts (SHAP waterfall)
- pdbe-molstar (3D viewer)

## File Structure

```
varis/m6_platform/
├── api/
│   ├── main.py              — FastAPI app, routes, CORS, error handlers
│   ├── models.py            — Pydantic response models
│   ├── database.py          — SQLAlchemy ORM models, session management
│   ├── worker.py            — ThreadPoolExecutor, job management
│   ├── investigation.py     — Build InvestigationResponse from VariantRecord
│   └── validation.py        — Input validation (security + domain)
├── frontend/
│   ├── package.json
│   ├── vite.config.js
│   ├── index.html
│   ├── src/
│   │   ├── App.jsx
│   │   ├── main.jsx
│   │   ├── pages/
│   │   │   ├── SearchPage.jsx
│   │   │   ├── InvestigationPage.jsx
│   │   │   └── JobStatusPage.jsx
│   │   ├── components/
│   │   │   ├── MolstarViewer.jsx
│   │   │   ├── ReliabilityStrip.jsx
│   │   │   ├── ShapWaterfall.jsx
│   │   │   ├── EvidencePanel.jsx
│   │   │   ├── PredictionBadge.jsx
│   │   │   └── ProvenanceFooter.jsx
│   │   └── api/
│   │       └── client.js    — fetch wrappers for all endpoints
│   └── tailwind.config.js
└── README.md
```
