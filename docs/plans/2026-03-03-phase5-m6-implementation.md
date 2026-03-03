# Phase 5: M6 VarisDB Platform Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Build the VarisDB web platform — FastAPI backend with SQLAlchemy database, React frontend with Tailwind CSS, SHAP waterfall chart (Recharts), and Mol* 3D protein viewer — so users can search, submit, and visually explore variant investigations.

**Architecture:** Backend: FastAPI with Pydantic response models, SQLAlchemy ORM (SQLite dev / PostgreSQL production), ThreadPoolExecutor for async pipeline jobs. Frontend: Vite + React + Tailwind, Recharts for SHAP waterfall, pdbe-molstar for 3D structure. Backend serves the built frontend in production.

**Tech Stack:** Python: FastAPI, SQLAlchemy, uvicorn, Pydantic. JavaScript: React 18, Vite, Tailwind CSS 4, Recharts, pdbe-molstar. Testing: pytest + FastAPI TestClient (backend), Vitest + React Testing Library (frontend).

---

### Task 1: Pydantic Response Models

**Files:**
- Create: `varis/m6_platform/api/models.py`
- Create: `tests/test_m6_api.py`

**Step 1: Write failing tests**

Create `tests/test_m6_api.py`:

```python
"""Tests for M6: VarisDB API."""
import pytest
from pydantic import ValidationError


class TestPydanticModels:
    """Tests for API response models."""

    def test_investigation_response_valid(self):
        from varis.m6_platform.api.models import InvestigationResponse
        data = {
            "variant_id": "BRCA1_p.Arg1699Trp",
            "gene": "BRCA1",
            "hgvs": "p.Arg1699Trp",
            "structure": {"source": "alphafold"},
            "features": [],
            "prediction": {"score": 0.91, "classification": "likely_pathogenic"},
            "explanation": [],
            "provenance": {
                "data_sources": ["ClinVar"],
                "modules_completed": ["M1"],
            },
        }
        resp = InvestigationResponse(**data)
        assert resp.variant_id == "BRCA1_p.Arg1699Trp"

    def test_feature_item_model(self):
        from varis.m6_platform.api.models import FeatureItem
        item = FeatureItem(name="conservation_score", value=0.95, units=None, evidence_tag=None, available=True)
        assert item.available is True

    def test_explanation_item_model(self):
        from varis.m6_platform.api.models import ExplanationItem
        item = ExplanationItem(feature="conservation_score", value=0.95, shap=0.24)
        assert item.shap == 0.24

    def test_job_status_model(self):
        from varis.m6_platform.api.models import JobStatusResponse
        job = JobStatusResponse(job_id="abc-123", status="running", variant_id="BRCA1_p.Arg1699Trp", current_step="M3: Structural Analysis")
        assert job.status == "running"
```

**Step 2: Implement models.py**

Create `varis/m6_platform/api/models.py` with all Pydantic models from the design doc: `StructureSection`, `FeatureItem`, `PredictionSection`, `ExplanationItem`, `ProvenanceSection`, `InvestigationResponse`, `VariantSummary`, `SearchResponse`, `JobStatusResponse`, `InvestigationRequest`, `StatsResponse`.

**Step 3: Run tests, commit**

```
feat(m6): implement Pydantic response models

StructureSection, FeatureItem, PredictionSection, ExplanationItem,
ProvenanceSection, InvestigationResponse, JobStatusResponse.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 2: SQLAlchemy Database Layer

**Files:**
- Replace: `varis/m6_platform/api/database.py`
- Modify: `tests/test_m6_api.py`

**Step 1: Write failing tests**

Add to `tests/test_m6_api.py`:

```python
class TestDatabase:
    """Tests for database.py — SQLAlchemy ORM."""

    @pytest.fixture
    def db_session(self, tmp_path):
        """Create an in-memory SQLite session for testing."""
        from varis.m6_platform.api.database import init_db, get_session
        engine, SessionLocal = init_db(f"sqlite:///{tmp_path}/test.db")
        session = SessionLocal()
        yield session
        session.close()

    def test_save_and_get_variant(self, db_session, fully_populated_record):
        """Save a VariantRecord and retrieve it."""
        from varis.m6_platform.api.database import save_variant_record, get_variant_record
        variant_id = save_variant_record(db_session, fully_populated_record)
        assert variant_id is not None
        retrieved = get_variant_record(db_session, variant_id)
        assert retrieved is not None
        assert retrieved["gene_symbol"] == "BRCA1"

    def test_get_nonexistent_variant(self, db_session):
        """Unknown variant_id returns None."""
        from varis.m6_platform.api.database import get_variant_record
        result = get_variant_record(db_session, "NONEXISTENT")
        assert result is None

    def test_search_by_gene(self, db_session, fully_populated_record):
        """Search finds variant by gene name."""
        from varis.m6_platform.api.database import save_variant_record, search_variants
        save_variant_record(db_session, fully_populated_record)
        results = search_variants(db_session, "BRCA1")
        assert len(results) >= 1

    def test_create_and_get_job(self, db_session):
        """Create a job and retrieve its status."""
        from varis.m6_platform.api.database import create_job, get_job, update_job_status
        job_id = create_job(db_session, "BRCA1_p.Arg1699Trp")
        assert job_id is not None
        job = get_job(db_session, job_id)
        assert job["status"] == "queued"
        update_job_status(db_session, job_id, "running", current_step="M1: Ingestion")
        job = get_job(db_session, job_id)
        assert job["status"] == "running"
        assert job["current_step"] == "M1: Ingestion"

    def test_idempotent_variant_exists(self, db_session, fully_populated_record):
        """Check if variant already exists."""
        from varis.m6_platform.api.database import save_variant_record, variant_exists
        save_variant_record(db_session, fully_populated_record)
        assert variant_exists(db_session, fully_populated_record.variant_id) is True
        assert variant_exists(db_session, "NONEXISTENT") is False
```

**Step 2: Implement database.py**

Replace `varis/m6_platform/api/database.py` with SQLAlchemy ORM:
- `VariantRow` and `JobRow` ORM models (dialect-agnostic, no PostgreSQL types)
- `init_db(database_url)` → (engine, SessionLocal)
- `get_session()` → session (for FastAPI dependency injection)
- `save_variant_record(session, variant_record)` → variant_id
- `get_variant_record(session, variant_id)` → dict | None
- `search_variants(session, query, limit=20)` → list[dict]
- `variant_exists(session, variant_id)` → bool
- `create_job(session, variant_id)` → job_id (UUID)
- `get_job(session, job_id)` → dict | None
- `update_job_status(session, job_id, status, current_step=None, error_message=None)`
- `mark_stale_jobs_failed(session)` — on startup, mark running → failed

`DATABASE_URL` defaults to `sqlite:///data/varis.db`, read from env var.

**Step 3: Run tests, commit**

```
feat(m6): implement SQLAlchemy database layer

Dialect-agnostic ORM for variants and jobs tables. SQLite dev,
PostgreSQL production via DATABASE_URL. Idempotency check.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 3: Input Validation

**Files:**
- Create: `varis/m6_platform/api/validation.py`
- Modify: `tests/test_m6_api.py`

**Step 1: Write failing tests**

```python
class TestValidation:
    """Tests for input validation."""

    def test_valid_hgvs(self):
        from varis.m6_platform.api.validation import validate_variant_input
        result = validate_variant_input("BRCA1", "p.Arg1699Trp")
        assert result["valid"] is True

    def test_security_length_cap(self):
        from varis.m6_platform.api.validation import validate_variant_input
        result = validate_variant_input("A" * 300, "p.Arg1699Trp")
        assert result["valid"] is False
        assert "length" in result["error"].lower()

    def test_security_bad_chars(self):
        from varis.m6_platform.api.validation import validate_variant_input
        result = validate_variant_input("BRCA1; DROP TABLE", "p.Arg1699Trp")
        assert result["valid"] is False

    def test_friendly_error_message(self):
        from varis.m6_platform.api.validation import validate_variant_input
        result = validate_variant_input("BRCA1", "invalid_format")
        assert result["valid"] is False
        assert "Expected" in result["error"]
```

**Step 2: Implement validation.py**

Two-layer: security (length + charset) then domain (HGVS pattern recognition).

**Step 3: Run tests, commit**

```
feat(m6): implement input validation with security and domain layers

Length cap, charset filter, HGVS format recognition. Friendly errors.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 4: Investigation Response Builder

**Files:**
- Create: `varis/m6_platform/api/investigation.py`
- Modify: `tests/test_m6_api.py`

**Step 1: Write failing tests**

```python
class TestInvestigationBuilder:
    """Tests for building InvestigationResponse from VariantRecord."""

    def test_build_from_full_record(self, fully_populated_record):
        """Full record produces complete response."""
        from varis.m6_platform.api.investigation import build_investigation_response
        resp = build_investigation_response(fully_populated_record)
        assert resp.variant_id is not None
        assert resp.structure.source == "alphafold"
        assert resp.prediction.score is not None
        assert len(resp.provenance.data_sources) > 0
        assert len(resp.provenance.modules_completed) > 0

    def test_build_features_list(self, fully_populated_record):
        """Features list includes availability flags."""
        from varis.m6_platform.api.investigation import build_investigation_response
        resp = build_investigation_response(fully_populated_record)
        assert len(resp.features) > 0
        for f in resp.features:
            assert hasattr(f, "available")

    def test_build_from_partial_record(self, m1_completed_record):
        """Partial record (M1 only) still produces valid response."""
        from varis.m6_platform.api.investigation import build_investigation_response
        resp = build_investigation_response(m1_completed_record)
        assert resp.variant_id is not None
        assert resp.prediction.score is None  # M5 not run yet
        assert len(resp.explanation) == 0  # No SHAP without M5

    def test_explanation_pre_sorted(self, fully_populated_record):
        """SHAP items come pre-sorted by |shap| descending."""
        from varis.m6_platform.api.investigation import build_investigation_response
        fully_populated_record.shap_top_features = [
            {"feature": "a", "value": 1.0, "shap": 0.1},
            {"feature": "b", "value": 2.0, "shap": -0.3},
            {"feature": "c", "value": 3.0, "shap": 0.2},
        ]
        resp = build_investigation_response(fully_populated_record)
        shap_abs = [abs(e.shap) for e in resp.explanation]
        assert shap_abs == sorted(shap_abs, reverse=True)
```

**Step 2: Implement investigation.py**

`build_investigation_response(variant_record) → InvestigationResponse` — maps VariantRecord fields to the 5-section Pydantic response.

**Step 3: Run tests, commit**

```
feat(m6): implement investigation response builder

Maps VariantRecord to 5-section API response: structure, features,
prediction, explanation, provenance.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 5: Background Worker

**Files:**
- Create: `varis/m6_platform/api/worker.py`
- Modify: `tests/test_m6_api.py`

**Step 1: Write failing tests**

```python
from unittest.mock import patch, MagicMock
import time

class TestWorker:
    """Tests for background worker."""

    def test_submit_job(self, tmp_path):
        """Submit creates a queued job."""
        from varis.m6_platform.api.database import init_db, create_job, get_job
        from varis.m6_platform.api.worker import InvestigationWorker
        engine, SessionLocal = init_db(f"sqlite:///{tmp_path}/test.db")
        session = SessionLocal()
        worker = InvestigationWorker(SessionLocal, max_workers=1)
        job_id = worker.submit("BRCA1", "p.Arg1699Trp", session)
        assert job_id is not None
        job = get_job(session, job_id)
        assert job["status"] in ("queued", "running")
        worker.shutdown()
        session.close()

    def test_job_completes(self, tmp_path):
        """Mocked pipeline → job succeeds."""
        from varis.m6_platform.api.database import init_db, get_job
        from varis.m6_platform.api.worker import InvestigationWorker
        engine, SessionLocal = init_db(f"sqlite:///{tmp_path}/test.db")
        session = SessionLocal()
        with patch("varis.m6_platform.api.worker.run_pipeline") as mock_pipeline:
            mock_record = MagicMock()
            mock_record.variant_id = "BRCA1_p.Arg1699Trp"
            mock_record.to_json.return_value = "{}"
            mock_record.gene_symbol = "BRCA1"
            mock_record.clinvar_id = None
            mock_record.classification = "likely_pathogenic"
            mock_record.ensemble_version = "v2026.03"
            mock_record.pipeline_version = "v1.0"
            mock_pipeline.return_value = mock_record
            worker = InvestigationWorker(SessionLocal, max_workers=1)
            job_id = worker.submit("BRCA1", "p.Arg1699Trp", session)
            time.sleep(2)  # Wait for worker thread
            job = get_job(SessionLocal(), job_id)
            assert job["status"] == "succeeded"
            worker.shutdown()
        session.close()
```

**Step 2: Implement worker.py**

`InvestigationWorker` class:
- `__init__(SessionLocal, max_workers=2)` — creates ThreadPoolExecutor
- `submit(gene, hgvs, session)` → job_id — creates job row, submits to executor
- `_run_investigation(job_id, gene, hgvs)` — the worker thread: updates status as it runs M1→M5, saves result
- `shutdown()` — graceful shutdown
- `mark_stale_jobs(session)` — on startup, mark running → failed

**Step 3: Run tests, commit**

```
feat(m6): implement background worker with DB-backed job tracking

ThreadPoolExecutor with max_workers=2. Updates job status through
M1→M5 pipeline steps. Marks stale jobs on startup.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 6: FastAPI Application and Routes

**Files:**
- Replace: `varis/m6_platform/api/main.py`
- Modify: `tests/test_m6_api.py`

**Step 1: Write failing tests**

```python
from fastapi.testclient import TestClient

class TestAPIEndpoints:
    """Tests for FastAPI routes."""

    @pytest.fixture
    def client(self, tmp_path):
        from varis.m6_platform.api.main import create_app
        app = create_app(database_url=f"sqlite:///{tmp_path}/test.db")
        return TestClient(app)

    def test_health_check(self, client):
        resp = client.get("/health")
        assert resp.status_code == 200
        assert resp.json()["status"] == "ok"

    def test_get_investigation_not_found(self, client):
        resp = client.get("/api/v1/investigations/NONEXISTENT")
        assert resp.status_code == 404

    def test_submit_investigation(self, client):
        resp = client.post("/api/v1/investigations", json={"gene": "BRCA1", "hgvs": "p.Arg1699Trp"})
        assert resp.status_code in (200, 201, 202)
        data = resp.json()
        assert "job_id" in data or "variant_id" in data

    def test_submit_invalid_input(self, client):
        resp = client.post("/api/v1/investigations", json={"gene": "", "hgvs": ""})
        assert resp.status_code == 422

    def test_search_empty(self, client):
        resp = client.get("/api/v1/variants?q=BRCA1")
        assert resp.status_code == 200
        assert "results" in resp.json()

    def test_stats(self, client):
        resp = client.get("/api/v1/variants/stats")
        assert resp.status_code == 200

    def test_get_job_not_found(self, client):
        resp = client.get("/api/v1/jobs/nonexistent-id")
        assert resp.status_code == 404

    def test_cors_headers(self, client):
        resp = client.options("/health", headers={"Origin": "http://localhost:5173"})
        # Should include CORS headers for allowed origin
        assert resp.status_code in (200, 204)
```

**Step 2: Implement main.py**

Replace `varis/m6_platform/api/main.py` with full FastAPI app:
- `create_app(database_url=None)` — configurable for testing
- All routes with proper Pydantic request/response models
- CORS middleware (localhost:5173 default, configurable via env var)
- Exception handlers for 404, 422, 429, 503
- Rate limiting on POST /investigations (in-memory counter)
- Lifespan event to init DB and mark stale jobs on startup

**Step 3: Run tests, commit**

```
feat(m6): implement FastAPI application with all routes

GET/POST investigations, jobs, search, stats, health. CORS,
input validation, structured error responses, rate limiting.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 7: Frontend Scaffold (Vite + React + Tailwind)

**Files:**
- Create: `varis/m6_platform/frontend/` (entire scaffold)

**Step 1: Initialize Vite project**

```bash
cd varis/m6_platform/frontend
npm create vite@latest . -- --template react
npm install
npm install react-router-dom recharts
npm install -D tailwindcss @tailwindcss/vite
```

**Step 2: Configure Tailwind and Vite**

Set up `vite.config.js` with Tailwind plugin and API proxy to FastAPI (port 8000).
Set up `tailwind.config.js` and import Tailwind in main CSS.

**Step 3: Create base App.jsx with routing**

```jsx
// src/App.jsx
import { BrowserRouter, Routes, Route } from "react-router-dom";
import SearchPage from "./pages/SearchPage";
import InvestigationPage from "./pages/InvestigationPage";
import JobStatusPage from "./pages/JobStatusPage";

export default function App() {
  return (
    <BrowserRouter>
      <div className="min-h-screen bg-gray-50">
        <header className="bg-white shadow-sm border-b">
          <div className="max-w-7xl mx-auto px-4 py-3">
            <h1 className="text-xl font-bold text-gray-900">
              Varis<span className="text-blue-600">DB</span>
            </h1>
          </div>
        </header>
        <Routes>
          <Route path="/" element={<SearchPage />} />
          <Route path="/variant/:variantId" element={<InvestigationPage />} />
          <Route path="/jobs/:jobId" element={<JobStatusPage />} />
        </Routes>
      </div>
    </BrowserRouter>
  );
}
```

**Step 4: Create API client**

```javascript
// src/api/client.js
const API_BASE = "/api/v1";

export async function getInvestigation(variantId) { ... }
export async function submitInvestigation(gene, hgvs) { ... }
export async function getJobStatus(jobId) { ... }
export async function searchVariants(query) { ... }
export async function getStats() { ... }
```

**Step 5: Create placeholder pages** (SearchPage, InvestigationPage, JobStatusPage)

Each page has minimal content ("Coming soon" or basic layout) — components filled in subsequent tasks.

**Step 6: Verify `npm run dev` works**

**Step 7: Commit**

```
feat(m6): scaffold React frontend with Vite + Tailwind

Router with 3 pages, API client, Tailwind setup.
Proxy to FastAPI backend on port 8000.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 8: Search Page

**Files:**
- Create: `varis/m6_platform/frontend/src/pages/SearchPage.jsx`

**Step 1: Implement SearchPage**

- Search bar (gene + HGVS input)
- Submit button that calls `POST /investigations`
- If variant already exists: redirect to `/variant/{id}`
- If new: redirect to `/jobs/{job_id}`
- Recent investigations list (from GET /variants?q=)
- Database stats panel (from GET /stats)
- Tailwind styling: clean, minimal, centered layout

**Step 2: Verify in browser**

**Step 3: Commit**

```
feat(m6): implement search page with variant submission

Search bar, recent investigations, stats panel. Submits to
POST /investigations, redirects to job or variant page.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 9: Job Status Page

**Files:**
- Create: `varis/m6_platform/frontend/src/pages/JobStatusPage.jsx`

**Step 1: Implement JobStatusPage**

- Polls `GET /jobs/{id}` every 3 seconds
- Shows M1→M5 step indicator with status icons (✓ complete, ● running, ○ pending)
- On succeeded: redirect to `/variant/{variant_id}`
- On failed: show error message
- Tailwind styling: centered card with step timeline

**Step 2: Verify in browser**

**Step 3: Commit**

```
feat(m6): implement job status page with step progress

Polls job status, shows M1→M5 progress. Redirects on completion.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 10: Investigation Page — Layout and PredictionBadge

**Files:**
- Create: `varis/m6_platform/frontend/src/pages/InvestigationPage.jsx`
- Create: `varis/m6_platform/frontend/src/components/PredictionBadge.jsx`
- Create: `varis/m6_platform/frontend/src/components/ProvenanceFooter.jsx`

**Step 1: Implement InvestigationPage layout**

- Fetches `GET /investigations/{variant_id}`
- Two-column layout: left (Mol* placeholder), right (panels stacked vertically)
- Header: gene, HGVS, ClinVar ID
- PredictionBadge: score circle, classification label, confidence interval, model agreement
- ProvenanceFooter: data sources, model version, timestamp

**Step 2: Verify in browser**

**Step 3: Commit**

```
feat(m6): implement investigation page layout with prediction badge

Two-column layout. Score circle, classification, confidence interval,
model agreement badge, provenance footer.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 11: Reliability Strip and Evidence Panel

**Files:**
- Create: `varis/m6_platform/frontend/src/components/ReliabilityStrip.jsx`
- Create: `varis/m6_platform/frontend/src/components/EvidencePanel.jsx`

**Step 1: Implement ReliabilityStrip**

- pLDDT badge (green ≥90, yellow ≥70, red <70)
- Feature availability icons: ✓ (available) / ✗ (unavailable) with tooltip showing null_reason
- Compact horizontal strip

**Step 2: Implement EvidencePanel**

- List of evidence tags with descriptions
- Header: "Suggested computational evidence"
- Each tag shows ACMG analog and criteria met
- Styled as cards/chips

**Step 3: Verify in browser**

**Step 4: Commit**

```
feat(m6): implement reliability strip and evidence panel

pLDDT badge, feature availability icons with tooltips.
Evidence tags labeled as computational suggestions.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 12: SHAP Waterfall Chart

**Files:**
- Create: `varis/m6_platform/frontend/src/components/ShapWaterfall.jsx`

**Step 1: Implement ShapWaterfall**

- Recharts horizontal `BarChart`
- Each bar: feature name (left), SHAP value bar (colored: positive=red, negative=blue)
- Pre-sorted by |SHAP| descending (data comes pre-sorted from API)
- Final score annotation with confidence interval
- Animated bars on mount
- Tooltip showing feature value and SHAP contribution

**Step 2: Verify in browser with test data**

**Step 3: Commit**

```
feat(m6): implement SHAP waterfall chart with Recharts

Horizontal bar chart sorted by |SHAP|. Positive=red, negative=blue.
Final score with confidence interval. Animated.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 13: Mol* 3D Viewer

**Files:**
- Create: `varis/m6_platform/frontend/src/components/MolstarViewer.jsx`

**Step 1: Install pdbe-molstar**

```bash
npm install pdbe-molstar
```

Or use the web component via CDN script tag in `index.html`.

**Step 2: Implement MolstarViewer**

- Wrapper around pdbe-molstar web component
- Props: `pdbPath`, `residueIndex`, `coordinateConfidence`, `plddtBucket`
- Color by pLDDT (green/yellow/red)
- Highlight mutation residue:
  - confidence "exact"/"high" → red sphere highlight
  - confidence "low" → yellow banner: "Residue mapping uncertain"
  - confidence "failed" → red banner: "Cannot locate residue", NO highlight
- pLDDT bucket "low" → "Low structure confidence" badge
- Loading state while structure loads

**Step 3: Verify in browser with BRCA1 AlphaFold structure**

**Step 4: Commit**

```
feat(m6): implement Mol* 3D viewer with confidence-aware highlighting

pdbe-molstar wrapper. pLDDT coloring. Residue highlight gated by
coordinate_mapping_confidence. Warning banners for low/failed mapping.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 14: Wire Everything Together + Config Updates

**Files:**
- Modify: `varis/config.py`
- Modify: `varis/m6_platform/frontend/src/pages/InvestigationPage.jsx`

**Step 1: Add platform config to config.py**

```python
DATABASE_URL = os.getenv("DATABASE_URL", "sqlite:///data/varis.db")
CORS_ORIGINS = os.getenv("CORS_ORIGINS", "http://localhost:5173").split(",")
MAX_INVESTIGATION_WORKERS = int(os.getenv("MAX_INVESTIGATION_WORKERS", "2"))
JOB_TIMEOUT_SECONDS = int(os.getenv("JOB_TIMEOUT_SECONDS", "600"))
RATE_LIMIT_PER_MINUTE = int(os.getenv("RATE_LIMIT_PER_MINUTE", "10"))
```

**Step 2: Wire InvestigationPage to all components**

Connect MolstarViewer, ReliabilityStrip, ShapWaterfall, EvidencePanel, PredictionBadge, ProvenanceFooter together on the InvestigationPage.

**Step 3: Add npm scripts**

Update `package.json`:
```json
{
  "scripts": {
    "dev": "vite",
    "build": "vite build",
    "preview": "vite preview"
  }
}
```

**Step 4: Verify full flow in browser**

Start backend: `uvicorn varis.m6_platform.api.main:app --reload --port 8000`
Start frontend: `cd varis/m6_platform/frontend && npm run dev`

**Step 5: Commit**

```
feat(m6): wire all components together, add platform config

DATABASE_URL, CORS_ORIGINS, worker config. Full investigation
page with all panels connected.

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```

---

### Task 15: Final Verification and Cleanup

**Step 1: Run backend tests**

```
pytest tests/test_m6_api.py -v
pytest tests/ -v --tb=short
```

**Step 2: Verify frontend builds**

```
cd varis/m6_platform/frontend && npm run build
```

**Step 3: Search for stale references**

- Old endpoint paths (`/api/v1/investigate/`, `/api/v1/variants/investigate`)
- Old `acmg_*` references in M6 code
- Old `map_acmg_codes` imports

**Step 4: Commit cleanup**

```
chore: Phase 5 cleanup — verify build and remove stale references

Co-Authored-By: Claude Opus 4.6 (1M context) <noreply@anthropic.com>
```
