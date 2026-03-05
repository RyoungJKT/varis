# M7 Evolution — Loops 2 & 3 Design

## Overview

Complete M7 Evolution by implementing Loop 2 (Tool Scout) and Loop 3 (Auto-Integrator). Both follow the Varis fallback philosophy: never crash, graceful degradation, log everything.

## Loop 2: Tool Scout

**Purpose:** Periodically scan public sources for new bioinformatics tools relevant to the pipeline.

### Sources (all free, no auth)

| Source | Method | What we look for |
|--------|--------|------------------|
| PyPI | XMLRPC/JSON API | Packages tagged with protein stability, variant effect, structural biology |
| GitHub | Search API (unauthenticated) | Repos by topic/keyword related to pipeline capabilities |
| bioRxiv | Content API | Recent papers with software tools for DDG, conservation, structure analysis |

### Flow

1. **`scan_sources()`** — HTTP queries (httpx) to each source, returns list of candidate dicts: `{name, source, url, description, version, keywords}`
2. **`score_candidate(candidate)`** — Rule-based relevance scoring: keyword overlap with pipeline capabilities (DDG, conservation, SASA, contacts, secondary structure, solvent accessibility), recency bonus, popularity bonus (GitHub stars, PyPI downloads)
3. **`deduplicate(candidates, log_db)`** — Skip candidates already logged as TOOL_DISCOVERY events in evolution log
4. **`run_scout_loop()`** — Orchestrate: scan → score → deduplicate → log high-scoring candidates as TOOL_DISCOVERY events with score and metadata

### Scoring rules

- Base score from keyword overlap (each matching keyword = +1 point)
- Pipeline-relevant keywords: `ddg, delta-g, stability, foldx, rosetta, evoef, conservation, sasa, solvent, dssp, secondary structure, variant effect, missense, pathogenicity, protein structure, alphafold, plddt, contacts, packing`
- Recency bonus: +2 if published/updated within last 90 days
- Popularity bonus: +1 if >100 GitHub stars or >1000 PyPI monthly downloads
- Threshold: score >= 3 to be logged as TOOL_DISCOVERY

### What Scout does NOT do

- No automatic installation (that's Loop 3)
- No LLM calls (offline/rule-based only)
- No wrapper generation

## Loop 3: Auto-Integrator

**Purpose:** Given a tool proposal (from scout or manual), attempt to install, run on benchmarks, and evaluate whether it adds value to the ML ensemble.

### Flow

1. **`attempt_install(package_name)`** — `pip install --dry-run` to check compatibility, then actual install in subprocess with timeout
2. **`probe_tool(package_name)`** — Import package, inspect its API (look for common function signatures like `predict`, `score`, `compute`, `run`), return `tool_info` dict with callable name + expected I/O
3. **`benchmark_new_feature(tool_info, benchmark_variants)`** — Run tool on benchmark variants via subprocess, compute a new feature column, retrain ensemble with extra feature, compare metrics to current production
4. **`attempt_integration(proposal)`** — Orchestrate: install → probe → benchmark → decide (INTEGRATE/REJECT) → log TOOL_INTEGRATION event

### Safety gates

- Dry-run install first (check for dependency conflicts)
- Subprocess execution with 60s timeout per variant
- Reject if any benchmark metric (ROC-AUC, PR-AUC) drops vs current production
- Never auto-modify pipeline source code — logs recommendation only
- Human decides whether to wire the tool into the pipeline

### Decision criteria

- INTEGRATE recommendation: ROC-AUC improves >= 0.005 AND no PR-AUC regression
- REJECT: any metric drops or tool fails on >20% of benchmark variants

## Shared additions

### CLI subcommands

```bash
python -m varis.m7_evolution scout              # Run tool discovery scan
python -m varis.m7_evolution integrate --tool X  # Attempt integration of tool X
```

### Evolution log events (already defined)

- `TOOL_DISCOVERY` — from Loop 2, logged with score + candidate metadata
- `TOOL_INTEGRATION` — from Loop 3, logged with decision + metrics delta

### Error handling

Both loops follow the Varis rule: catch all exceptions, log them, continue. A failed HTTP request to PyPI doesn't crash the scout — it logs a warning and moves to the next source. A failed tool install doesn't crash the integrator — it logs REJECT with reason.

## Testing strategy

- Mock all HTTP calls (httpx) — no network in tests
- Mock pip install subprocess calls
- Mock tool imports and function calls
- Test scoring with known candidates (high-relevance, low-relevance, duplicate)
- Test integration decision gate with mock metrics
- Test CLI subcommands parse correctly
